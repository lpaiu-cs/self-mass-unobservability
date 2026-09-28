"""Close the finite momentum complement with a certified 16-point Gauss rule."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
from pathlib import Path
import gzip,json,math,shlex,subprocess,sys
import numpy as np
import mpmath as mp
from mpmath import iv
import sympy as sp
import gr_response_momentum_tail as tail
import gr_plasma_interval as native

cusp=tail.cusp;rule=cusp.rule;ROOT=tail.ROOT;OUT=tail.OUT.parent/'gr-smooth-response'
CACHE=native.g.CACHE/'smooth-response';I=cusp.I;low=cusp.low;high=cusp.high;sha=rule.sha
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    tail.verify();assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    data=dict(np.load(cusp.OUT/'inputs.npz'));lines=[]
    for row in tail.records():
        i,j=row['position'],row['z_index'];lo,hi=str(data['root_intervals'][i])[1:-1].split(',')
        fields=[i*12+j,row['cell'],j,*[float(x).hex() for x in [data['beta'][i],data['Q'][i,j],data['Sref'][i]]],
            lo.strip(),hi.strip(),float(F(row['window_half_width'])).hex(),row['pcut_hex']]
        lines.append(' '.join(map(str,fields))+'\n')
    (OUT/'inputs.tsv').write_text(''.join(lines))
    paths=[ROOT/'verification/gr_smooth_response_certificate.py',ROOT/'verification/gr_smooth_response_certificate.cpp',
        tail.OUT/'manifest.json',tail.OUT/'bounds.jsonl.gz',cusp.OUT/'manifest.json',cusp.OUT/'inputs.npz',OUT/'inputs.tsv']
    save('plan.json',dict(classification='Proven',checkpoint='c212c6e',bits=128,rule_bits=256,order=16,cases=len(lines),
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in paths},processes=3,block_size=256,
        interval_headers=json.loads((native.OUT/'plan.json').read_text())['interval_headers'],
        finite_component_absolute_error_budget='1e-11',remainder_budget='1e-13',total_response_budget='3e-11',
        parameter_order=cusp.PARAMETERS,radius_half_width_ratio=4,maximum_depth=48,bound_exponent_floor=-128,
        control_positions=[0,1603,3205],control_z_indices=[0,5,8,11],control_digits=70,
        exact_rule='Shifted Legendre polynomials and rational moments prove orthogonality. Exact companion traces of Christoffel weights prove degrees 0..31 and the degree-32 norm. Isolate roots at 256 bits and enclose positive weights at 256 bits.',
        domain='The exact binary pcut from the certified tail, minus the already certified cusp interval [Q*(1-h),Q*(1+h)]. Binary parameters and dyadic subdivisions are checked for exact point arithmetic. No domain gaps or repeated cusp contribution.',
        complex_bound='At center c and radius R=4h, require R<abs(c-Q), R<c+Q, g2=1+max(c-R,0)^2-R^2>0, phase=(c+R)*R/(beta*sqrt(g2))<=1. Then Re gamma>=sqrt(g2), Re t>=(g2-1)/(beta*(sqrt(g2)+1)), abs(t)<=(c+R)^2/(beta*(sqrt(g2)+1)).',
        kernel_bound='With g=sqrt(g2), pmax=c+R: abs(A)<=pmax^2/(g*Sref), abs(B)<=pmax*(abs(1+c^2-Q^2)+R*(2c+R))/(2Q*g*Sref). Bound the two analytic logs separately by abs(log(c+Q))-log(1-R/(c+Q)) and abs(log(abs(c-Q)))-log(1-R/abs(c-Q)). Use the same six complex Fermi derivative envelopes as the cusp.',
        remainder='Cauchy gives M32<=32!*M/(R-h)^32. The 16-point Gauss remainder is <=2h*norm16*(2/3)^32*M. Bisect until analytic constraints and all local error budgets pass. Allocate 1e-13 times panel_width/pcut, so the sum is at most 1e-13. Whole-panel omission uses its full signed analytic envelope with the same budget.',
        complete_response='Sum finite enclosures, previously certified cusp enclosures and signed tail bounds. Audit the binary midpoint-to-exported-endpoint errors as exact rationals. This certifies continuous momentum at the stated finite beta,Q samples and certified eta intervals; not continuous Q or full physical EOS.',
        physical_EOS_certified=False,outer_Q_certified=False))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for rel,digest in p['interval_headers'].items():assert sha(Path(rel))==digest,rel
    return p


def gaussian_rule():
    plan=bindings();iv.prec=plan['rule_bits'];x=sp.symbols('x');n=plan['order']
    polys=[sp.Poly(sp.legendre(k,2*x-1),x) for k in range(n+1)];norms=[]
    for k,p in enumerate(polys):
        for j in range(k):assert sum(p.nth(l)/sp.Rational(l+j+1) for l in range(k+1))==0
        square=p*p;norm=sum(square.nth(l)/sp.Rational(l+1) for l in range(2*k+1));assert norm==sp.Rational(1,2*k+1);norms.append(norm)
    final=polys[-1].monic();kernel=sum(p.as_expr()**2/v for p,v in zip(polys[:-1],norms[:-1]));wp=sp.Poly(sp.invert(kernel,final.as_expr(),x),x)
    C=sp.zeros(n)
    for k in range(1,n):C[k,k-1]=1
    for k in range(n):C[k,n-1]=-final.nth(k)
    W=sum((wp.nth(k)*C**k for k in range(n)),sp.zeros(n));norm=sp.Rational(1,(2*n+1)*sp.binomial(2*n,n)**2)
    errors=[sp.factor(sp.Rational(1,k+1)-sp.trace(W*C**k)) for k in range(2*n+1)]
    assert all(e==0 for e in errors[:-1]) and errors[-1]==norm
    isolated=sp.intervals(final,eps=sp.Rational(1,2**256));assert len(isolated)==n;nodes=[];rows=[]
    for (a,b),multiplicity in isolated:
        assert multiplicity==1 and 0<a<=b<1
        node=iv.mpf([I(a).a,I(b).b]);weight=1/sum(rule.evaluate(p,node)**2/I(v) for p,v in zip(polys[:-1],norms[:-1]));assert low(weight)>0
        nodes.append((node,weight));rows.append(dict(isolation=[str(a),str(b)],node=cusp.interval_text(node),weight=cusp.interval_text(weight)))
    # Independent analytic control; exp(8x) makes the degree-32 error resolvable.
    partial=sum(F(8**k,math.factorial(k)*(k+1)) for k in range(161));first=F(8**161,math.factorial(161)*162);upper=partial+first/(1-F(8,162))
    ref=iv.mpf([I(partial).a,I(upper).b]);quad=sum(w*iv.exp(8*v) for v,w in nodes);error=ref-quad;bound=I(norm)*8**32/math.factorial(32)
    assert low(error)>=low(bound) and high(error)<=high(bound*iv.exp(8))
    save('rule.json',dict(classification='Proven',passed=True,order=n,monic_polynomial=str(final.as_expr()),norm_exact=str(norm),
        degree_0_to_32_errors_exact=list(map(str,errors)),nodes=rows,exponential_control=dict(passed=True,error=cusp.interval_text(error)),
        theorem='Degree-31 Hermite interpolation at the 16 simple roots and repeated Rolle gives the degree-32 remainder norm16*M32/32!. Exact rational moment identities and positive interval weights are independently checked.'))
    encoded=[]
    for row in rows:encoded.append(' '.join([*row['node'][1:-1].split(','),*row['weight'][1:-1].split(',')]))
    encoded.append(' '.join(cusp.interval_text(I(norm))[1:-1].split(',')));(OUT/'rule.tsv').write_text('\n'.join(encoded)+'\n')
    save('rule-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in [OUT/'rule.json',OUT/'rule.tsv']}))
    print('PASS exact 16-point Gauss rule and independent exponential control',flush=True)


def build():
    bindings();flags=shlex.split(subprocess.check_output([str(native.CAPD/'build-request15-mp/bin/capd-config'),'--cflags','--libs'],text=True))
    cmd=['g++',str(ROOT/'verification/gr_smooth_response_certificate.cpp'),f'-I{native.DEPS}/include',f'-I{native.DEPS}/include/x86_64-linux-gnu',f'-L{native.DEPS}/lib/x86_64-linux-gnu',*flags,'-O2','-o',str(CACHE/'smooth')]
    done=subprocess.run(cmd,capture_output=True,text=True);(OUT/'build.log').write_text(done.stdout+done.stderr);save('build.json',dict(command=cmd,returncode=done.returncode));assert done.returncode==0,done.stderr
    linked=subprocess.check_output(['ldd',str(CACHE/'smooth')],text=True);(OUT/'linked-libraries.txt').write_text(linked)
    libs={word:sha(Path(word)) for line in linked.splitlines() for word in line.split() if word.startswith('/') and Path(word).is_file()}
    save('runtime.json',dict(binary=str(CACHE/'smooth'),binary_sha256=sha(CACHE/'smooth'),libraries=libs))


def runtime():
    bindings()
    r=json.loads((OUT/'runtime.json').read_text());assert sha(Path(r['binary']))==r['binary_sha256']
    for name,digest in r['libraries'].items():assert sha(Path(name))==digest,name
    for rel,digest in json.loads((OUT/'rule-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    return r['binary']


def block(first,count):
    path=OUT/f'block-{first:05d}.jsonl';assert not path.exists() and not path.with_suffix('.jsonl.gz').exists()
    done=subprocess.run([str(CACHE/'smooth'),str(OUT/'inputs.tsv'),str(OUT/'rule.tsv'),str(path),str(first),str(count)],capture_output=True,text=True)
    if done.returncode:(OUT/f'block-{first:05d}-error.log').write_text(done.stdout+done.stderr)
    assert done.returncode==0,done.stderr
    with gzip.open(path.with_suffix('.jsonl.gz'),'xb',compresslevel=1) as stream:stream.write(path.read_bytes())
    path.unlink();return first


def records():
    for path in sorted(OUT.glob('block-*.jsonl.gz')):
        with gzip.open(path,'rt',encoding='utf-8') as stream:
            for line in stream:yield json.loads(line)


def run():
    plan=bindings();runtime();finished=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        work=[pool.submit(block,i,min(plan['block_size'],plan['cases']-i)) for i in range(0,plan['cases'],plan['block_size'])]
        for future in as_completed(work):
            finished.append(future.result());save('progress.json',dict(completed_blocks=sorted(finished)));print('FINITE MOMENTUM BLOCKS',len(finished),flush=True)
    aggregate()


def aggregate():
    plan=bindings();iv.prec=128;maximum=F(0);finite_max=F(0);panels=0;omitted=0;count=0
    with gzip.open(OUT/'complete-response.jsonl.gz','xt',encoding='utf-8',compresslevel=1) as stream:
        for row,window,remainder in zip(records(),cusp.records(),tail.records(),strict=True):
            i,j=divmod(count,12);assert row['case']==count and (window['position'],window['z_index'])==(i,j)==(remainder['position'],remainder['z_index'])
            assert row['cell']==window['cell']==remainder['cell'];total=[];finite_errors=[];errors=[];approx=[]
            for pair,old,bound in zip(row['enclosures'],window['enclosures'],remainder['upper_bounds'],strict=True):
                lo,hi=map(F,pair);p=F.from_float(float((lo+hi)/2));finite_errors.append(max(abs(p-lo),abs(hi-p)))
                value=cusp.interval('['+','.join(pair)+']')+cusp.interval(old)+cusp.symmetric(F(bound));text=cusp.interval_text(value);a,b=cusp.endpoints(text);point=float((a+b)/2)
                total.append(text);approx.append(point);errors.append(max(abs(a-F.from_float(point)),abs(b-F.from_float(point))))
            score=max(errors);fs=max(finite_errors);maximum=max(maximum,score);finite_max=max(finite_max,fs);panels+=row['panels'];omitted+=row['bounded_omissions'];count+=1
            assert max(map(F,row['remainder_upper']))<=F(plan['remainder_budget'])
            stream.write(json.dumps(dict(position=i,cell=row['cell'],z_index=j,enclosures=total,approximation=approx,
                maximum_error_upper=str(score),maximum_finite_error_upper=str(fs),passed=score<F(plan['total_response_budget']) and fs<F(plan['finite_component_absolute_error_budget'])),separators=(',',':'))+'\n')
    assert count==plan['cases']
    save('result.json',dict(classification='Proven',passed=maximum<F(plan['total_response_budget']) and finite_max<F(plan['finite_component_absolute_error_budget']),cases=count,
        maximum_error_upper=str(maximum),maximum_finite_error_upper=str(finite_max),panels=panels,bounded_omissions=omitted,
        display_only=dict(maximum_error=float(maximum),maximum_finite_error=float(finite_max)),continuous_momentum_certified=True,continuous_Q_certified=False,physical_EOS_certified=False))


def controls():
    plan=bindings();mp.mp.dps=plan['control_digits'];data=dict(np.load(cusp.OUT/'inputs.npz'));selected={(i,j) for i in plan['control_positions'] for j in plan['control_z_indices']}
    chosen={divmod(r['case'],12):r for r in records() if divmod(r['case'],12) in selected};checks=[]
    def M(x):
        q=F(x);return mp.mpf(q.numerator)/q.denominator
    for end in tail.records():
        i,j=end['position'],end['z_index']
        if (i,j) not in chosen:continue
        row=chosen[(i,j)];beta,Q,Sref=map(lambda x:mp.mpf(float(x)),[data['beta'][i],data['Q'][i,j],data['Sref'][i]])
        eta=sum(map(M,cusp.endpoints(str(data['root_intervals'][i]))))/2;h=M(end['window_half_width']);pcut=mp.mpf(end['pcut']);domains=[(mp.mpf(0),min(pcut,Q*(1-h)))]
        if pcut>Q*(1+h):domains.append((Q*(1+h),pcut))
        values=[]
        for k in range(6):
            def integrand(p):
                if p==0:return mp.mpf(0)
                gamma=mp.sqrt(1+p*p);t=p*p/(beta*(gamma+1));q=1/(1+mp.exp(t-eta));v=1-q
                kernel=(p*p/gamma+p*(1+p*p-Q*Q)/(2*Q*gamma)*mp.log(abs((p+Q)/(p-Q))))/Sref
                return kernel*[q,q*v,q*v*(1-2*q),t*q*v,t*q*v*(1-2*q),(t*t*(1-2*q)-t)*q*v][k]
            value=mp.mpf(0)
            for a,b in domains:
                cuts=sorted(set([a,b]+[mp.sqrt(beta)*v for v in map(mp.mpf,['.1','1','3','8','15']) if a<mp.sqrt(beta)*v<b]));value+=mp.quad(integrand,cuts)
            values.append(value)
        passed=all(M(a)<=v<=M(b) for (a,b),v in zip(row['enclosures'],values));checks.append(dict(cell=row['cell'],z_index=j,passed=passed,reference=list(map(str,values))));print('SMOOTH INDEPENDENT CONTROL',row['cell'],j,passed,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=all(c['passed'] for c in checks),controls=checks));assert len(checks)==12 and all(c['passed'] for c in checks)


def finalize():
    runtime()
    for name in ['rule.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    plan=bindings();runtime()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for name in ['rule.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    count=0;maximum=F(0)
    with gzip.open(OUT/'complete-response.jsonl.gz','rt',encoding='utf-8') as stream:
        for line in stream:
            row=json.loads(line);assert (row['position'],row['z_index'])==divmod(count,12) and row['passed']
            errors=[]
            for text,p in zip(row['enclosures'],row['approximation'],strict=True):
                a,b=cusp.endpoints(text);errors.append(max(abs(a-F.from_float(p)),abs(b-F.from_float(p))))
            assert max(errors)==F(row['maximum_error_upper'])<F(plan['total_response_budget']);maximum=max(maximum,max(errors));count+=1
    r=json.loads((OUT/'result.json').read_text());assert count==r['cases']==plan['cases'] and maximum==F(r['maximum_error_upper'])
    print('PASS full momentum response and six partials at certified roots and finite Q samples; continuous Q/full EOS/GR remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
