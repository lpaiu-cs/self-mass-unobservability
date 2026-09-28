"""Certify reusable stellar H coefficients and actual product response pilots."""
from fractions import Fraction as F
from pathlib import Path
import json,shlex,subprocess,sys
import numpy as np
import mpmath as mp
from mpmath import iv
import gr_legendre_log_moments as fraction

rule=fraction.rule;ROOT=rule.ROOT;OUT=rule.OUT.parent/'gr-response-product-defined';native=rule.highq.moments.native;CACHE=native.g.CACHE/'response-product-defined'
I=rule.I;low=rule.low;high=rule.high;cusp=rule.cusp;sha=rule.sha

def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    fraction.verify();assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir();p=json.loads((rule.OUT/'rule.json').read_text())
    lines=[' '.join(x.strip() for x in text[1:-1].split(','))+'\n' for text in p['nodes']]
    lines+=[' '.join(row)+'\n' for row in p['chebyshev_to_legendre']];(OUT/'rule.tsv').write_text(''.join(lines))
    data=dict(np.load(rule.highq.moments.OUT/'inputs.npz'));state=dict(np.load(cusp.THERMO/'states.npz'));source=dict(np.load(cusp.OUT/'inputs.npz'));positions=[0,1603,3205]
    lines=[]
    for i in positions:
        lines.append(' '.join(map(str,[i,int(data['cells'][i]),float(data['beta'][i]).hex(),float(data['Sref'][i]).hex(),float(state['eta'][i]).hex(),float(data['pcut'][i]).hex()]))+'\n')
        (OUT/f'q-{i:04d}.tsv').write_text(''.join(f'{j} {float(Q).hex()}\n' for j,Q in enumerate(source['Q'][i])))
    (OUT/'inputs.tsv').write_text(''.join(lines))
    files=[ROOT/'verification/gr_response_product_native.py',ROOT/'verification/gr_response_product_native.cpp',OUT/'rule.tsv',OUT/'inputs.tsv',
        fraction.OUT/'manifest.json',rule.OUT/'manifest.json',rule.OUT/'tail.json',cusp.THERMO/'states.npz',cusp.OUT/'inputs.npz',rule.highq.moments.OUT/'inputs.npz']+sorted(OUT.glob('q-*.tsv'))
    save('plan.json',dict(classification='Proven',checkpoint='470196f',coefficient_bits=128,response_bits=256,positions=positions,order=32,
        uniform_response_interpolation_budget='1e-18',response_midpoint_budget='1e-15',finite_reference_gate='2e-9',H_control_digits=80,H_control_x=['-1/2','0','1/2'],
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},interval_headers=json.loads((native.OUT/'plan.json').read_text())['interval_headers'],
        target='At the stored binary64 eta reference and exact beta/Sref/P inputs, certify the degree-31 piecewise H/Sref interpolation and its five eta/tau jets on all p in [0,P]. Reuse its coefficients at all twelve fixed Q samples. Add the already certified entire H-tail bound.',
        coefficient_domain='Use exact dyadic panels, R=4h, Re gamma>=sqrt(1+max(c-R,0)^2-R^2)>0 and phase<=1. With a=eta-tmin, qmax=min(1,exp(a)) and Lmax=max(a,0)+exp(min(a,0)), |log(1+exp(eta-t))|<=Lmax by its integral of the occupation. The saved six H jet envelopes yield the Chebyshev remainder 2M/6^32. Bound omissions by the entire envelope. Enforce P*max_panel_error<=1e-18 for each jet.',
        coefficient_rounding='Evaluate the interval Chebyshev nodes, explicit H jets, discrete cosine coefficients and exact rational Legendre connection at 128 bits. Store exact dyadic panel boundaries and outward coefficient endpoints. The interpolation remainder is separate from these coefficient enclosures.',
        response_domain='Evaluate the analytic logarithmic product moments at 256 bits. The continued fraction encloses its unknown terminal ratio in [0,rho] and uses 160 extra levels. Q and panel boundaries remain exact dyadics; derived rational z coordinates are enclosed outward. The response point interval includes polynomial evaluation, P*max_interpolation_error and the entire omitted thermal tail.',
        boundary='Only three selected reference states and twelve Q samples. Eta is the stored reference center, not the exact neutral root. The earlier certified neutral-center/Hessian displacement must be added when this is used for physical fields. No full outer integral or physical EOS claim.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for path,digest in p['interval_headers'].items():assert sha(Path(path))==digest,path
    return p


def build():
    bindings();flags=shlex.split(subprocess.check_output([str(native.CAPD/'build-request15-mp/bin/capd-config'),'--cflags','--libs'],text=True))
    cmd=['g++',str(OUT/'effective-source.cpp'),f'-I{native.DEPS}/include',f'-I{native.DEPS}/include/x86_64-linux-gnu',f'-L{native.DEPS}/lib/x86_64-linux-gnu',*flags,'-O2','-o',str(CACHE/'product')]
    done=subprocess.run(cmd,capture_output=True,text=True);(OUT/'build.log').write_text(done.stdout+done.stderr);save('build.json',dict(command=cmd,returncode=done.returncode));assert done.returncode==0,done.stderr
    linked=subprocess.check_output(['ldd',str(CACHE/'product')],text=True);(OUT/'linked-libraries.txt').write_text(linked)
    libs={word:sha(Path(word)) for line in linked.splitlines() for word in line.split() if word.startswith('/') and Path(word).is_file()}
    save('runtime.json',dict(binary=str(CACHE/'product'),binary_sha256=sha(CACHE/'product'),libraries=libs))


def runtime():
    bindings();r=json.loads((OUT/'runtime.json').read_text());assert sha(Path(r['binary']))==r['binary_sha256']
    for path,digest in r['libraries'].items():assert sha(Path(path))==digest,path
    return r['binary']


def moment_controls():
    binary=runtime();target=OUT/'native-moments.jsonl';assert not target.exists();done=subprocess.run([binary,'moments',str(target)],capture_output=True,text=True)
    (OUT/'moment-control.log').write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr
    actual={r['z']:r for r in map(json.loads,target.read_text().splitlines())};previous=json.loads((rule.OUT/'controls.json').read_text());checks=[]
    for row in previous['checks']:
        a,b=map(F,actual[row['z']]['moments'][row['degree']]);ref=F(row['reference']);passed=a-F('1e-280')<=ref<=b+F('1e-280');assert passed,(row['z'],row['degree'])
        checks.append(dict(z=row['z'],degree=row['degree'],passed=passed))
    save('moment-controls.json',dict(classification='Counterexample candidate',passed=True,checks=checks))


def panels(position):
    for line in (OUT/f'h-{position:04d}.tsv').read_text().splitlines():
        words=line.split();skip=int(words[2]);assert len(words)==(9 if skip else 393)
        coefficients=None if skip else [[(F(words[9+12*n+2*k]),F(words[10+12*n+2*k])) for k in range(6)] for n in range(32)]
        yield dict(a=F(words[0]),b=F(words[1]),skip=bool(skip),errors=list(map(F,words[3:9])),coefficients=coefficients)


def run():
    p=bindings();binary=runtime();iv.prec=256;state=dict(np.load(cusp.THERMO/'states.npz'));data=dict(np.load(rule.highq.moments.OUT/'inputs.npz'));qdata=dict(np.load(cusp.OUT/'inputs.npz'))
    tails=json.loads((rule.OUT/'tail.json').read_text())['rows'];results=[]
    for i in p['positions']:
        polynomial=OUT/f'h-{i:04d}.tsv';assert not polynomial.exists();cmd=[binary,'build',str(OUT/'inputs.tsv'),str(OUT/'rule.tsv'),str(polynomial),str(i)]
        done=subprocess.run(cmd,capture_output=True,text=True);(OUT/f'h-{i:04d}-native.log').write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr
        summary=json.loads(done.stdout);pp=list(panels(i));assert len(pp)==summary['panels'] and pp[0]['a']==0 and pp[-1]['b']==F.from_float(float(data['pcut'][i]))
        assert all(a['b']==b['a'] for a,b in zip(pp,pp[1:]));maximum=[max(x['errors'][j] for x in pp) for j in range(6)]
        assert max(maximum)*F.from_float(float(data['pcut'][i]))<F(p['uniform_response_interpolation_budget'])
        target=OUT/f'response-{i:04d}.jsonl';assert not target.exists();cmd=[binary,'response',str(OUT/'inputs.tsv'),str(OUT/'rule.tsv'),str(target),str(i),str(polynomial),str(OUT/f'q-{i:04d}.tsv')]
        done=subprocess.run(cmd,capture_output=True,text=True);(OUT/f'response-{i:04d}-native.log').write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr
        timing=json.loads(done.stdout);rows=[];max_error=F(0)
        for j,line in enumerate(target.read_text().splitlines()):
            row=json.loads(line);assert row['position']==i and row['z_index']==j and row['cell']==int(data['cells'][i]);Q=I(qdata['Q'][i,j]);enclosures=[];points=[];errors=[]
            for (a,b),bound in zip(row['finite_response'],tails[i]['upper_coefficients']):
                finite=iv.mpf([I(a).a,I(b).b]);total=finite+cusp.symmetric(high(I(bound)/I(max(F(1),low(Q*Q)))))
                text=cusp.interval_text(total);lo,hi=cusp.endpoints(text);point=float((lo+hi)/2);error=max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point)))
                enclosures.append(text);points.append(point);errors.append(str(error));max_error=max(max_error,error)
            rows.append(dict(z_index=j,enclosures=enclosures,approximation=points,error_upper=errors))
        assert len(rows)==12
        values=rule.highq.neutral.thermo.response(np.full(12,float(state['eta'][i])),np.full(12,float(data['beta'][i])),qdata['Q'][i],np.full(12,float(data['Sref'][i])),2e-13)
        score=float(np.max(np.abs(values-np.array([r['approximation'] for r in rows]).T))/max(np.max(np.abs(values)),1e-300))
        results.append(dict(position=i,cell=int(data['cells'][i]),passed=max_error<F(p['response_midpoint_budget']),rows=rows,
            maximum_midpoint_error=str(max_error),uniform_H_error=list(map(str,maximum)),panels=summary['panels'],bounded_omissions=summary['bounded_omissions'],
            coefficient_seconds=summary['seconds'],response_seconds=timing['seconds'],reference_control=dict(classification='Counterexample candidate',relative_vector_difference=score,passed=score<float(p['finite_reference_gate']))))
        print('REUSABLE H RESPONSE',i,summary['panels'],summary['seconds'],timing['seconds'],float(max_error),score,flush=True)
    save('result.json',dict(classification='Proven',passed=all(x['passed'] for x in results),results=results,
        exact_neutral_root_certified=False,full_outer_integral_certified=False,physical_EOS_certified=False))


def H_controls():
    p=bindings();iv.prec=256;mp.mp.dps=p['H_control_digits'];state=dict(np.load(cusp.THERMO/'states.npz'));data=dict(np.load(rule.highq.moments.OUT/'inputs.npz'));checks=[]
    def M(x):
        q=F(x);return mp.mpf(q.numerator)/q.denominator
    for i in p['positions']:
        pp=list(panels(i));eta=mp.mpf(float(state['eta'][i]));beta=mp.mpf(float(data['beta'][i]));Sref=mp.mpf(float(data['Sref'][i]))
        for index in [0,len(pp)//2,len(pp)-1]:
            panel=pp[index]
            for v in map(F,p['H_control_x']):
                momentum=M((panel['a']+panel['b'])/2+(panel['b']-panel['a'])/2*v);gamma=mp.sqrt(1+momentum*momentum);t=momentum*momentum/(beta*(gamma+1));q=1/(1+mp.exp(t-eta));w=1-q;L=mp.log1p(mp.exp(eta-t))
                actual=[q/gamma+2*beta*L,q*w/gamma+2*beta*q,q*w*(1-2*q)/gamma+2*beta*q*w,t*q*w/gamma+2*beta*(L+t*q),
                    t*q*w*(1-2*q)/gamma+2*beta*(q+t*q*w),(t*t*(1-2*q)-t)*q*w/gamma+2*beta*(L+t*q+t*t*q*w)]
                legendre=[I(1),I(v)]
                for n in range(1,31):legendre.append(((2*n+1)*I(v)*legendre[n]-n*legendre[n-1])/(n+1))
                enclosed=[I(0)]*6
                if not panel['skip']:
                    for n,values in enumerate(panel['coefficients']):
                        for k,(a,b) in enumerate(values):enclosed[k]+=legendre[n]*iv.mpf([I(a).a,I(b).b])
                enclosed=[x+cusp.symmetric(e) for x,e in zip(enclosed,panel['errors'])]
                inside=all(M(low(x))<=y/Sref<=M(high(x)) for x,y in zip(enclosed,actual));assert inside,(i,index,v)
                checks.append(dict(position=i,panel=index,x=str(v),passed=inside,reference=[str(x/Sref) for x in actual]))
    save('H-controls.json',dict(classification='Counterexample candidate',passed=True,components=6*len(checks),checks=checks))


def finalize():
    runtime()
    for name in ['moment-controls.json','H-controls.json','result.json']:assert json.loads((OUT/name).read_text())['passed'],name
    assert all(x['reference_control']['passed'] for x in json.loads((OUT/'result.json').read_text())['results'])
    save('manifest.json',dict(sha256={path.relative_to(ROOT).as_posix():sha(path) for path in OUT.iterdir() if path.is_file()}));verify()


def verify():
    p=bindings();runtime()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for name in ['moment-controls.json','H-controls.json','result.json']:assert json.loads((OUT/name).read_text())['passed'],name
    r=json.loads((OUT/'result.json').read_text());assert [x['position'] for x in r['results']]==p['positions']
    for result in r['results']:
        assert result['reference_control']['passed']
        for row in result['rows']:
            for text,point,bound in zip(row['enclosures'],row['approximation'],row['error_upper'],strict=True):
                lo,hi=cusp.endpoints(text);assert max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point)))==F(bound)<F(p['response_midpoint_budget'])
    print('PASS three actual H interpolation tables and 36 full momentum product responses at stored reference centers; physical neutral displacement/outer integration remain separate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
