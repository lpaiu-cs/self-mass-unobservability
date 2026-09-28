"""Certify reusable truncated moment jets for the high-wavenumber expansion."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
from pathlib import Path
import gzip,json,shlex,shutil,subprocess,sys
import numpy as np
import mpmath as mp
import sympy as sp
import gr_response_momentum_tail as tail
import gr_plasma_interval as native

cusp=tail.cusp;ROOT=cusp.ROOT;OUT=cusp.OUT.parent/'gr-highq-moments';CACHE=native.g.CACHE/'highq-moments'
RULE=cusp.OUT.parent/'gr-smooth-response-defined';sha=cusp.rule.sha
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    tail.verify();assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    data=dict(np.load(cusp.OUT/'inputs.npz'));cuts=[r for r in tail.records() if r['z_index']==0];assert len(cuts)==3206
    np.savez_compressed(OUT/'inputs.npz',cells=data['cells'],beta=data['beta'],Sref=data['Sref'],root_intervals=data['root_intervals'],pcut=np.array([r['pcut'] for r in cuts]))
    lines=[]
    for i,row in enumerate(cuts):
        assert row['position']==i and row['cell']==data['cells'][i]
        lo,hi=str(data['root_intervals'][i])[1:-1].split(',');lines.append(' '.join(map(str,[i,row['cell'],float(data['beta'][i]).hex(),float(data['Sref'][i]).hex(),lo.strip(),hi.strip(),row['pcut_hex']]))+'\n')
    (OUT/'inputs.tsv').write_text(''.join(lines));shutil.copy2(RULE/'rule.tsv',OUT/'rule.tsv')
    files=[ROOT/'verification/gr_highq_moments.py',ROOT/'verification/gr_highq_moments.cpp',tail.OUT/'manifest.json',tail.OUT/'bounds.jsonl.gz',
        cusp.OUT/'manifest.json',cusp.OUT/'inputs.npz',RULE/'rule-manifest.json',RULE/'rule.json',OUT/'inputs.npz',OUT/'inputs.tsv',OUT/'rule.tsv']
    save('plan.json',dict(classification='Proven',checkpoint='caec43c',bits=128,cells=3206,moment_orders=list(range(1,18)),parameter_order=cusp.PARAMETERS,
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},interval_headers=json.loads((native.OUT/'plan.json').read_text())['interval_headers'],
        processes=3,block_size=32,component_error_budget='1e-12',remainder_budget='1e-13',high_Q_series_order=16,high_Q_minimum_ratio=4,
        control_positions=[0,1603,3205],control_moment_orders=[1,9,17],control_digits=70,
        target='For each certified neutral eta interval, certify mu_m=integral_0^pcut (p/pcut)^(2m)*f(p)/(gamma*Sref) dp and its five eta/tau partials, m=1..17. pcut is the exact upward binary cutoff from the earlier momentum-tail certificate and is held fixed during differentiation.',
        quadrature='Reuse the exact 16-point rule. Adaptive dyadic panels with R=4h require g2>0 and phase<=1, then use Cauchy M32<=32!*M/(R-h)^32. The analytic integrand majorant is qmax*(pmax/pcut)^(2m)/(gmin*Sref) times the six existing occupation bounds. Allocate 1e-13*width/pcut to each moment component and include all directed evaluation rounding.',
        future_series='For Q>=4pcut, S(Q)/Sref on the truncated momentum domain is approximated by sum_{m=1}^16 (pcut/Q)^(2m)*[mu_m/(2m-1)+2pcut^2*mu_(m+1)/((2m-1)*(2m+1))]. The positive geometric remainder, omitted momentum tail and nonlinear self-integrand propagation must still be certified separately.',
        physical_EOS_certified=False,high_Q_response_certified=False))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for rel,digest in p['interval_headers'].items():assert sha(Path(rel))==digest,rel
    return p


def symbolic():
    m=sp.symbols('m',integer=True,positive=True)
    assert sp.simplify(1/(2*m-1)-1/(2*m+1)-2/((2*m-1)*(2*m+1)))==0
    x,p=sp.symbols('x p',positive=True)
    h=sum((p*p*x)**n/sp.Rational(2*n+1) for n in range(6));j=sum(2*(p*p*x)**n/sp.Rational((2*n-1)*(2*n+1)) for n in range(1,7))
    kernel=p*p*x*(h+j/x)
    for k in range(1,6):assert sp.expand(kernel).coeff(x,k)==p**(2*k)/sp.Rational(2*k-1)+2*p**(2*k+2)/sp.Rational((2*k-1)*(2*k+1))
    save('symbolic.json',dict(classification='Proven',passed=True,
        positive_series='From h(u)=sum u^(2n)/(2n+1) and j(u)=sum_{n>=1} 2u^(2n)/((2n-1)*(2n+1)), K_p=gamma^-1*sum_{m>=1} Q^(-2m)*[p^(2m)/(2m-1)+2p^(2m+2)/((2m-1)*(2m+1))]. This converges absolutely for p<Q.',
        normalized_moments='Holding pcut and Sref fixed, multiplying each occupation derivative by (p/pcut)^(2m)/gamma produces exactly the same six parameter jets of mu_m. The coefficient formula holds for all six truncated response components.',
        scope='Moment identities and the truncated interval certificate. This stage does not yet assert an error budget for the complete high-Q response or its nonlinear outer integral.'))


def build():
    bindings();symbolic();flags=shlex.split(subprocess.check_output([str(native.CAPD/'build-request15-mp/bin/capd-config'),'--cflags','--libs'],text=True))
    cmd=['g++',str(ROOT/'verification/gr_highq_moments.cpp'),f'-I{native.DEPS}/include',f'-I{native.DEPS}/include/x86_64-linux-gnu',f'-L{native.DEPS}/lib/x86_64-linux-gnu',*flags,'-O2','-o',str(CACHE/'moments')]
    done=subprocess.run(cmd,capture_output=True,text=True);(OUT/'build.log').write_text(done.stdout+done.stderr);save('build.json',dict(command=cmd,returncode=done.returncode));assert done.returncode==0,done.stderr
    linked=subprocess.check_output(['ldd',str(CACHE/'moments')],text=True);(OUT/'linked-libraries.txt').write_text(linked)
    libs={word:sha(Path(word)) for line in linked.splitlines() for word in line.split() if word.startswith('/') and Path(word).is_file()}
    save('runtime.json',dict(binary=str(CACHE/'moments'),binary_sha256=sha(CACHE/'moments'),libraries=libs))


def runtime():
    bindings();r=json.loads((OUT/'runtime.json').read_text());assert sha(Path(r['binary']))==r['binary_sha256']
    for path,digest in r['libraries'].items():assert sha(Path(path))==digest,path
    return r['binary']


def block(first,count):
    path=OUT/f'block-{first:04d}.jsonl';assert not path.exists() and not path.with_suffix('.jsonl.gz').exists()
    done=subprocess.run([str(CACHE/'moments'),str(OUT/'inputs.tsv'),str(OUT/'rule.tsv'),str(path),str(first),str(count)],capture_output=True,text=True)
    if done.returncode:(OUT/f'block-{first:04d}-error.log').write_text(done.stdout+done.stderr)
    assert done.returncode==0,done.stderr
    with gzip.open(path.with_suffix('.jsonl.gz'),'xb',compresslevel=1) as stream:stream.write(path.read_bytes())
    path.unlink();return first


def records():
    for path in sorted(OUT.glob('block-*.jsonl.gz')):
        with gzip.open(path,'rt',encoding='utf-8') as stream:
            for line in stream:yield json.loads(line)


def run():
    p=bindings();runtime();finished=[]
    with ProcessPoolExecutor(max_workers=p['processes']) as pool:
        work=[pool.submit(block,i,min(p['block_size'],p['cells']-i)) for i in range(0,p['cells'],p['block_size'])]
        for future in as_completed(work):
            finished.append(future.result());save('progress.json',dict(completed_blocks=sorted(finished)));print('CERTIFIED MOMENT BLOCKS',len(finished),flush=True)
    maximum=F(0);count=0;panels=0;omitted=0
    for row in records():
        assert row['position']==count and len(row['moments'])==17;count+=1;panels+=row['panels'];omitted+=row['bounded_omissions']
        assert F(row['maximum_remainder_upper'])<F(p['remainder_budget'])
        for moment in row['moments']:
            assert len(moment)==6
            for a,b in moment:
                lo,hi=F(a),F(b);mid=F.from_float(float((lo+hi)/2));assert lo<=hi;maximum=max(maximum,abs(mid-lo),abs(hi-mid))
    assert count==p['cells']
    save('result.json',dict(classification='Proven',passed=maximum<F(p['component_error_budget']),cells=count,moment_components=count*17*6,
        maximum_midpoint_error_exact=str(maximum),panels=panels,bounded_omissions=omitted,display_only=dict(maximum_error=float(maximum)),
        complete_high_Q_response_certified=False,outer_integral_certified=False,physical_EOS_certified=False))


def controls():
    p=bindings();mp.mp.dps=p['control_digits'];data=dict(np.load(OUT/'inputs.npz'));checks=[]
    def M(x):
        q=F(x);return mp.mpf(q.numerator)/q.denominator
    for row in records():
        i=row['position']
        if i not in p['control_positions']:continue
        eta=sum(map(M,cusp.endpoints(str(data['root_intervals'][i]))))/2;beta=mp.mpf(float(data['beta'][i]));Sref=mp.mpf(float(data['Sref'][i]));pcut=mp.mpf(float(data['pcut'][i]))
        for order in p['control_moment_orders']:
            values=[]
            for k in range(6):
                def f(v):
                    momentum=pcut*v;gamma=mp.sqrt(1+momentum*momentum);t=momentum*momentum/(beta*(gamma+1));q=1/(1+mp.exp(t-eta));z=1-q
                    occupation=[q,q*z,q*z*(1-2*q),t*q*z,t*q*z*(1-2*q),(t*t*(1-2*q)-t)*q*z][k]
                    return pcut*v**(2*order)*occupation/(gamma*Sref)
                values.append(mp.quad(f,[0,mp.mpf('.25'),mp.mpf('.5'),1]))
            inside=all(M(a)<=v<=M(b) for (a,b),v in zip(row['moments'][order-1],values))
            checks.append(dict(cell=row['cell'],moment_order=order,passed=inside,reference=list(map(str,values))));print('MOMENT INDEPENDENT CONTROL',row['cell'],order,inside,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=all(c['passed'] for c in checks),checks=checks));assert len(checks)==9 and all(c['passed'] for c in checks)


def finalize():
    runtime()
    for name in ['symbolic.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    p=bindings();runtime()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for name in ['symbolic.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    count=0;maximum=F(0)
    for row in records():
        assert row['position']==count;count+=1
        for moment in row['moments']:
            for a,b in moment:
                lo,hi=F(a),F(b);mid=F.from_float(float((lo+hi)/2));maximum=max(maximum,abs(mid-lo),abs(hi-mid))
    r=json.loads((OUT/'result.json').read_text());assert count==r['cells']==3206 and maximum==F(r['maximum_midpoint_error_exact'])<F(p['component_error_budget'])
    print('PASS 3206 x17 x6 truncated thermal moment enclosures; high-Q tail/nonlinear outer propagation remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
