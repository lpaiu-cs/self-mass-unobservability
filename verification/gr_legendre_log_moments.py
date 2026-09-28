"""Bound all product-integration log moments with one interval continued fraction."""
from fractions import Fraction as F
import json,statistics,sys,time
import sympy as sp
from mpmath import iv
import gr_response_product_rule as rule

ROOT=rule.ROOT;OUT=rule.OUT.parent/'gr-legendre-log-moments';I=rule.I;low=rule.low;high=rule.high;cusp=rule.cusp;sha=rule.sha

def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def all_moments(z,order=32,padding=160):
    z=F(z)
    if abs(z)<F(5,4):return [rule.logarithmic_moment(n,z) for n in range(order+1)]
    a=I(abs(z));rho=1/(a+iv.sqrt(a*a-1));r=iv.mpf([0,I(high(rho)).b]);ratios=[None]*(order+1)
    for n in range(order+1+padding,0,-1):
        denominator=(2*n+1)*a-(n+1)*r;assert low(denominator)>0;r=n/denominator
        if n-1<=order:ratios[n-1]=r
    Q=[iv.log((a+1)/(a-1))/2]
    for r in ratios:Q.append(Q[-1]*r)
    out=[(a+1)*iv.log(a+1)-(a-1)*iv.log(a-1)-2]
    for n in range(1,order+1):out.append(2*(Q[n+1]-Q[n-1])/(2*n+1)*(-1 if z<0 and n%2 else 1))
    return out


def prepare():
    rule.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_legendre_log_moments.py',rule.OUT/'manifest.json',rule.OUT/'rule.json',rule.OUT/'controls.json']
    save('plan.json',dict(classification='Proven',checkpoint='396c94e',bits=256,maximum_degree=32,far_threshold='5/4',padding=160,
        comparison_allowance='1e-280',benchmark_z=['2','10','1000000'],benchmark_repetitions=5,
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        scope='Same real logarithmic product moments. This certifies and measures their batch evaluator; actual stellar H interpolation coefficients and outer product integration remain separate.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def run():
    p=bindings();iv.prec=p['bits'];t,r=sp.symbols('t r',real=True);z=(r*r+1)/(2*r)
    assert sp.cancel(r-(1-t*t)/(2*(z-t))-(t-r)**2/(2*(z-t)))==0
    checks=[];old=json.loads((rule.OUT/'controls.json').read_text());cache={}
    for row in old['checks']:
        if row['z'] not in cache:cache[row['z']]=all_moments(row['z'],p['maximum_degree'],p['padding'])
        value=cache[row['z']][row['degree']];text=cusp.interval_text(value);lo,hi=cusp.endpoints(text);ref=F(row['reference']);margin=F(p['comparison_allowance'])
        passed=lo-margin<=ref<=hi+margin;assert passed,(row['degree'],row['z'])
        a,b=cusp.endpoints(row['enclosure']);assert max(a,lo)<=min(b,hi)
        checks.append(dict(degree=row['degree'],z=row['z'],enclosure=text,reference_inclusion_with_allowance=passed,old_interval_overlap=True))
    timings=[]
    for z in p['benchmark_z']:
        before=[];after=[]
        for _ in range(p['benchmark_repetitions']):
            start=time.perf_counter();[rule.logarithmic_moment(n,z) for n in range(p['maximum_degree']+1)];before.append(time.perf_counter()-start)
            start=time.perf_counter();all_moments(z,p['maximum_degree'],p['padding']);after.append(time.perf_counter()-start)
        timings.append(dict(classification='Counterexample candidate',z=z,series_seconds=before,fraction_seconds=after,median_speedup=statistics.median(before)/statistics.median(after)))
    save('result.json',dict(classification='Proven',passed=True,checks=checks,benchmarks=timings,
        positive_ratio='Rodrigues and n integrations by parts give Qn(z)=2^(-n-1)*integral_-1^1 (1-t^2)^n/(z-t)^(n+1)dt for z>1. Hence Q_(n+1)/Qn is a positive weighted average of (1-t^2)/(2(z-t)). If rho=z-sqrt(z^2-1), rho-(1-t^2)/(2(z-t))=(t-rho)^2/(2(z-t))>=0. Thus every terminal ratio lies in [0,rho].',
        recurrence='The exact Legendre recurrence gives r_(n-1)=n/[(2n+1)z-(n+1)r_n]. Initialize the distant terminal ratio with the entire rigorous interval [0,rho], then propagate backward with directed interval arithmetic and positive denominators. No convergence estimate or unproved zero terminal condition is used. The fixed padding controls width, while validity holds for every nonnegative padding.',
        normalization='Q0=log((z+1)/(z-1))/2 fixes the normalization. Multiply by the enclosed ratios to obtain all Qn, then all Jn from the already certified identity. Negative z follows Jn(-z)=(-1)^n Jn(z). For |z|<5/4 retain the existing polynomial recurrence with exact endpoint limits.',
        control_scope='All 54 stored independent 320-digit controls agree with the new intervals within the original 1e-280 numerical allowance and overlap the earlier rigorous intervals. The benchmark compares all 33 moments together at 256 bits and is a local finite runtime measurement only.',
        actual_stellar_interpolation_certified=False,full_outer_integral_certified=False))
    save('manifest.json',dict(sha256={path.relative_to(ROOT).as_posix():sha(path) for path in OUT.iterdir() if path.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and len(r['checks'])==54
    assert all(x['reference_inclusion_with_allowance'] and x['old_interval_overlap'] for x in r['checks'])
    print('PASS interval continued-fraction log moments; local median speedups:',[x['median_speedup'] for x in r['benchmarks']],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
