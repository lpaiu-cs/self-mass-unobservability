"""Audit actual node identity and neutral-center response compatibility."""
from fractions import Fraction as F
import gzip,json,sys
import sympy as sp
from mpmath import iv
import gr_outer_path_comparison as paths

ROOT=paths.ROOT;OUT=paths.OUT.parent/'gr-outer-node-comparison';cusp=paths.product.cusp;sha=paths.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def run():
    assert not OUT.exists();paths.verify();OUT.mkdir();iv.prec=256
    pdir=paths.product.OUT;ddir=paths.direct.OUT
    pr=json.loads((pdir/'cell-0000-result.json').read_text());delta=F(pr['eta_center_distance_upper'])
    exponential=iv.exp(cusp.I(delta));factor=cusp.high(exponential);inverse=cusp.low(iv.exp(-cusp.I(delta)))
    f=sp.symbols('f',real=True);polynomial=1-6*f+6*f*f
    assert sp.expand(polynomial+sp.Rational(1,2)-6*(f-sp.Rational(1,2))**2)==0
    assert sp.expand(1-polynomial-6*f*(1-f))==0
    files=[ROOT/'verification/gr_outer_node_comparison.py',paths.OUT/'manifest.json',pdir/'cell-0000-result.json',
        pdir/'cell-0000-points.jsonl.gz',ddir/'cell-0000-points.jsonl.gz']
    save('plan.json',dict(classification='Proven',position=0,expected_nodes=pr['nodes'],bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        eta_center_distance_upper=str(delta),exp_delta_upper=str(factor),exp_minus_delta_lower=str(inverse),
        theorem='At fixed beta,Q,Sref the positive momentum kernel is independent of eta. For f in [0,1], |f_eta|,|f_etaeta|,|f_etaetaeta|<=f and f(eta+d)<=exp(|d|)*f(eta). Thus R0(center) lies in exp([-delta,delta])*R0(root); R_eta and R_etaeta shift by at most delta*exp(delta)*R0(root). The etaetaeta polynomial has range [-1/2,1] by the two exact identities checked here.',
        gate='Both actual Q values must be the same exact dyadic at every index. The three center-adjusted direct response intervals must intersect the product intervals. No tolerance is fitted to the results.',
        boundary='A necessary finite compatibility check on two interval algorithms. The tau derivatives and all possible states are not tested here; shared endpoints/domains remain shared.'))
    minimum=[None]*3;disjoint=[0]*3;count=0
    with gzip.open(ddir/'cell-0000-points.jsonl.gz','rt') as direct,gzip.open(pdir/'cell-0000-points.jsonl.gz','rt') as product:
        for i,(a,b) in enumerate(zip(direct,product,strict=True)):
            a=json.loads(a);b=json.loads(b);assert a['node']==b['node']==i and F(a['Q_exact'])==F(b['Q_exact'])
            d=[cusp.endpoints(x) for x in a['response'][:3]];p=[cusp.endpoints(x) for x in b['response'][:3]];assert d[0][0]>0
            shift=delta*factor*d[0][1];enclosed=[(inverse*d[0][0],factor*d[0][1])]+[(lo-shift,hi+shift) for lo,hi in d[1:]]
            for k,((lo,hi),(pl,ph)) in enumerate(zip(enclosed,p,strict=True)):
                width=min(hi,ph)-max(lo,pl);minimum[k]=width if minimum[k] is None else min(minimum[k],width);disjoint[k]+=width<0
            count+=1
    result=dict(classification='Counterexample candidate',passed=count==pr['nodes'] and not any(disjoint),nodes=count,components=3*count,
        disjoint_counts=disjoint,minimum_intersection_widths=list(map(str,minimum)),physical_EOS_certified=False)
    save('result.json',result);assert result['passed'],result
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    p=json.loads((OUT/'plan.json').read_text());m=json.loads((OUT/'manifest.json').read_text())
    for rel,digest in {**p['bindings'],**m['sha256']}.items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['nodes']==p['expected_nodes'] and r['components']==3*r['nodes'] and r['disjoint_counts']==[0,0,0]
    assert min(map(F,r['minimum_intersection_widths']))>=0
    print('PASS',r['nodes'],'exact common Q nodes and',r['components'],'center-adjusted independent response interval comparisons',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
