"""Extend the frozen true-neutral path and node checks to completed state 1603."""
from fractions import Fraction as F
import gzip,json,sys
from mpmath import iv
import gr_outer_node_comparison as prior

paths=prior.paths;ROOT=paths.ROOT;sha=paths.sha;cusp=prior.cusp
OUT=paths.OUT.parent/'gr-outer-second-path-comparison';POSITION=1603


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def inputs():
    stem=f'cell-{POSITION:04d}'
    return stem,paths.direct.OUT,paths.product.OUT


def prepare():
    assert not OUT.exists();prior.verify();paths.product.verify()
    text,_=paths.direct.source();assert text==(paths.direct.OUT/'candidate.py').read_text()
    obj=paths.direct.module(text);obj.bindings();obj.runtime()
    stem,ddir,pdir=inputs();files=[ROOT/'verification/gr_outer_second_path_comparison.py',
        ROOT/'verification/gr_outer_node_comparison.py',ROOT/'verification/gr_outer_path_comparison.py',
        prior.OUT/'plan.json',prior.OUT/'manifest.json',paths.product.OUT/'manifest.json',
        paths.direct.OUT/'preflight-manifest.json']
    for folder in [ddir,pdir]:
        current=list(folder.glob(stem+'-*'));assert len(current)==6;files+=current
    a=json.loads((ddir/(stem+'-result.json')).read_text());b=json.loads((pdir/(stem+'-result.json')).read_text())
    assert a['passed'] and b['passed'] and a['position']==b['position']==POSITION
    OUT.mkdir();save('plan.json',dict(classification='Proven',position=POSITION,expected_nodes=19024,
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        field_gate='All seven complete true-neutral physical-field intervals intersect, with each original absolute error below 2e-7.',
        node_gate='All 19024 Q values agree exactly. At each node, the three eta responses intersect after the previously proved neutral-center shift. No fitted tolerance.',
        theorem_source=(prior.OUT/'plan.json').relative_to(ROOT).as_posix(),
        scope='Two momentum representations share the outer Gauss rule, analytic domains, IR/high-Q bounds and field definitions. This finite compatibility check does not provide a full physical EOS or GR trajectory.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return plan


def calculate():
    plan=bindings();stem,ddir,pdir=inputs();iv.prec=256
    a=json.loads((ddir/(stem+'-result.json')).read_text());b=json.loads((pdir/(stem+'-result.json')).read_text())
    field=paths.compare(a,b)
    for item in [a,b]:
        for text,point,error in zip(item['complete_field_enclosures'],item['field_approximation'],item['field_error_upper'],strict=True):
            lo,hi=cusp.endpoints(text);assert max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point)))==F(error)<F('2e-7')
    delta=F(b['eta_center_distance_upper']);factor=cusp.high(iv.exp(cusp.I(delta)));inverse=cusp.low(iv.exp(-cusp.I(delta)))
    minimum=[None]*3;disjoint=[0]*3;count=0
    with gzip.open(ddir/(stem+'-points.jsonl.gz'),'rt') as direct,gzip.open(pdir/(stem+'-points.jsonl.gz'),'rt') as product:
        for i,(x,y) in enumerate(zip(direct,product,strict=True)):
            x=json.loads(x);y=json.loads(y);assert x['node']==y['node']==i and F(x['Q_exact'])==F(y['Q_exact'])
            d=[cusp.endpoints(t) for t in x['response'][:3]];p=[cusp.endpoints(t) for t in y['response'][:3]];assert d[0][0]>0
            shift=delta*factor*d[0][1];enclosed=[(inverse*d[0][0],factor*d[0][1])]+[(lo-shift,hi+shift) for lo,hi in d[1:]]
            for k,((lo,hi),(pl,ph)) in enumerate(zip(enclosed,p,strict=True)):
                width=min(hi,ph)-max(lo,pl);minimum[k]=width if minimum[k] is None else min(minimum[k],width);disjoint[k]+=width<0
            count+=1
    bad=json.loads(json.dumps(b));bad['complete_field_enclosures'][0]='[1, 2]';assert not paths.compare(a,bad)['passed']
    return dict(classification='Counterexample candidate',passed=field['passed'] and count==plan['expected_nodes'] and not any(disjoint),
        field_comparison=field,nodes=count,components=3*count,disjoint_counts=disjoint,minimum_intersection_widths=list(map(str,minimum)),
        eta_center_distance_upper=str(delta),exp_delta_upper=str(factor),exp_minus_delta_lower=str(inverse),
        disjoint_field_negative_control_passed=True,physical_EOS_certified=False,nonlinear_GR_evolution=False)


def run():
    assert not (OUT/'result.json').exists();result=calculate();save('result.json',result);assert result['passed']
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    plan=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result['passed'] and result['nodes']==plan['expected_nodes']
    assert result['components']==3*result['nodes'] and result['disjoint_counts']==[0,0,0]
    assert result['field_comparison']['passed'] and min(map(F,result['minimum_intersection_widths']))>=0
    assert not result['physical_EOS_certified'] and not result['nonlinear_GR_evolution']
    print('PASS state1603: seven true-neutral fields and',result['components'],'exact-Q, center-adjusted response comparisons',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
