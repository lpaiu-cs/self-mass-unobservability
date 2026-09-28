"""Match the original face radius before evaluating the frozen thermal gradient."""
from fractions import Fraction as F
from functools import reduce
from operator import mul
import json,subprocess,sys
import numpy as np
from mpmath import iv
import gr_heat_boundary_compatibility as original

ROOT=original.ROOT;OUT=original.OUT.parent/'gr-heat-boundary-radius'
sha=original.sha;I=original.I;low=original.low;high=original.high;exact=original.exact;ld=original.ld


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def box(a,b):return iv.mpf([I(a).a,I(b).b])


def basis(x,s,weights,interval=False):
    lift=I if interval else lambda q:q
    values=[];derivatives=[]
    for j,w in enumerate(weights):
        value=lift(w)*reduce(mul,(x-lift(s[k]) for k in range(len(s)) if k!=j),lift(1))
        derivative=value*sum((1/(x-lift(s[k])) for k in range(len(s)) if k!=j),lift(0))
        values.append(value);derivatives.append(derivative)
    return values,derivatives


def profile(cell,n):
    metric=original.heat.original.metric
    reference=json.loads((metric.OUT/'plan.json').read_text());path=ROOT/reference['cell_sources'][str(cell)]
    data=metric.node_file(str(path.with_name(path.stem+f'-nodes-{n}.npz')))
    j=list(data['cells']).index(cell);order=np.argsort(data['radius_cm'][j]);a=data['eos'][j,order]
    r=list(map(exact,data['radius_cm'][j,order]));t=list(map(exact,data['lnT'][j,order]))
    c=ld(metric.g.c.gr.C)*100;C=exact(ld(data['C_X'][j])*c*c)
    H=[C+exact(row[2])+exact(row[1])/exact(row[0]) for row in a];assert min(H)>0
    s=original.endpoint_rule(n,1 if cell==0 else 0)[0]
    weights=[1/reduce(mul,(s[j]-s[k] for k in range(n) if k!=j),F(1)) for j in range(n)]
    grid=np.load(metric.g.OUT/'gr-increment-structure/path-4.npz')
    target=exact(ld(grid['radius_m'][0 if cell==0 else 5735])*100)
    return s,weights,r,t,H,target


def prepare():
    assert not OUT.exists();OUT.mkdir();prior=json.loads((original.OUT/'plan.json').read_text())
    files=[ROOT/'verification/gr_heat_boundary_radius.py',ROOT/'verification/gr_heat_boundary_compatibility.py',
        original.OUT/'manifest.json',original.OUT/'result.json']
    prior['bindings'].update({p.relative_to(ROOT).as_posix():sha(p) for p in files})
    prior.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        root_bracket_powers=list(range(50,7,-1)),bisections=96,
        root_certificate='Exact rational opposite signs of r(s)-r_face, plus an outward strictly positive r_s bound on the whole initial bracket, prove one root in that bracket. Bisect the exact polynomial 96 times; enclose the thermal derivative on the entire retained root interval.',
        preservation='Use exactly the previous dyadic nodes and nodal thermodynamic numbers. The same polynomial is evaluated at a new parameter value; no node, radius, temperature or heat flux is fitted. Record if reaching the original face requires continuing that polynomial outside its original parameter cell [0,1].',
        boundary='This removes the endpoint-coordinate mismatch within the declared polynomial only. Its continuation, native EOS errors and the true physical profile remain uncertified. No global polynomial-root uniqueness or GR trajectory claim.')
    save('plan.json',prior)
    checks=[]
    for n in prior['nodes']:
        s=original.endpoint_rule(n,1)[0]
        weights=[1/reduce(mul,(s[j]-s[k] for k in range(n) if k!=j),F(1)) for j in range(n)]
        for x in [F(-1,1024),F(1025,1024)]:
            V,D=basis(x,s,weights)
            for k in range(n):
                assert sum((v*q**k for v,q in zip(V,s)),F(0))==x**k
                assert sum((d*q**k for d,q in zip(D,s)),F(0))==(0 if k==0 else k*x**(k-1))
        checks.append(dict(nodes=n,exact_polynomial_controls=True))
    save('controls.json',dict(classification='Proven',passed=True,rows=checks))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def evaluate(cell,n,end,bits,plan):
    iv.prec=bits;s,weights,r,t,H,target=profile(cell,n);dr=[q-r[0] for q in r]
    def value(x):return r[0]+sum((v*y for v,y in zip(basis(x,s,weights)[0],dr)),F(0))-target
    def radial_derivative(a,b):
        return sum((d*I(y) for d,y in zip(basis(box(a,b),s,weights,True)[1],dr)),I(0))
    e=F(end);found=None
    for power in plan['root_bracket_powers']:
        a,b=e-F(1,2**power),e+F(1,2**power)
        if any(a<=q<=b for q in s):break
        if value(a)<=0<=value(b):found=(a,b);break
    assert found is not None,(cell,n,'no root bracket before an interpolation node')
    a,b=found;initial_derivative=radial_derivative(a,b);assert low(initial_derivative)>0
    for _ in range(plan['bisections']):
        m=(a+b)/2;v=value(m)
        if v==0:a=b=m;break
        if v<0:a=m
        else:b=m
    va,vb=value(a),value(b);assert va<=0<=vb
    _,D=basis(box(a,b),s,weights,True);rs=sum((d*I(y) for d,y in zip(D,dr)),I(0));assert low(rs)>0
    potential=[I(v-t[0])-iv.log(I(h/H[0])) for v,h in zip(t,H)]
    gradient=sum((d*y for d,y in zip(D,potential)),I(0))/rs
    lo,hi=low(gradient),high(gradient)
    return dict(cell=cell,nodes=n,end=end,precision_bits=bits,original_face_radius_cm=str(target),
        initial_parameter_bracket=list(map(str,found)),initial_radial_derivative_lower=str(low(initial_derivative)),
        root_parameter_bounds=[str(a),str(b)],root_radius_residual_bounds=[str(va),str(vb)],
        unique_root_in_initial_bracket=True,parameter_continuation_outside_cell=a<0 or b>1,
        gradient_per_cm=[str(lo),str(hi)],zero_excluded=lo>0 or hi<0)


def run():
    plan=bindings();assert not (OUT/'result.json').exists();rows=[]
    for target in plan['endpoints']:
        for n in plan['nodes']:
            pair=[evaluate(target['cell'],n,target['end'],bits,plan) for bits in plan['bits']]
            assert pair[0]['root_parameter_bounds']==pair[1]['root_parameter_bounds']
            a,b=[list(map(F,r['gradient_per_cm'])) for r in pair]
            assert max(a[0],b[0])<=min(a[1],b[1])
            rows.append(dict(label=target['label'],evaluations=pair,precision_intersection=True))
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        physical_boundary_certified=False,native_or_continuous_derivative_error_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    assert json.loads((OUT/'controls.json').read_text())['passed']
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and len(r['rows'])==4
    for row in r['rows']:
        assert row['precision_intersection']
        for result in row['evaluations']:
            assert result==evaluate(result['cell'],result['nodes'],result['end'],result['precision_bits'],p)
    print('PASS original-radius root brackets and enclosed frozen thermal gradients; physical boundary remains separate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
