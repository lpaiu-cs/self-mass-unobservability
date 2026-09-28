"""Initial insulating-boundary compatibility and exact frozen-profile endpoints."""
from fractions import Fraction as F
from functools import reduce
from operator import mul
import json,subprocess,sys
import numpy as np
import sympy as sp
from mpmath import iv
from gr_logarithmic_gauss_rule import I,low,high
import gr_subcell_heat_drive_coordinates as heat

ROOT=heat.ROOT;OUT=heat.OUT.parent/'gr-heat-boundary-compatibility';sha=heat.sha;ld=np.longdouble


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def exact(x):return F(*ld(x).as_integer_ratio())


def endpoint_rule(number,end):
    points=heat.original.derivative_rule(number)[0];s=list(map(exact,points));e=F(end)
    values=[];derivatives=[]
    for j in range(number):
        q=reduce(mul,((e-s[k])/(s[j]-s[k]) for k in range(number) if k!=j),F(1))
        values.append(q);derivatives.append(q*sum((1/(e-s[k]) for k in range(number) if k!=j),F(0)))
    for k in range(number):
        assert sum((a*x**k for a,x in zip(values,s)),F(0))==e**k
        target=F(0) if k==0 else k*e**(k-1)
        assert sum((a*x**k for a,x in zip(derivatives,s)),F(0))==target
    assert sum(derivatives,F(0))==0
    return s,values,derivatives


def symbolic():
    w,h,N,tau,K,T,a,c=sp.symbols('w h N tau K T a c',positive=True)
    R,QF,vdot,qd,gradient=sp.symbols('R QF vdot qdot gradient',real=True)
    # At initial v=Q=0, a material boundary with Q identically zero has qdot=0.
    momentum=qd+w*vdot-R;transport=tau*qd+w*h*vdot-N*QF
    assert sp.simplify(transport.subs({qd:0,vdot:R/w})-(h*R-N*QF))==0
    assert sp.diff(h*R-N*QF,tau)==0
    qf=-K*T*gradient/(a*c)
    assert sp.simplify((h*R-N*QF).subs({R:0,QF:qf})-N*K*T*gradient/(a*c))==0
    return dict(classification='Proven',passed=True,
        assumptions='Classical differentiable material boundary, initial material v=0 and Q=0, no additional boundary heat channel, and Q=0 along that material worldtube for a time interval. Consequently Q_dot=0 at the initial event. Use the same complete heat-only quadratic law and momentum equation.',
        general_condition='If Q_dot+w*v_dot=R, the full initial heat equation at Q=0 is tau*Q_dot+w*h*v_dot=N*Q_F. Compatibility therefore requires h*R=N*Q_F.',
        hydrostatic_condition='At an initially hydrostatic boundary R=0, compatibility requires v_dot=0 and Q_F=0. For positive K,T,a,c this means partial_r ln(T*N)=0 at that boundary.',
        cannot_retune_tau='The compatibility condition is independent of tau because Q_dot=0 at the stipulated material boundary. Changing a positive finite relaxation time alone cannot repair a nonzero thermal force.',
        alternatives='A nonzero initial heat time derivative through a pressure-bearing atmosphere/radiation boundary, a consistently non-hydrostatic initial momentum force, or a separately constructed compatible thermal profile changes an explicit premise. None is silently imposed in this audit. A weak initial layer is not ruled out by this classical compatibility condition.',
        scope='Necessary initial-corner condition, not existence of a compatible atmosphere or a finite GR trajectory. A frozen interpolation endpoint test does not certify the physical stellar boundary.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_heat_boundary_compatibility.py',ROOT/'verification/gr_logarithmic_gauss_rule.py',
        ROOT/'verification/interval_records.py',
        heat.OUT/'manifest.json',heat.OUT/'plan.json',heat.OUT/'candidate.py',
        heat.original.closure.OUT/'symbolic.json',heat.original.closure.old.OUT/'symbolic.json']
    prior=json.loads((heat.OUT/'plan.json').read_text())
    save('plan.json',dict(classification='Counterexample candidate',
        checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        bindings=dict(prior['bindings'],**{p.relative_to(ROOT).as_posix():sha(p) for p in files}),
        endpoints=[dict(cell=0,end=1,label='outer surface'),dict(cell=5734,end=0,label='regular centre')],nodes=[8,16],bits=[192,256],
        model='Use exact dyadic stored Gauss parameters, radii, log temperatures and pressure-mode EOS outputs. Define nodal H=u+P/rho and the same frozen C_X*c^2. Interpolate r and log(T)-log(C+H) on those parameters. The common lapse normalization cancels exactly from its derivative. Evaluate the endpoint derivative with exact rational Lagrange weights and outward logarithm intervals.',
        controls='Prove all represented polynomial endpoint value/derivative identities in exact rational arithmetic, then evaluate the same frozen physical numbers at two interval precisions and require interval overlap. Report the endpoint radius defect relative to the original face and report every nonzero gradient without a fitted acceptance tolerance.',
        scope='Declared finite interpolants and the separate conditional boundary theorem. No uncertainty bound on the native EOS, nodal root error, Gauss-node error, true continuous profile, reference TOV residual or physical atmosphere. Centre results diagnose the unconstrained interpolant, not a physical singularity.'))
    save('symbolic.json',symbolic())


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def evaluate(cell,number,end,precision):
    iv.prec=precision
    reference=json.loads((heat.original.metric.OUT/'plan.json').read_text());path=ROOT/reference['cell_sources'][str(cell)]
    data=heat.original.metric.node_file(str(path.with_name(path.stem+f'-nodes-{number}.npz')))
    j=list(data['cells']).index(cell);order=np.argsort(data['radius_cm'][j]);eos=data['eos'][j,order]
    radius=list(map(exact,data['radius_cm'][j,order]));temperature=list(map(exact,data['lnT'][j,order]))
    c=ld(heat.original.metric.g.c.gr.C)*100;C=exact(ld(data['C_X'][j])*c*c)
    enthalpy=[C+exact(row[2])+exact(row[1])/exact(row[0]) for row in eos];assert min(enthalpy)>0
    points,V,D=endpoint_rule(number,end)
    radius_at=sum((x*y for x,y in zip(V,radius)),F(0));rs=sum((d*(r-radius[0]) for d,r in zip(D,radius)),F(0));assert rs>0
    phi=[I(t-temperature[0])-iv.log(I(H/enthalpy[0])) for t,H in zip(temperature,enthalpy)]
    derivative=sum((I(d)*p for d,p in zip(D,phi)),I(0))/I(rs)
    lo,hi=low(derivative),high(derivative);assert lo<=hi
    grid=np.load(heat.original.metric.g.OUT/'gr-increment-structure/path-4.npz')
    face=exact(ld(grid['radius_m'][0 if end==1 else 5735])*100)
    return dict(cell=cell,nodes=number,end=end,precision_bits=precision,
        gradient_per_cm=[str(lo),str(hi)],zero_excluded=lo>0 or hi<0,
        endpoint_radius_cm=str(radius_at),original_face_radius_cm=str(face),radius_endpoint_defect_cm=str(radius_at-face),
        radius_parameter_derivative_cm=str(rs),exact_polynomial_controls_passed=True)


def run():
    p=bindings();assert not (OUT/'result.json').exists();rows=[]
    for target in p['endpoints']:
        for n in p['nodes']:
            pair=[evaluate(target['cell'],n,target['end'],b) for b in p['bits']]
            a,b=[list(map(F,r['gradient_per_cm'])) for r in pair]
            assert max(a[0],b[0])<=min(a[1],b[1])
            rows.append(dict(label=target['label'],evaluations=pair,precision_intersection=True))
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        physical_boundary_incompatibility_certified=False,native_or_continuous_derivative_error_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())==symbolic()
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and len(r['rows'])==4
    for row in r['rows']:
        assert row['precision_intersection']
        for a in row['evaluations']:
            assert a==evaluate(a['cell'],a['nodes'],a['end'],a['precision_bits'])
    print('PASS conditional insulating-boundary identity and frozen-profile endpoint intervals; inspect gradient and radius defects',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
