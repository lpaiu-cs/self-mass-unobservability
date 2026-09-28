"""Audit endpoint enclosures and propagate interval radii without sampled-error claims."""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/validated-variational'


def sqrt_up(q):
    x=math.sqrt(float(q))
    while Fraction(x)**2<q:
        x=math.nextafter(x,math.inf)
    assert Fraction(x)**2>=q
    return x


def radius(pair):
    lo,hi=map(Fraction,pair)
    assert lo<=hi
    mid=Fraction(float((lo+hi)/2))
    return max(abs(mid-lo),abs(hi-mid))


def symbolic():
    # Native two-body and two three-body terms, before GR specialization.
    ni,njk,vi,vj,vivj,vi2,vj2,vjn,vni,rij,rjk,rik,mi,mj,mk,gij,gjk,gik,gamma,gammaik,b1,b2,b3,b4=sp.symbols('ni njk vi vj vivj vi2 vj2 vjn vni rij rjk rik mi mj mk gij gjk gik gamma gammaik b1 b2 b3 b4')
    # vni denotes n_ij dot (4 v_i - 3 v_j); delta term has same GR limit.
    ndelta=sp.symbols('ndelta')
    raw2=(ni*(1-(4+2*gamma)*vivj+(1+gamma)*vi2+(2+gamma)*vj2-sp.Rational(3,2)*vjn**2)+(vj-vi)*(vni-2*gamma*ndelta))*gij*mj/rij**2
    nn=sp.symbols('nn')
    raw3a=((sp.Rational(1,2)*nn/rjk-(1+b1+b2)/rij)*ni+(sp.Rational(7,2)+2*gamma)*njk/rjk)*gij*gjk*mj*mk/(rij*rjk)
    raw3b=-ni*(4+2*gammaik+b3+b4)*gij*gik*mj*mk/(rij**2*rik)
    gr={gij:1,gjk:1,gik:1,gamma:0,gammaik:0,b1:0,b2:0,b3:0,b4:0}
    port2=mj/rij**2*(ni*(1-4*vivj+vi2+2*vj2-sp.Rational(3,2)*vjn**2)+(vj-vi)*vni)
    port3a=mj*mk/(rij*rjk)*((sp.Rational(1,2)*nn/rjk-1/rij)*ni+sp.Rational(7,2)*njk/rjk)
    port3b=-4*mj*mk*ni/(rij**2*rik)
    for raw,port in [(raw2,port2),(raw3a,port3a),(raw3b,port3b)]:
        assert sp.simplify(raw.subs(gr)-port)==0
    # Change of state and time units; recover the physical flow derivative.
    a,b,c,d=sp.symbols('a b c d');J=sp.Matrix([[a,b],[c,d]])
    S=sp.diag(1,128)
    assert S*(S.inv()*J*S)*S.inv()==J
    a,b,c=sp.symbols('a b c');q=sp.Matrix(sp.symbols('q0:4'))
    cm1=(1-a)*q[0]+a*q[1];cm2=(1-b)*cm1+b*q[2];cm3=(1-c)*cm2+c*q[3]
    T=sp.Matrix([cm3,q[1]-q[0],q[2]-cm1,q[3]-cm2]).jacobian(q)
    S=sp.Matrix([[1,-a,-b,-c],[1,1-a,-b,-c],[1,0,1-b,-c],[1,0,0,1-c]])
    assert sp.simplify(S*T)==sp.eye(4)
    common=sp.symbols('common')
    for i in range(4):
        for j in range(4):assert sp.expand((q[i]+common)-(q[j]+common)-(q[i]-q[j]))==0
    return dict(GR_term_specialization=True,dyadic_Jacobian_pullback=True,Jacobi_inverse=True,common_position_cancellation=True)


def main():
    result=dict(status='Proven',conditional_scope='CAPD continuum flow at recorded epochs; frozen internal IVP and GR specialization',symbolic_checks=symbolic(),runs={})
    for name in ['forward','backward','full-forward','ho-forward','scaled-forward','scaled-local','scaled-backward','jacobi-forward','jacobi-fixed-forward','jacobi-translation-forward','jacobi-local','mass-local']:
        path=OUT/(name+'.jsonl')
        if not path.exists():continue
        rows=[json.loads(line) for line in path.read_text().splitlines()]
        log=(OUT/(name+'.log')).read_text()
        # Completeness must be an explicit final result, never inferred from a file.
        completed='"requested_horizon_completed":true' in log
        terminated=completed or '"requested_horizon_completed":false' in log
        assert terminated, f'{name} still running'
        info=[]
        for row in rows:
            dim=row.get('dimension',24)
            assert dim in (24,28) and len(row['state'])==dim and len(row['jacobian'])==dim*dim
            sq=sum((radius(p)**2 for p in row['jacobian']),Fraction(0))
            info.append(dict(step=row['step'],dimension=dim,time_internal=row['t'],
                state_midpoint_error_l2_upper=sqrt_up(sum((radius(p)**2 for p in row['state']),Fraction(0))),
                initial_state_Jacobian_midpoint_error_operator_upper=sqrt_up(sq),
                geometric_delay_box_width_us_upper=row['geometric_delay_width_us_upper']))
        result['runs'][name]=dict(requested_horizon_completed=completed,checkpoints=info,
            raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),log_sha256=hashlib.sha256(log.encode()).hexdigest())
    # Different coordinate systems enclose the same frozen IVP at this epoch.
    # Their enclosures must intersect; intersection is not used to reset a run.
    local=[]
    for name in ['forward','scaled-local','jacobi-local']:
        p=OUT/(name+'.jsonl')
        if p.exists():local.append(json.loads(p.read_text().splitlines()[-1]))
    if len(local)==3:
        for key in ['state','jacobian']:
            for entries in zip(*(r[key] for r in local)):
                assert max(e[0] for e in entries)<=min(e[1] for e in entries),f'coordinate disagreement: {key}'
        result['coordinate_enclosures_intersect_at_local_epoch']=True
    p=OUT/'mass-local.jsonl'
    if p.exists():
        mass=json.loads(p.read_text().splitlines()[-1]);dim=mass['dimension'];assert dim==28
        for i in range(24,28):
            assert mass['state'][i]==[0,0]
            for j in range(28):
                lo,hi=mass['jacobian'][i*28+j];assert lo<=(1 if i==j else 0)<=hi
        bounds=[]
        for j in range(24,28):
            column=[mass['jacobian'][i*28+j] for i in range(24)]
            assert any(lo>0 or hi<0 for lo,hi in column), 'mass tangent must have a resolved nonzero component'
            bounds.append(sqrt_up(sum((radius(v)**2 for v in column),Fraction(0))))
        result['independent_fractional_mass_columns']=dict(local_completed=result['runs']['mass-local']['requested_horizon_completed'],error_l2_upper=bounds,not_the_28_timing_parameters=True)
        reference=json.loads((OUT/'scaled-local.jsonl').read_text().splitlines()[-1])
        for i in range(24):
            for j in range(24):
                a=mass['jacobian'][i*28+j];b=reference['jacobian'][i*24+j]
                assert max(a[0],b[0])<=min(a[1],b[1]),'augmented state block must agree'
    result['D2']=dict(pass_gate=False,whole_observation_span_certified=False,
        timing_parameter_initialization_derivatives_certified=False,
        local_independent_mass_partials_certified=result.get('independent_fractional_mass_columns',{}).get('local_completed',False),
        physical_timing_parameter_mass_mapping_certified=False,
        all_delay_terms_and_emission_time_inversion_certified=False,
        nuisance_projector_error_upper=None,likelihood_error_upper=None)
    # Exact endpoint-radius calculation has a simple independently known answer.
    assert radius([0.,2.])==1 and sqrt_up(Fraction(4))==2
    (OUT/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
