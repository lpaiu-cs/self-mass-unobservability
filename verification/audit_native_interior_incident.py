"""Independent interior source and selected higher-field enclosure.

Counterexample candidate. The source interpolation and exact static operator
are conditional inputs, not an enclosure of physical EOS or background error.
"""
from pathlib import Path
import json, resource, time
import numpy as np
import mpmath as mp
from scipy.integrate import quad
from numpy.polynomial import legendre as leg
import couple_native_interior_incident as run


def box():
    errors=[];C=run.C;eta=run.infinity.incident.ETA;t=np.linspace(0,1,5);clock=np.linspace(0,1,9)
    p=np.r_[np.zeros(4),256*np.array([1.,-4.,6.,-4.,1.])];D=.5;v=np.minimum(clock/D,1.)
    A=D*np.polynomial.polynomial.polyval(v,np.polynomial.polynomial.polyint(p))
    B=D*D*np.polynomial.polynomial.polyval(v,np.polynomial.polynomial.polyint(np.r_[0.,p]))
    expected=(clock*A-B)/4
    for q in [4,8]:
        gx,_=leg.leggauss(q)
        z=dict(t=t,xfaces=np.array([-C,0.]),inverse=np.linalg.inv(leg.legvander(gx,q-1))[None],
            volume_density=np.full((len(t),q),1/C**2),outer_zeta=np.zeros(len(t)),external_zeta=np.zeros(len(t)),r0=1/eta)
        actual,boundary,external=run.integrate(z,clock)
        errors.append(float(np.max(abs(actual-expected))/max(abs(expected))))
        assert not np.any(boundary) and not np.any(external)
    assert max(errors)<1e-11,errors
    return errors


def flat(u):
    """Direct unreduced wave source on an interior domain, c=1."""
    D=.5;R=2.
    def pulse(t,x):
        s=(t+x)/D
        if not 0<s<1:return 0.,0.,0.
        return (256*s**4*(1-s)**4,1024*s**3*(1-s)**3*(1-2*s)/D,
                1024*s*s*(1-s)**2*(3-14*s+14*s*s)/D**2)
    def integral(fn,a,b):return quad(fn,a,b,epsabs=3e-12,epsrel=3e-12,limit=100)[0] if b>a else 0.
    def source(t,x,reduced):
        U,Ux,Uxx=pulse(t,x);Z=t*(1+x/R);Zx=t/R;Zt=1+x/R
        return -Zx*U/(R+x) if reduced else 2*Z*Uxx+Zx*(Ux-U/(R+x))+Zt*Ux
    cuts=np.unique(np.clip([-u/2,(D-u)/2,0.],-u/2,0.))
    def outer(x,reduced):return integral(lambda t:source(t,x,reduced),-x,min(D-x,u+x))
    direct=sum(integral(lambda x:outer(x,False),a,b) for a,b in zip(cuts[:-1],cuts[1:]))/2
    bulk=sum(integral(lambda x:outer(x,True),a,b) for a,b in zip(cuts[:-1],cuts[1:]))/2
    boundary=integral(lambda t:t*pulse(t,0.)[1],0,min(u,D))/2
    error=abs(direct-bulk-boundary)/max(abs(direct),1e-100);assert error<1e-10,(u,error)
    return dict(u=u,direct=direct,volume=bulk,boundary=boundary,relative=error)


def enclosure():
    run.initialize();z=np.load(run.OUT/'operator-128-g8.npz');r=z['radius']
    model=run.wave.base.Response();bg=model.coeff(r);a=bg['lapse']*np.sqrt(bg['b'])
    prior=run.read(run.exterior.OUT/'bound-result.json');assert prior['analytic_envelope_verified']
    frozen=run.read(run.infinity.ROOT/'def-native-global-scalar-closure/bound.json')
    mp.iv.dps=50;iv=mp.iv;I=iv.mpf
    def B(x):
        n,d=float(x).as_integer_ratio();return I(n)/I(d)
    def up(x):return float(np.nextafter(float(x.b),np.inf))
    inv=z['inverse'];q=inv.shape[-1];n=len(inv);dt=min(np.diff(z['t']))
    maxima=[np.max(abs(z[k]),axis=0).reshape(n,q) for k in ['zeta','volume_density','zeta_prime']]
    sums=[I(0),I(0),I(0)]
    for j in range(n):
        weights=[sum((B(abs(v)) for v in inv[j,:,k]),I(0)) for k in range(q)]
        width=B(z['xfaces'][j+1])-B(z['xfaces'][j])
        for h in range(3):
            sums[h]+=width*sum((weights[k]*B(maxima[h][j,k])*(B(a[j*q+k]) if h==2 else 1) for k in range(q)),I(0))
    Z,H,X=sums;cc=B(run.C);T=B(z['t'][-1]);D=T/2;M=B(z['M']);U0=B(run.infinity.incident.ETA)*B(z['r0'])
    Ut=U0/D*I(2048)/70;Ux=Ut/cc;Uxx=U0/(cc*D)**2*I(2048)/5
    Zt=2*Z/B(dt) # Whole declared time interpolation, with no derivative sampling.
    def source_norm(U,Ux,Uxx,Ut):return cc*T/(2*M)*(2*Z*Uxx+H*U+X*Ux+Zt*Ut/cc**2)
    primary=source_norm(U0,Ux,Uxx,Ut)
    b=prior['Born_source'];born=source_norm(B(b['U']),B(b['Ux']),B(b['Uxx']),B(b['Ut']))
    eta=B(frozen['global_potential_contraction']);assert up(eta)<1
    potential=eta/(1-eta)*(primary+born)
    e0=B(max(np.load(run.infinity.BACKGROUND)['epsilon']));total=(born+potential)/(1-e0)
    point=run.read(run.OUT/'result.json');target=B(abs(point['endpoint_previous_body']))
    return dict(classification='Counterexample candidate',analytic_envelope_verified=up(total/target)<.01,
        primary_free_field_norm=up(primary),finite_Born_input_interaction_bound=up(born),
        additional_static_potential_bound=up(potential),selected_higher_terms_bound=up(total),
        bound_over_stored_body=up(total/target),global_potential_contraction=up(eta),
        endpoint_with_selected_terms_interval=[float(np.nextafter(point['endpoint_combined_selected']-up(total),-np.inf)),
                                             float(np.nextafter(point['endpoint_combined_selected']+up(total),np.inf))],
        scope='Only the interior metric acting on the declared continuous primary plus finite Born source and their additional initial-static-potential returns. The domain boundary is matched. Not the full exterior/interior reciprocal Einstein problem or an error enclosure of the previous body, source interpolation or actual supplied field.',
        full_goal_complete=False)


def main():
    out=run.OUT;assert not (out/'audit-plan.json').exists()
    files=[Path(__file__),Path(run.__file__),out/'plan.json',out/'result.json',out/'operator-128-g8.npz',
           out/'coefficients.json',run.exterior.OUT/'bound-result.json',run.infinity.BACKGROUND,
           run.infinity.ROOT/'def-native-global-scalar-closure/bound.json']
    run.write(out/'audit-plan.json',dict(classification='Counterexample candidate',budget_seconds=60,
        claim='Independently test the characteristic pulse integration and interior boundary sign, then bound the missing finite-Born and further fixed-potential terms for this current interior operator.',
        method='Exact box/pulse antiderivative at4/8 source degrees; independent unreduced flat wave source versus volume plus boundary.50-digit outward intervals use exact saved binary nodal values and inverse-Legendre coefficients; all time slopes are bounded by2*sup/dtmin.',
        gate='Original independent1e-10 checks. Separately require selected higher terms below1percent of the stored body; this does not certify body, interpolation, EOS or full mixed constraints.',
        reuse='No new fluid/EOS states, physical mesh, source times or potential iterations. Reuse previous continuous finite-Born source and same initial-static global coefficient bound; rebuild norms for the actual current interior source.',
        stop='Preserve a failed bound or audit; do not automatically increase quadrature, horizon, sampling or calculation budget.',
        bindings={str(p):run.sha(p) for p in files}))
    start=time.monotonic();cpu=time.process_time();run.infinity.incident.native.deadline(60)
    rows=[flat(u) for u in [.2,.5,1.]];checks=box();bound=enclosure()
    run.write(out/'higher-bound.json',bound)
    assert bound['analytic_envelope_verified'],bound
    data=run.read(out/'result.json');coeff=run.read(out/'coefficients.json')
    assert data['passed'] and max(r['interface_relative'] for r in coeff['rows'])<1e-12
    assert max(r['ray_invariant'] for r in coeff['rows'])<1e-10
    for p,h in run.read(out/'plan.json')['bindings'].items():assert run.sha(p)==h,p
    result=dict(classification='Counterexample candidate',passed=True,flat_unreduced_source=rows,
        exact_box_pulse_relative=checks,actual_interface_matched=True,current_source_bindings_verified=True,
        selected_Born_and_potential_terms_bounded=True,full_mixed_GR_closed=False,full_goal_complete=False,
        seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    run.write(out/'audit.json',result);print(json.dumps(dict(audit=result,bound=bound)),flush=True)


if __name__=='__main__':
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));start=time.monotonic();cpu=time.process_time();error=None
    assert not (run.OUT/'audit-receipt.json').exists()
    try:main()
    except Exception as exc:error=repr(exc);raise
    finally:run.write(run.OUT/'audit-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=run.sha(__file__)))
