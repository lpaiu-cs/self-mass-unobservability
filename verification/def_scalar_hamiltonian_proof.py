"""Conditional algebra and saved-endpoint replay; no trajectories are rerun."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import def_scalar_regular_radius as model

ROOT, ld, PI = model.ROOT, model.ld, model.PI
OUT = model.OUT / 'execution-fixed'


def algebra():
    model.prior.symbolic()
    r,h,p,x,xr,m,mp,mpp = sp.symbols('r h p x xr m mp mpp', real=True)
    K=r*r*h*p*p/2
    L=r*(r-h)*h*(x-h*xr/2)**2/4
    R=r*(r+h)*h*(x+h*xr/2)**2/4
    C=1-K/r-2*L/(r-h/2)
    D=1+K/r+2*R/(r+h/2)
    logratio=sp.series(sp.log(C)-sp.log(D),h,0,3).removeO().expand()
    assert sp.simplify(logratio+r*h*(p*p+x*x))==0
    ml,mr=m-h*mp/2+h*h*mpp/8,m+h*mp/2+h*h*mpp/8
    constraint=mr-ml-K*(1-(ml+mr)/r)-L*(1-2*ml/(r-h/2))-R*(1-2*mr/(r+h/2))
    series=sp.series(constraint,h,0,3).removeO()
    assert sp.simplify(series-h*(mp-r*r*(1-2*m/r)*(p*p+x*x)/2))==0
    # Exact dual-cell gradient energy for linear u=r*phi; centre u~r has zero gradient.
    a,b,fa,fb=sp.symbols('a b fa fb', positive=True)
    intercept=a*b*(fa-fb)/(b-a)
    assert sp.simplify(intercept**2*(1/a-1/b)-a*b*(fb-fa)**2/(b-a))==0
    # Heat stress identity at fixed B=a*rho*W and J=a^2*S.
    e,P,Q,v,de,dP,dQ,dv,ell,alpha,T,dphi=sp.symbols('e P Q v de dP dQ dv ell alpha T dphi')
    W2=1/(1-v*v);w=e+P
    E=(e+P*v*v+2*Q*v)*W2
    S=(w*v+Q*(1+v*v))*W2
    stress=(e*v*v+P+2*Q*v)*W2
    differential=lambda F:sum(sp.diff(F,z)*dz for z,dz in [(e,de),(P,dP),(Q,dQ),(v,dv)])
    assert sp.simplify(differential(E)-v*differential(S)-de-v*dQ-(S+Q)*dv)==0
    first_law=-w*(ell+W2*v*dv)-alpha*T*dphi
    dE_fixed_J=-2*v*S*ell+first_law+v*dQ+(S+Q)*dv
    remainder=sp.factor(dE_fixed_J+(E+stress)*ell+alpha*T*dphi)
    assert sp.simplify(remainder-(v*dQ+2*Q*W2*dv+2*Q*v*ell))==0
    assert sp.simplify(remainder.subs({Q:0,dQ:0}))==0
    rest=remainder.subs({v:0,dQ:0,dv:-2*Q*ell/w})
    assert sp.simplify(rest+4*Q*Q*ell/w)==0
    return dict(classification='Proven',passed=True,
        checks=['polarized constraint','interior mass and lapse expansion through h^2',
                'dual-cell gradient energy','fixed baryon and a^2*S perfect-fluid identity',
                'nonzero-heat obstruction at v=0'],
        scope='Conditional algebra only: smooth interior and isentropic fixed-composition matter first law. No global nonlinear PDE bound or nonequilibrium heat Hamiltonian.')


def replay():
    manifest=json.loads((OUT/'manifest.json').read_text())
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,sha in {**manifest['sha256'],**plan['bindings']}.items():
        assert model.digest(ROOT/rel)==sha,rel
    saved=json.loads((OUT/'result.json').read_text());assert saved['passed']
    paths={};drifts=[]
    for file in sorted(OUT.glob('*.npz')):
        kind,n,k=file.stem.split('-');n,k=int(n),int(k);gravity=kind=='gr'
        s=model.Scalar(n,gravity);y=np.load(file)['state'];paths[gravity,n,k]=y
        m0=s.state(model.prior.initial(s,'.05'))['mf'][-1]
        drift=float(abs(s.state(y)['mf'][-1]-m0)/m0)
        row=json.loads(file.with_suffix('.json').read_text())
        assert drift==row['relative_final_mass_drift'];drifts.append(drift)
    assert len(paths)==8
    norm=model.prior.weighted_norm
    errors=[]
    for n in [32,64,128]:
        s=model.Scalar(n,False);profile=ld('.05')*np.sinc(s.r)
        exact=np.array([profile*np.cos(PI*ld('1.7')),-PI*profile*np.sin(PI*ld('1.7'))])
        errors.append(norm(s,paths[False,n,4*n]-exact))
    flat_orders=np.log2(np.array(errors[:-1])/np.array(errors[1:]))
    times=[norm(model.Scalar(64),paths[True,64,k]-paths[True,64,2*k]) for k in [128,256]]
    spaces=[]
    for n in [32,64]:
        fine=paths[True,2*n,8*n]
        spaces.append(norm(model.Scalar(n),paths[True,n,4*n]-(fine[:,::2]+fine[:,1::2])/2))
    time_orders=np.log2(times[0]/times[1]);space_orders=np.log2(spaces[0]/spaces[1])
    for key,value in [('flat_exact_errors',errors),('flat_exact_orders',flat_orders),
                      ('nonlinear_time_orders',time_orders),('nonlinear_space_orders',space_orders)]:
        assert np.array_equal(np.asarray(value,dtype=float),saved[key]),key
    assert max(drifts)==saved['maximum_relative_mass_drift']
    return dict(classification='Counterexample candidate',passed=True,paths=8,
        maximum_relative_mass_drift=max(drifts),flat_orders=flat_orders.astype(float).tolist(),
        time_orders=time_orders.astype(float).tolist(),space_orders=space_orders.astype(float).tolist(),
        scope='Source hashes, eight saved endpoints, constraints, flat exact errors and convergence orders replayed. Intermediate trajectories and their step diagnostics are bound, not independently reconstructed.')


if __name__=='__main__':
    result=dict(algebra=algebra(),regular_operator=model.check(),saved_endpoints=replay(),
                proof_source_sha256=model.digest(Path(__file__)))
    target=model.OUT/'proof-and-replay.json'
    assert not target.exists()
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
