"""Check the low-momentum boundary before accepting a fitted memory time."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import def_electron_correlated_response as model


def main():
    h=model.h;out=model.base.OUT/'infrared-boundary';assert not out.exists();out.mkdir()
    source=model.base.OUT/'correlated-regular'
    paths=[Path(__file__),Path(model.__file__),source/'plan.json',source/'response.npz',source/'result.json']
    h.write(out/'plan.json',dict(classification='Proven',
        sources={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        claim='Determine whether the first memory moment exists for the declared elastic off-Fermi extension. Do not refine an actually divergent integral.',
        budget=dict(hard_seconds=30,new_kinetic_production=0,new_stellar_time_steps=0,new_native_calls=0)))
    p,q,a,b=sp.symbols('p q a b',positive=True)
    integrand=q*(-sp.expm1(-b*p*p*q) if hasattr(sp,'expm1') else 1-sp.exp(-b*p*p*q))/(2*(q+a/(p*p))**2)
    leading=sp.limit(integrand/p**6,p,0,dir='+')
    assert sp.simplify(leading-b*q*q/(2*a*a))==0
    angular_leading=sp.integrate(leading,(q,0,1));assert angular_leading==b/(6*a*a)
    # p^4 dp is the nonrelativistic endpoint of the relativistic current
    # measure. The collision rate is p^-3 times this angular p^6 factor.
    eps=sp.symbols('eps',positive=True)
    assert sp.limit(sp.integrate(1/p**2,(p,eps,1)),eps,0,dir='+')==sp.oo
    fractional=sp.integrate(p/(1+p**3),(p,0,sp.oo))
    assert sp.simplify(fractional-2*sp.pi/(3*sp.sqrt(3)))==0
    proof=dict(classification='Proven',passed=True,
        angular_limit='Lambda/p^6 -> b/(6*a^2) for u=a/p^2,w=b*p^2, v/c=O(p).',
        collision_limit='nu_ei(p)=C*p^3+o(p^3), C>0.',
        measures='W dE is proportional to p^4 dp near p=0 at any positive T. The DC variance density W*tau*(z-mean)^2 dE scales as p dp and is integrable; its first time moment with another tau scales as dp/p^2 and diverges when z(0)!=mean.',
        response='The same elastic extension has a finite frequency response but a nonanalytic leading correction proportional to -s^(2/3); integral x/(1+x^3) dx = 2*pi/(3*sqrt(3)). A finite tau_eff cannot represent its derivative at zero.',
        scope='Only the declared off-Fermi elastic continuation of the fitted ion potential. This is not a theorem that the actual finite-temperature plasma has an infrared divergence: inelastic and electron-electron collisions were omitted.')
    h.write(out/'symbolic.json',proof)
    r=np.load(source/'response.npz');d,data=model.base.inputs();idx=r['indices'];T=data['T'][idx]
    t=model.k*T/(model.m_e*model.c**2);pdim=np.sqrt(t[:,None]*r['x']*(2+t[:,None]*r['x']))
    slopes=-np.log(r['tau'][:,1]/r['tau'][:,0])/np.log(pdim[:,1]/pdim[:,0])
    separation=abs(-r['eta']-r['moments'][1]/r['moments'][0])
    assert np.min(separation)>0
    result=dict(classification='Counterexample candidate',saved_cohort_cells=len(idx),
        first_two_node_collision_power_range=[float(slopes.min()),float(slopes.max())],
        zero_momentum_thermal_weight_separation_range=[float(separation.min()),float(separation.max())],
        DC_integral_compatible=True,correlated_tau_eff_accepted=False,frequency_controls_scaled_by_that_tau_accepted=False,
        decision='Withdraw the finite tau_eff interpretation of correlated-regular/result.json; retain it as a failed truncated-integral result. Static ei+ee compatibility is independent and remains accepted. Do not increase energy resolution to cure an infinite first moment.',
        next='Construct a finite-temperature electron-electron/inelastic collision operator with the correct conservation and detailed balance before dynamic use. A published DC ee rate is insufficient to choose that operator.',
        full_dynamic_charge_solved=False)
    h.write(out/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
