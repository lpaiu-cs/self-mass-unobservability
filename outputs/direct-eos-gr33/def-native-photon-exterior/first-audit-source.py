"""Integral conservation, null rays, stress moments and metric reconstruction."""
from pathlib import Path
import json
import signal
import time
import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid
import def_native_photon_exterior as p


def moments(e,r,t):
    mu,w,delay=e.rays(r,128,48)
    lum=p.history(t-delay)
    integ=p.history(t-delay,True)@w
    _,n,b,c,v=e.metric(r)
    factor=1/(4*np.pi*r*r*n*n)
    energy=factor*((lum/mu)@w)
    radial=factor*((lum*mu)@w)
    flux=factor*(lum@w)
    return energy,radial,flux,integ,(n,b,c,v)


def main():
    assert not (p.OUT/'audit.json').exists()
    p.write(p.OUT/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Independently integrate redshifted photon energy, check local energy and radial momentum conservation, flat exact ray travel time, analytic steady hemispheric moments, and reconstruct the radial metric perturbations from the solved scalar field.',
        gates=dict(flat_delay_absolute=1e-10,steady_moments_relative=1e-7,
            energy_inventory_relative=1e-5,local_conservation_relative=1e-5,
            exact_debit_credit_absolute=1e-10,physical_cutoff_relative=.01),
        budget=dict(hard_seconds=30,native_calls=0,new_evolution_paths=0),
        bindings={str(x.relative_to(p.ROOT)):p.digest(x) for x in
                  [Path(__file__),Path(p.__file__),p.OUT/'result.json',p.OUT/'angular.npz',p.OUT/'thin.npz']}))
    signal.alarm(30);start=time.monotonic()
    result=json.loads((p.OUT/'result.json').read_text());assert result['passed']
    e=p.Exterior()
    points=np.array([1.002,1.04,1.3,2.,4.,6.])
    flat=p.Exterior(flat=True)
    mu,w,delay=flat.rays(points,64,48)
    x,_=np.polynomial.legendre.leggauss(64);mu0=(x+1)/2
    expected=np.sqrt(points[:,None]**2-1+mu0**2)-mu0
    flat_error=float(max(abs(delay-expected).ravel()))
    assert flat_error<1e-10
    r=np.r_[1.,np.geomspace(1.0001,9,120)]
    mu,w,_=e.rays(r,128,48)
    n=e.metric(r)[1];minimum=np.sqrt(np.maximum(0,1-(n/(e.Nb*r))**2))
    expectedE=2/(1+minimum)
    expectedP=2*(1+minimum+minimum**2)/(3*(1+minimum))
    steady_error=float(max(np.max(abs((1/mu)@w/expectedE-1)),np.max(abs(mu@w/expectedP-1))))
    assert steady_error<1e-7,steady_error
    # A different radial integration of the local E gives stored photon energy.
    grid=1+8*np.linspace(0,1,1601)**2
    mu,w,delay=e.rays(grid,128,48)
    _,n,b,c,_=e.metric(grid)
    inventory=[]
    for t in [.5,1.5,3.,6.]:
        local=(p.history(t-delay)/mu)@w/c
        measured=simpson(local,x=grid)
        emitted=float(p.history(t,True))
        crossed=float(p.history(t-delay[-1],True)@w)
        expected=emitted-crossed
        error=abs(measured-expected)/max(expected,1e-30)
        assert error<1e-5,(t,error)
        inventory.append(dict(time=t,emitted=emitted,stored=measured,escaped=crossed,relative_error=error))
    # Numerical derivatives test conservation, rather than just repeating its
    # cumulative-energy identity. Orthornormal E,P,F use c=1 in these units.
    dr=1e-5;dt=1e-5;local_errors=[]
    for t in [.5,1.5,3.,6.]:
        E,P,F,H,(n,b,c,v)=moments(e,points,t)
        Em,Pm,Fm,Hm,_=moments(e,points-dr,t)
        Ep,Pp,Fp,Hp,_=moments(e,points+dr,t)
        Et,Pt,Ft,Ht,_=moments(e,points,t+dt)
        Eb,Pb,Fb,Hb,_=moments(e,points,t-dt)
        mass=e.metric(points)[0]
        nuprime=mass/(points**2*b)+points*v*v/2
        radial=(Pp-Pm)/(2*dr)+nuprime*(E+P)+(3*P-E)/points+(Ft-Fb)/(2*dt*c)
        energy=(Hp-Hm)/(2*dr)+4*np.pi*points**2*n/np.sqrt(b)*E
        temporal=(Ht-Hb)/(2*dt)-4*np.pi*points**2*n*n*F
        norms=[float(max(abs(radial))/max(np.max(abs(E/points)),1e-30)),
               float(max(abs(energy))/max(np.max(abs(4*np.pi*points**2*n/np.sqrt(b)*E)),1e-30)),
               float(max(abs(temporal))/max(np.max(abs(4*np.pi*points**2*n*n*F)),1e-30))]
        assert max(norms)<1e-5,(t,norms)
        local_errors.append(dict(time=t,radial_momentum=norms[0],radial_energy=norms[1],time_energy=norms[2]))
    saved=np.load(p.OUT/'angular.npz');r=saved['r'];psi=saved['field']
    E,P,F,H,(n,b,c,v)=moments(e,r,6.)
    J=-H*np.sqrt(b)/n
    dphi=e.scalar_scale*psi
    dphi_r=np.gradient(dphi,r,edge_order=2)
    dm=r*r*b*v*dphi+e.q*J
    dnuprime=dm/(r*r*b*b)+r*v*dphi_r+4*np.pi*r*e.q*P/b
    primitive=cumulative_trapezoid(dnuprime,r,initial=0)
    dnu=primitive-primitive[-1]
    dlambda=dm/(r*b)
    np.savez(p.OUT/'metric.npz',r=r,dphi=dphi,dphi_time=e.scalar_scale*saved['velocity'],
        delta_m_over_R=dm,delta_logN=dnu,delta_lambda=dlambda,delta_logN_prime=dnuprime,
        photon_E=e.q*E,photon_Pr=e.q*P,photon_F=e.q*F,J=e.q*J)
    # Compare physical output units as well as the producer's scaled solutions.
    a=np.load(p.OUT/'angular.npz')['traces'];d=np.load(p.OUT/'thin.npz')['traces']
    ea,ed=e,p.Exterior('thin')
    # Same physical coordinate time; extraction points remain 2 and 3 surface
    # radii in each model, so this is a cutoff robustness check, not fixed-r data.
    physical=np.column_stack([np.interp(a[:,0]*ea.R/ed.R,d[:,0],d[:,i])*ed.scalar_scale for i in range(1,4)])
    reference=a[:,1:]*ea.scalar_scale
    cutoff=(np.max(abs(physical-reference),axis=0)/np.max(abs(reference),axis=0)).tolist()
    assert max(cutoff)<.01,cutoff
    total=result['seconds']+json.loads((p.OUT/'pilot.json').read_text())['total_seconds']+time.monotonic()-start
    assert total<240
    audit=dict(classification='Counterexample candidate',passed=True,
        flat_exact_ray_delay=flat_error,steady_hemisphere_relative=steady_error,
        radial_energy_integrals=inventory,local_conservation=local_errors,
        physical_cutoff_comparison_relative=cutoff,
        metric_maximum_delta_m_over_R=float(max(abs(dm))),metric_maximum_delta_logN=float(max(abs(dnu))),
        geometric_mass_at_r9=float(dm[-1]),
        energy_bookkeeping='At t=6 all pulse photons remain inside r=9. Body debit plus redshifted photon inventory is zero. This enclosed ADM-type balance is not the retarded Bondi mass loss; no eternal steady radiation bath is introduced.',
        boundary_force_use='The saved R*dphi_prime(surface) is the additive exterior forcing for Dirichlet data zero. The full response requires the homogeneous exterior boundary operator acting on the interior-evolved surface value.',
        remaining_gas_pressure_fraction=ea.Pgas/ea.Prad,
        full_physical_boundary_matched=False,final_stellar_charge_solved=False,full_goal_complete=False,
        seconds=time.monotonic()-start,total_compute_seconds=total)
    p.write(p.OUT/'audit.json',audit);signal.alarm(0)
    print(json.dumps(audit),flush=True)


if __name__=='__main__':main()
