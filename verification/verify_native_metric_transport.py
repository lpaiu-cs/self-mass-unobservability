"""Independent quadrature of the lapse constraint and actual boundary work."""
from pathlib import Path
import json
import signal
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline
import def_native_metric_transport as task

OUT=task.OUT


def main():
    assert not (OUT/'audit.json').exists();start=time.monotonic();signal.alarm(40)
    task.write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Independently integrate the same Einstein lapse variation with8-point volume and96-angle exterior quadrature; verify moving energy/momentum support at the actual final surface.',
        gates=dict(lapse_relative=.02,moving_mass_jump_relative=1e-7,closure=1e-8),
        budget_seconds=40,new_evolution_paths=0,
        boundary='An energy/momentum ledger with prescribed pressure does not specify the gravitating exterior support.',
        bindings={str(p.relative_to(task.old.ROOT)):task.old.photons.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'result.json',OUT/'p4-64.npz']}))
    m=task.Model();d=np.load(OUT/'p4-64.npz');state=tuple(d[k] for k in ['q','v','e','d'])
    m.history_t=d['emission_times'].tolist();m.history_e=d['emission_energy_deviation'].tolist();m.history_d=d['emission_flux_deviation'].tolist()
    m.history_z=d['surface'][:,0].tolist();m.history_zd=d['surface'][:,2].tolist();t=m.history_t[-1]
    m.closure(state,t);reference=m.last_lapse_native.copy();surface_reference=float(m.last_lapse_surface)
    cuts=np.unique(np.r_[m.cells,m.edges,m.native,m.bg.x]);cuts=cuts[cuts<=1]
    gx,gw=np.polynomial.legendre.leggauss(8)
    r=(cuts[:-1,None]+np.diff(cuts)[:,None]*(gx+1)/2).ravel();weights=(np.diff(cuts)[:,None]*gw/2).ravel()
    p=m.bg.sample(r);b=1-2*p['m']/r;A=np.exp(-8*p['phi']**2);alpha=-4*p['phi'];V,D=m.evaluation(r)
    q,v,e,f=state;z=V[0]@q;zp=D[0]@q;psi=V[1]@q;psip=D[1]@q
    _,loss0,J0=m.source_values(p,m.f0);_,losse,Je=m.source_values(p,e)
    loss=t*loss0+losse;J=t*J0+Je;ad=np.interp(r,m.native,m.thermo[:,4])
    nuprime=p['m']/(r*r*b)+4*np.pi*r*A*p['p']/b+r*p['v']**2/2
    dlz=r*r*p['v']**2/2-4*np.pi*r*r*A*p['p']/b-p['m']/(r*b)
    dp=-p['gamma']*p['p']*(r*zp+(3+dlz+3*alpha*r*p['v'])*z+(r*p['v']+3*alpha)*psi+J/(r*b))-ad*loss
    dp+=r*z*(p['e']+p['p'])*(nuprime+alpha*p['v'])
    dm=r*r*b*p['v']*psi-4*np.pi*r**3*A*(p['e']+p['p'])*z+J
    derivative=(1+8*np.pi*r*r*A*p['p'])*dm/(r*r*b*b)+4*np.pi*r*A*(dp+4*alpha*p['p']*psi)/b+r*p['v']*psip
    cumulative=np.r_[np.longdouble(0),np.cumsum(weights*derivative,dtype=np.longdouble)]
    bulk8=surface_reference-(cumulative[-1]-cumulative[np.searchsorted(r,m.native)])
    bulk_error=float(max(abs(bulk8-reference))/max(abs(reference)))
    # Rebuild exterior quadrature independently; same saved physical history.
    ec=m.cells[m.cells>=1];re=(ec[:-1,None]+np.diff(ec)[:,None]*(gx+1)/2).ravel();we=(np.diff(ec)[:,None]*gw/2).ravel()
    rays=task.moving.Rays(m.bg,re,96,64);u=t-rays.delay;positive=u>0;clipped=np.maximum(u,0)
    heat=CubicHermiteSpline(m.history_t,m.history_e,m.history_d);move=CubicHermiteSpline(m.history_t,m.history_z,m.history_zd)
    H=(m.f0[-1]*clipped+heat(clipped))*positive;L=(m.f0[-1]+heat(clipped,1))*positive
    hm,_,pm=rays.moments(t,move);factor=task.old.G/(task.old.C**4*m.bg.R)
    J=-factor*(H@rays.weights+m.f0[-1]*hm)*np.sqrt(rays.b)/rays.n
    P=factor*((L*rays.mu)@rays.weights+m.f0[-1]*pm)/(4*np.pi*re*re*rays.n**2)
    ps=m.bg.sample(np.array([1.]));Vs=m.surfaceV[1];Ve,_=m.evaluation(re)
    surface8=ps['v'][0]*(Vs@q)[0]-2*np.sum(we*rays.v/rays.b*(Ve[1]@q))
    surface8-=np.sum(we*(J/(re*re*rays.b**2)+4*np.pi*re*P/rays.b))
    surface_error=float(abs(surface8-surface_reference)/abs(surface8))
    full8=bulk8-surface_reference+surface8
    total_error=float(max(abs(full8-reference))/max(abs(full8)))
    # Test the surface J jump from retarded moving rays, without subtracting
    # the much larger common emitted-energy debit.
    surface_rays=task.moving.Rays(m.bg,np.array([1.]),48,32)
    h=surface_rays.moments(t,move)[0][0];bs=1-2*ps['m'][0];Ns=ps['N'][0]
    dJ=-factor*m.f0[-1]*h*np.sqrt(bs)/Ns;zeta=d['surface'][-1,0]
    env=np.load(task.old.prior.OUT/'final-envelope.npz');geo=task.old.G*m.bg.R**2/task.old.C**4;A4=np.exp(-8*ps['phi'][0]**2)
    Pr=float(env['Prad'][-1])*geo;Pg=float(env['Pgas'][-1])*geo
    radiation_expected=-4*np.pi*A4*(3*Pr+Pr)*zeta
    mass_error=float(abs(dJ-radiation_expected)/abs(radiation_expected))
    support_mass=4*np.pi*A4*Pg*zeta
    expected=-4*np.pi*A4*(Pr+Pg+3*Pr)*zeta+support_mass
    assert abs(expected-radiation_expected)/abs(radiation_expected)<1e-14
    rows=[json.loads((OUT/f'p4-{n}.json').read_text()) for n in [8,32,64]]
    maxclosure=max(row['max_closure_relative'] for row in rows)
    assert maxclosure<1e-8
    result=json.loads((OUT/'result.json').read_text());a=np.load(OUT/'p4-32.npz');time_errors={}
    for name in ['temperature','velocity','scalar']:
        error=float(np.max(abs(a[name]-d[name][::2]))/np.max(abs(d[name])))
        assert abs(error-result['time_relative'][name])<1e-14;time_errors[name]=error
    prior=json.loads((task.prior.OUT/'p4-64.json').read_text())['history'][-1]
    current=rows[-1]['history'][-1]
    data=dict(classification='Counterexample candidate',
        independent_lapse_passed=bool(max(bulk_error,surface_error,total_error)<.02),
        bulk_quadrature_relative=bulk_error,surface_quadrature_angular_relative=surface_error,whole_lapse_relative=total_error,
        surface_lapse_independent=float(surface8),surface_lapse_production=surface_reference,
        moving_mass_jump_relative=mass_error,moving_mass_jump_passed=mass_error<1e-7,
        local_support_mass_geom_over_R=float(support_mass),maximum_closure_relative=maxclosure,
        maximum_closure_iterations=max(row['max_closure_iterations'] for row in rows),time_relative=time_errors,
        old_surface_luminosity_deviation=prior['outgoing_luminosity_relative'],new_surface_luminosity_deviation=current['outgoing_luminosity_relative'],
        luminosity_deviation_ratio=current['outgoing_luminosity_relative']/prior['outgoing_luminosity_relative'],
        seconds=time.monotonic()-start,accounted_numerical_seconds=result['seconds']+time.monotonic()-start,
        external_support_gravity_closed=False,radiation_corrected_scalar_junction_closed=False,spatial_failure_resolved=False,final_charge_solved=False,full_goal_complete=False)
    task.write(OUT/'audit.json',data);signal.alarm(0);print(json.dumps(data),flush=True)


if __name__=='__main__':main()
