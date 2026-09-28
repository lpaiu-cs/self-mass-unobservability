"""Recover material and metric states from the saved coupled time endpoint."""
from pathlib import Path
import json
import numpy as np
import def_reactive_cauchy as coupled

h=coupled.h
OUT=coupled.OUT/'readout'


def main():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(coupled.__file__),coupled.OUT/'plan.json',coupled.OUT/'result.json',
           coupled.OUT/'fine-64.npz',coupled.OUT/'fine-rays.npz']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        claim='Reconstruct baryon density, pressure, temperature, composition, radial motion and metric constraints from the actual saved linear time solution. Independently compare the conservative mass defect to direct ray crossings.',
        scope='Postprocessing only; no added time paths, rays or native calls. Lapse normalized at the causally unperturbed 2Rstar boundary; this is not an asymptotic outgoing charge extraction.',
        gates=dict(surface_and_outer_mass_source_relative=2e-13,linear_pressure_operator_replay=1e-11),
        budget=dict(hard_seconds=30,new_time_steps=0,native_calls=0)))
    plan=json.loads((coupled.OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert h.digest(h.ROOT/rel)==sha,rel
    ray=dict(np.load(coupled.OUT/'fine-rays.npz'));radiation=coupled.Radiation(ray);bg=coupled.Background(radiation,2)
    saved=np.load(coupled.OUT/'fine-64.npz');assert np.array_equal(saved['grid'],bg.grid)
    y=saved['response'];ydot=saved['velocity'];r=bg.grid;N,a=radiation.geometry.metric(r)
    z,eta,f,V=y.T;xi=r*z;b=1/a**2;p,e,phi,v,ga=[bg.nodes[k] for k in ['p','e','phi','v','gamma']]
    rr,loss,E,P,J=radiation.source(1.,bg.nodes,bg.node_projection);A4=np.exp(-8*phi*phi);w=e+p
    Dm=r*r*b*v*f-(4*np.pi*r*r*A4*p+r*r*b*v*v/2)*xi+J
    dp=p*(eta-ga*rr);de=w*eta/ga-loss;drho=eta/ga;scalar=f-xi*v
    dl=np.divide(Dm,r,out=np.zeros_like(r),where=r>0)
    dl-=np.divide(bg.nodes['m']*xi,r*r,out=np.zeros_like(r),where=r>0);dl/=b
    fn,_=coupled.reactive.symbolic()
    mid=bg.mid;vals=fn(*[mid[k] for k in ['r','m','p','e','phi','v']],-4.)
    _,g,F=vals[:3];dg=vals[9:15]
    # Delta nu_prime = Delta g - beta*Phi*f - alpha*V.
    # Eulerian increment subtracts the advected background gradient. Evaluate
    # that gradient from the same background equations, not differencing data.
    M=vals[0];alpha=-4*mid['phi'];wm=mid['e']+mid['p']
    gprime=dg[0]+dg[1]*M+dg[2]*(-wm*g)+dg[4]*mid['v']+dg[5]*F
    nusecond=gprime+4*mid['v']**2-alpha*F
    ym=(y[:-1]+y[1:])/2;rm=mid['r'];xim=rm*ym[:,0]
    rrm,lm,Em,Pm,Jm=radiation.source(1.,mid,bg.mid_projection)
    bm=1-2*mid['m']/rm;A4m=np.exp(-8*mid['phi']**2)
    Dmm=rm*rm*bm*mid['v']*ym[:,2]-(4*np.pi*rm*rm*A4m*mid['p']+rm*rm*bm*mid['v']**2/2)*xim+Jm
    increments=[xim,Dmm,mid['p']*(ym[:,1]-mid['gamma']*rrm),wm*ym[:,1]/mid['gamma']-lm,ym[:,2],ym[:,3]]
    Dg=sum(q*u for q,u in zip(dg,increments))+4*np.pi*rm*Pm/bm
    euler_nuprime=Dg+4*mid['v']*ym[:,2]-alpha*ym[:,3]-xim*nusecond
    integral=np.r_[0,np.cumsum(np.diff(r)*euler_nuprime)];euler_nu=integral-integral[-1]
    # Recover material values at the original native centres, not the vacuum
    # coordinate labels outside the star or the zero-temperature endpoint.
    d=bg.native;s=bg.reactions;rn=d['radius_cm']/(100*radiation.geometry.R)
    interp=lambda q:np.interp(rn,r,q)
    density=interp(drho);temperature=d['thermo'][:,4]*density+interp(bg.nodes['theta_rate'])*radiation.geometry.tc
    composition=d['A'][:,None]*d['N'][:,None]*s['dxdt']*radiation.geometry.tc
    physical_v=a/N*r*ydot[:,0]*h.gr.C
    direct=[];inventory=radiation.inventory(1.);net=inventory[:,0]-radiation.power*radiation.geometry.tc
    face=-np.r_[0,np.cumsum(net)];total=ray['emitter_power'].sum()*radiation.geometry.tc
    for key,position in [('surface_crossings',1.),('outer_crossings',2.)]:
        q=ray[key];em=ray['emitter_power'][q[:,0].astype(int)]
        escaped=float((em*q[:,1])@np.maximum(1-q[:,2],0))*radiation.geometry.tc
        face_id=int(np.flatnonzero(ray['edges']==position)[0]);error=abs(face[face_id]-escaped)/total
        assert error<2e-13,error
        node=int(np.flatnonzero(r==position)[0]);expected=-h.gr.G*1e-7*escaped/(h.gr.C**4*radiation.geometry.R*N[node]*a[node])
        direct.append(dict(radius_over_R=position,escaped_erg=escaped,mass_defect_geom_m=float(J[node]*radiation.geometry.R),
                           independent_expected_geom_m=float(expected*radiation.geometry.R),energy_norm_relative=float(error)))
    # Independent zero-source operator replay: the new inertia sign and lifted
    # variables reduce to the earlier adiabatic harmonic operator at rr=0.
    inside={k:q[mid['r']<1] for k,q in mid.items()};omega=.37
    original,_,_,_=coupled.surface.old.operators(inside,omega)
    new,B,_,_=coupled.operators(inside,fn)
    expected=new[:,:,:4]-omega*omega*B
    operator_error=float(np.max(abs(original-expected))/max(abs(original).max(),1e-100))
    assert operator_error<1e-11,operator_error
    np.savez_compressed(OUT/'endpoint.npz',radius_over_R=r,Delta_m_over_R=Dm,Delta_lambda=dl,
        Euler_delta_nu=euler_nu,Delta_e=de,Delta_p=dp,Delta_log_rho=drho,Euler_delta_phi=scalar,
        velocity_m_s=physical_v,radiation_E=E,radiation_Pr=P,mass_defect_over_R=J,
        native_radius_over_R=rn,native_Delta_log_T=temperature,native_Delta_X=composition)
    massweights=d['dm']/d['dm'].sum()
    result=dict(classification='Counterexample candidate',passed=True,
        adiabatic_operator_replay_relative=operator_error,direct_ray_mass_checks=direct,
        endpoint_seconds=radiation.geometry.tc,native_max_abs_log_rho=float(abs(density).max()),
        native_max_abs_log_T=float(abs(temperature).max()),native_mass_RMS_log_T=float(np.sqrt(massweights@(temperature**2))),
        native_max_abs_Delta_X=float(abs(composition).max()),
        native_velocity_mass_RMS_m_s=float(np.sqrt(massweights@(interp(physical_v)**2))),
        native_scalar_mass_RMS=float(np.sqrt(massweights@(interp(scalar)**2))),
        max_abs_Euler_delta_nu=float(abs(euler_nu).max()),max_abs_Delta_lambda=float(abs(dl).max()),
        perturbation_smallness_is_nonlinear_error_bound=False,mechanical_spatial_convergence=False,
        lapse_is_postprocessed_constraint=True,full_dynamic_charge_solved=False,
        runtime=dict(production_seconds=16.809798546,process_wall_seconds=20.09,peak_RSS_KiB=402348,
                     CPU_workers=1,BLAS_threads=1,native_calls=0,production_time_steps=240))
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
