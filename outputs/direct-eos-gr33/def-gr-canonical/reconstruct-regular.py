"""Independent physical pressure/metric reconstruction of the saved trajectory."""
from pathlib import Path
import signal
import time
import numpy as np
import def_gr_energy_evolution as evolution

import def_gr_canonical_regular as regular
fem=evolution.fem;fem.canonical=regular;OUT=fem.OUT/'regular'
assert not (OUT/'reconstruction.json').exists()
fem.write(OUT/'reconstruction-plan.json',dict(classification='Counterexample candidate',
    claim='Reconstruct physical pressure, density, temperature and the same GR mass constraint from the actually propagated endpoint. Check baryon pressure identity independently and inspect source quadrature alignment.',
    budget=dict(hard_seconds=60,new_time_paths=0,new_EOS_calls=0),
    source=evolution.digest(Path(__file__)),endpoint=evolution.digest(OUT/'fine-512.npz')))
signal.alarm(60);start=time.monotonic();model=fem.Model(evolution.BANK/'fine-bank.npz');d=dict(np.load(OUT/'fine-512.npz'))
q=d['q'];assert np.array_equal(d['grid'],model.grid)
nodal=np.r_[q,0.][model.indices];dx=np.diff(model.grid);local=np.array([(1-1/np.sqrt(3))/2,(1+1/np.sqrt(3))/2])
value=((1-local)[None,:,None]*nodal[:-1,None,:]+local[None,:,None]*nodal[1:,None,:]).reshape(-1,2)
derivative=np.repeat((nodal[1:]-nodal[:-1])/dx[:,None],2,axis=0)
src=np.column_stack([a@d['heat_energy'] for a in fem.source_points(model.heat,model.points)])
reference=model.heat.source(1.,model.points)[[0,1,4]].T
source_errors=[float(np.max(abs(src[:,j]-reference[:,j]))/max(np.max(abs(reference[:,j])),1e-100)) for j in range(3)]
assert max(source_errors)<1e-12
data=model.canonical_data;A=data[:,:4].reshape(-1,2,2);B=data[:,4:8].reshape(-1,2,2)
g=np.einsum('nij,nj->ni',data[:,16:22].reshape(-1,2,3),src)
R=data[:,28:32].reshape(-1,2,2);S=data[:,32:36].reshape(-1,2,2)
P=np.einsum('nij,nj->ni',B,derivative-np.einsum('nij,nj->ni',A,value)-g)
physical=P-np.einsum('nij,nj->ni',S,value)/2
point=model.points;inside=point['r']<1;r=point['r'][inside];z=value[inside,0];psi=value[inside,1]
eta=physical[inside,0]/R[inside,0,0];rr,loss,J=src[inside].T
m,p,e,phi,Phi,Gamma=[point[k][inside] for k in ['m','p','e','phi','v','gamma']]
b=1-2*m/r;A4=np.exp(-8*phi**2);alpha=-4*phi;xi=r*z;f=psi+xi*Phi
Dm=r*r*b*Phi*f-(4*np.pi*r*r*A4*p+r*r*b*Phi*Phi/2)*xi+J
dl=(Dm/r-m*xi/r**2)/b
eta_baryon=-Gamma*(r*derivative[inside,0]+3*z+dl+3*alpha*f+rr)
condition=abs(eta)+abs(Gamma)*(abs(r*derivative[inside,0])+3*abs(z)+abs(dl)+3*abs(alpha*f)+abs(rr))
error=float(np.max(abs(eta-eta_baryon)/np.maximum(condition,1e-100)))
assert error<1e-12,error
eta_ad=eta+Gamma*rr;density=eta_ad/Gamma
h=fem.base.task.h;geo=h.gr.G*.1*model.original.radiation.geometry.R**2/h.gr.C**4
rho=np.interp(r,model.heat.rnative,model.heat.d['raw'][::-1,0])
cvT=np.interp(r,model.heat.rnative,model.heat.d['thermo'][::-1,5])
temperature=point['adiabatic_T_rho'][inside]*density-loss/(rho*geo*cvT)
assert np.all(np.isfinite(temperature))
edges=model.heat.edges;ids=np.clip(np.searchsorted(edges,model.grid),0,len(edges)-1)
nearest=np.minimum(abs(model.grid-edges[ids]),abs(model.grid-edges[np.maximum(ids-1,0)]))
crossings=np.searchsorted(model.grid,edges,side='right')-1
unaligned=[int(i) for i in range(1,len(edges)-1) if not np.any(abs(model.grid-edges[i])<2e-14)]
np.savez_compressed(OUT/'physical-endpoint.npz',radius=r,zeta=z,Eulerian_scalar=psi,Lagrangian_scalar=f,
    physical_pressure=eta,adiabatic_pressure=eta_ad,delta_log_density=density,delta_log_temperature=temperature,
    delta_mass=Dm,delta_lambda=dl,source_rho_ref=rr,source_loss=loss,source_J=J)
result=dict(classification='Counterexample candidate',physical_pressure_metric_reconstruction_checked=True,
    source_relative_errors=source_errors,independent_baryon_pressure_backward_error=error,
    maximum_abs=dict(physical_pressure=float(abs(eta).max()),delta_log_density=float(abs(density).max()),
        delta_log_temperature=float(abs(temperature).max()),delta_mass=float(abs(Dm).max()),delta_lambda=float(abs(dl).max())),
    heat_faces_not_aligned_to_spatial_nodes=len(unaligned),heat_internal_faces=len(edges)-2,
    unaligned_active_faces=int(np.isin(model.heat.face_ids,unaligned).sum()),
    seconds=time.monotonic()-start,original_failure_resolved=False,full_dynamic_charge_solved=False,
    limitation='Algebraic pressure/baryon/metric consistency does not bound spatial interpolation or prove the unsolved full nonlinear/outer-current/photon problem.')
fem.write(OUT/'reconstruction.json',result);signal.alarm(0);print('RECONSTRUCT',result,flush=True)
