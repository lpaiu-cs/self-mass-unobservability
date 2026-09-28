"""Saved-source decomposition only; never reclassifies the failed join audit."""
from pathlib import Path
from types import MethodType
import inspect
import json
import signal
import time
import textwrap
import numpy as np
import def_native_material_join as m
import verify_native_material_join as v

out=m.OUT
m.write(out/'diagnostic-plan.json',dict(classification='Counterexample candidate',seconds=20,new_fluid_steps=0,
    question='Does the upwind acoustic pressure-jump mass flux dominate the new baryon charge? Does the failed endpoint neutral sum fit independently bounded floating-point cancellation?',
    method='Evaluate only the saved17 states. Integrate the explicit acoustic dissipative mass flux, retain the same retarded rest-charge operator, and compare linear/cubic histories. Re-form endpoint divergences in long double without changing the original executed RHS.',
    stop='Keep the original failed audit and time-converged negative charge. A reconstructed contribution is a diagnostic, not a new accepted physical trajectory.'))
signal.signal(signal.SIGALRM,m.old.optical.timeout);signal.alarm(20);start=time.monotonic()
model=m.Coupled();z=np.load(out/'coupled-128.npz');data=dict(np.load(out/'source-128.npz'));t=data['t']
ids=[np.argmin(abs(z['snapshot_t']-tt)) for tt in t];flux=[];adv=[]
for tt,i in zip(t,ids):
    model.h=z['snapshot_h'][i];model.Pi=z['snapshot_Pi'][i];model.mass=model.mass0-model.h[1:]+model.h[:-1];model.set_material(tt)
    p=model.bulk.eos.gas(z['snapshot_theta'][i],z['snapshot_eta'][i])[0];dp=p-model.f0['p0'];Z=np.sqrt(model.cx*model.face_rho*model.face_K)
    f=model.area_gas*model.face_rho*(np.r_[dp[0],dp]-np.r_[dp,dp[-1]])/(2*Z);f[0]=f[-1]=0;flux.append(f)
    vel=model.velocity()*m.C;vv=(np.r_[0.,vel]+np.r_[vel,0.])/2;vv[0]=vv[-1]=0
    adv.append(model.area_gas*model.face_rho*vv)
rows=[]
for kind in ['linear','cubic']:
    h=v.previous.prior.green.polynomial(t,np.array(flux),kind).antiderivative()(t)
    d=dict(data);d['nonrest_trace_erg']=np.zeros_like(data['nonrest_trace_erg']);d['baryon_g']=np.zeros_like(data['baryon_g']);d['baryon_g'][:,:model.bulk.n]=-np.diff(h,axis=1)
    _,q,_=v.previous.direct(d,model)
    rows.append(dict(kind=kind,endpoint_diffusive_rest_charge=float(q[-1]-q[0]),endpoint_diffusive_mass_norm_g=float(abs(-np.diff(h[-1])).sum())))
# Recover the exact endpoint flux arrays from the already executed expression.
src=textwrap.dedent(inspect.getsource(m.Coupled.material_rhs));src=m.replace(src,'return [mass,force,-np.diff(energy),-np.diff(neutral)]','return mass,force,energy,neutral')
ns=dict(vars(m));exec(compile(src,__file__,'exec'),ns)
model.recover_material(model.material_state(z['u'],z['eta']),m.old.END,z['theta'])
f=model.flow;V=f.primitive(z['U']);_,R=f.reconstruct(V,m.old.END);model.join_flux(R[:,0],V[3,0])
mass,force,energy,neutral=MethodType(ns['material_rhs'],model)(z['u'],z['theta'],z['eta'])
audit=[]
for key,F in [('mass',mass),('energy',energy),('neutral',neutral)]:
    D=-np.diff(F);res=float(D.sum()+F[-1]-F[0]);bound=float(16*np.finfo(float).eps*(abs(F).sum()+abs(D).sum()))
    LD=np.longdouble;fine=-np.diff(F.astype(LD));relative=float(abs(fine.sum(dtype=LD)+LD(F[-1])-LD(F[0]))/max(abs(F[-1]),1.))
    audit.append(dict(quantity=key,original_residual=res,roundoff_bound=bound,within_roundoff=bool(abs(res)<=bound),longdouble_port_relative=relative))
actual=json.loads((out/'result.json').read_text())['endpoint_direct_relative']
result=dict(classification='Counterexample candidate',diagnostic_completed=True,endpoint_direct=actual,pressure_jump_components=rows,
    diffusive_rest_over_direct=float(rows[0]['endpoint_diffusive_rest_charge']/actual),arithmetic=audit,seconds=time.monotonic()-start,
    original_join_audit_passed=False,physical_sign_established=False,new_fluid_steps=0)
m.write(out/'diagnostic.json',result);print(json.dumps(result),flush=True)
