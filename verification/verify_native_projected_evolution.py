"""Actual projected-state sources; no inherited old-background charge bound."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect
import json
import signal
import sys
import time
import numpy as np
import def_native_projected_evolution as run
import verify_native_material_join as joined
import verify_native_interior_feedback as shared

OUT=run.OUT;write=run.write;sha=run.sha;C=run.C


def geometry(m):
    return run.Geometry(m) if isinstance(m.bg,run.Background) else run.green.Geometry(m)


green=SimpleNamespace(**dict(vars(run.green),Geometry=geometry))
charge=run.previous.feedback.prior
charge_ns=dict(vars(charge),green=green)
# Integrated source contents are distributed with the actual proper-volume
# measure inside their cells. Normalize locally so packet/baryon contents are
# exact for both the new geometry and the saved old midpoint-volume model.
integration=inspect.getsource(charge.integrate)
integration=run.previous.replace(integration,'w,delay,*_=geom(r-model.m.RJ)',
    'w,delay,_,B,_=geom(r-model.m.RJ)\n    measure=(wq*r*r*B).reshape(len(edges)-1,order)\n    wq=(measure/measure.sum(1)[:,None]).ravel()')
exec(compile(integration,__file__,'exec'),charge_ns)
charge_ns['rays']=FunctionType(charge.rays.__code__,charge_ns,argdefs=charge.rays.__defaults__)
prior=SimpleNamespace(**charge_ns)
source=run.previous.replace(joined.source,'bp0=b.eos.base.gas(np.zeros(b.n),np.zeros(b.n))[0]',"bp0=model.f0['p0']")
ns=dict(vars(shared),run=run,OUT=OUT,prior=prior)
exec(compile(source,__file__,'exec'),ns)
for name in ['direct','readout']:
    fn=getattr(shared,name);ns[name]=FunctionType(fn.__code__,ns,argdefs=fn.__defaults__)


def prepare():
    assert not (OUT/'readout-plan.json').exists()
    plan=json.loads((run.previous.OUT/'readout-plan.json').read_text())
    plan.update(claim='Read actual projected-background material/photon trajectories and recompute direct charge and the separately conditional outgoing-photon normalization.',
        geometry='Use corrected mass,lapse,phi,Phi throughout the represented source and the exact vacuum beginning at the outer photon face. Within a cell, distribute its integrated source by normalized Jordan proper volume, not uniform coordinate radius. Apply the same proper-volume readout to the saved Phase119 sources as an independent old-background comparator.',
        limits='Frozen corrected initial metric, first-order density/inventory EOS, finite microphysics and resolution. Initial moment matching is not evolved metric/scalar closure. No old GR lower bound is transferred; direct-plus-photon is not the final physical charge.',
        budget=dict(saved_readout_seconds=60,native_endpoint_seconds=25,native_calls=160,new_fluid_steps=0),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),run.INPUT/'balanced-20.npz',run.previous.OUT/'source-128.npz',run.previous.OUT/'wave-128.npz']} )
    write(OUT/'readout-plan.json',plan)
    (OUT/'expanded-source.py').write_text(source);(OUT/'expanded-integrate.py').write_text(integration)


def readout():
    ns['readout']()
    start=time.monotonic();old_model=run.previous.Coupled();d=dict(np.load(run.previous.OUT/'source-128.npz'))
    t,q,parts=ns['direct'](d,old_model)
    old=np.load(run.previous.OUT/'wave-128.npz');new=np.load(OUT/'wave-128.npz')
    # Compare at the new observer times without changing either physical run.
    qq=np.interp(new['t'],t,q-q[0]);olduniform=np.interp(new['t'],old['t'],old['direct_relative'])
    np.savez_compressed(OUT/'old-background-proper-volume.npz',t=t,direct=q,components=parts)
    scale=max(abs(new['direct_relative']))
    write(OUT/'comparison.json',dict(classification='Counterexample candidate',
        projected_endpoint=float(new['direct_relative'][-1]),old_same_volume_rule_endpoint=float(qq[-1]),
        old_uniform_coordinate_endpoint=float(olduniform[-1]),
        projected_vs_old_same_rule=float(max(abs(new['direct_relative']-qq))/scale),
        old_volume_rule_change=float(max(abs(qq-olduniform))/max(abs(qq))),
        new_initial_state_actually_evolved=True,new_GR_lower_bound_certified=False,
        final_charge_solved=False,seconds=time.monotonic()-start))


def audit():
    assert not (OUT/'audit.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,run.old.optical.timeout);signal.alarm(25)
    m=run.Coupled();b=m.bulk;z=np.load(OUT/'coupled-128.npz');m.Pi=z['Pi'];m.h=z['h'];m.mass=m.mass0-m.h[1:]+m.h[:-1];m.set_material(run.old.END)
    # Direct native endpoint checks exclude the explicitly first-order
    # advected-inventory term, just as the inherited endpoint audit does.
    native=run.previous.chem.old.Native(cap=160);theta=z['theta'];eta=z['eta'];xi=b.eos.xi.copy();b.eos.xi[:]=0
    p,u,*_=b.eos.gas(theta,eta);errors=[];anchor=b.eos.base.d
    for j in range(b.n):
        run.previous.chem.setup(native,anchor,j)
        density=np.log1p(b.eos.density_shift[j])+np.log1p(b.eos.x[j])
        state=native.state(density,np.log(anchor['T'][j])+b.eos.theta0[j]+theta[j],anchor['y0'][j]*(1+eta[j]))
        errors.append([abs(p[j]/state['raw'][1]-1),abs(u[j]/state['raw'][2]-1)])
    b.eos.xi=xi
    f=m.flow;V=f.primitive(z['U']);_,right=f.reconstruct(V,run.old.END);flux=m.join_flux(right[:,0],V[3,0])
    converted=flux*4*np.pi*m.m.RJ**2*f.eos.rho0*np.array([C,C*C,C**3,C*f.eos.nH])
    face=float(np.max(abs(converted-m.mflux)/np.maximum(abs(m.mflux),1.)))
    row=dict(classification='Counterexample candidate',passed=bool(np.max(errors)<.002 and face<1e-12),
        native_endpoint_pressure_energy_relative=float(np.max(errors)),native_calls=native.ion.calls,
        shared_material_face_units=face,shared_photon_face_area_identical=bool(b.area[-1]==m.area[0]),
        corrected_initial_geometry_used_in_material_and_photon_operators=True,
        full_neutral_trajectory_ledger=False,uniform_EOS_derivative_bound=False,full_GR_feedback=False,final_charge_solved=False,seconds=time.monotonic()-start)
    write(OUT/'audit.json',row);print(json.dumps(row),flush=True);signal.alarm(0);assert row['passed']


if __name__=='__main__':globals()[sys.argv[1]]()
