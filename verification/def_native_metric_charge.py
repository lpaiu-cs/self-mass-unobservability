"""Apply repaired gas histories to explicit matter-induced GR scalar sources.

Counterexample candidate: prescribed gas, local acoustic bulk response, direct
scalar Green term, metric stress source and the distributed layer mass source.
This is not the iterated metric/scalar/fluid solution or a total charge bound.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect
import json
import signal
import time
import numpy as np
import def_native_metric_release as evolution
import def_native_release_charge as prior

OUT=evolution.OUT/'charge'
C,G=prior.C,prior.G


def adapter(n):
    flow=evolution.Flow(n);m=flow.base;m.eos=flow.eos;m.primitive=flow.primitive
    return m


class MetricGeometry(prior.Geometry):
    def metric(self,x):
        w,a,B,re,phi=super().metric(x)
        Phi=self.bg.sample(re)['v']/self.m.R
        return a*Phi,a,B,re,phi


def readouts():
    env=dict(vars(prior),OUT=OUT,prior=SimpleNamespace(Flow=adapter,OUT=evolution.OUT))
    direct=FunctionType(prior.readout.__code__,env,argdefs=prior.readout.__defaults__)
    metric=dict(env,OUT=OUT/'metric-stress',Geometry=MetricGeometry)
    source=inspect.getsource(prior.acoustic)
    assert source.count('u+p/rr-3*gamma*p/rr')==1
    source=source.replace('u+p/rr-3*gamma*p/rr','u+p/rr-gamma*p/rr')
    exec(compile(source,__file__,'exec'),metric)
    source=inspect.getsource(prior.readout)
    assert source.count('rho*u-3*p')==1
    source=source.replace('rho*u-3*p','rho*u-p')
    exec(compile(source,__file__,'exec'),metric)
    return direct,metric['readout']


def mass_source(n,utimes):
    """Forced J constraint, using actual inner Killing energy and baryon ports.

    The algebraic response r^2*b*Phi*f and f-dependent part of J are not
    included in this first source application. They require the wave iterate.
    """
    m=adapter(n);geom=prior.Geometry(m);saved=np.load(evolution.OUT/f'cells-{n}.npz')
    h=saved['history'];t=np.r_[0,h[:,0]]
    states=np.concatenate([saved['initial'][None],saved['snapshots']]).astype(np.longdouble)
    scale=np.longdouble(4*np.pi*m.RJ**2*m.eos.rho0)
    delta=(states-states[0])*m.vol*scale
    portB=np.r_[0,h[:,3]].astype(np.longdouble)*scale
    portK=np.r_[0,h[:,4]-h[:,5]].astype(np.longdouble)*scale
    w,d,a,B,re=geom(m.x);p=m.bg.sample(re);r=re*m.R
    b=1-2*p['m']/re;A=a/p['N'];Phi=p['v']/m.R
    Eg=p['e']/m.R**2;Pg=p['p']/m.R**2;alpha=-4*geom.metric(m.x)[4]
    dr=m.dx*B*np.sqrt(b)/A
    # kappa=mu*sqrt(b)/N. Its tiny variation is retained with expm1;
    # forming 1+variation before multiplying a large rest mass loses it.
    logk_rate=-4*np.pi*r*A**4*(Eg+Pg)/b
    logk=np.cumsum(logk_rate*dr)-logk_rate*dr/2
    eps=np.expm1(logk).astype(np.longdouble)
    energy=delta[:,2]+np.longdouble(m.a0*m.eos.cx)*delta[:,0]
    inner=-portK-np.longdouble(m.a0*m.eos.cx)*portB
    prefix=np.cumsum(energy,axis=1)-energy/2
    weighted=np.cumsum(energy*eps,axis=1)-energy*eps/2
    J=(G/C**2)*(np.sqrt(b)/p['N'])*np.exp(-logk)*(inner[:,None]+prefix+weighted)
    coefficient=2*(p['N']**2*b)*Phi/(r*b*b)*(1+4*np.pi*r*r*A**4*(Pg-Eg))
    coefficient-=8*np.pi*p['N']**2/b*alpha*A**4*(Eg-3*Pg)
    Hj=prior.polynomial(t,J).antiderivative();q=[]
    for tt in utimes:
        q.append(C/2*np.sum((m.dx*B/a)*coefficient*prior.paired(Hj,tt+d),dtype=np.longdouble))
    q=np.asarray(q,float)
    # Audit the integrated constraint independently with right-face values.
    increments=energy*(1+eps)
    cumulative=inner[:,None]+np.cumsum(increments,axis=1)
    residual=np.diff(np.column_stack([inner,cumulative]),axis=1)-increments
    relative=float(np.max(abs(residual))/max(np.max(abs(increments)),1e-100))
    # Full rest-mass conservation uses the inferred positive unresolved tail.
    # Evaluate the outer monopole from its centered identity, never subtract
    # the entire large rest inventory to find the tiny scattering work.
    missing=portB-np.sum(delta[:,0],axis=1)
    tailK=np.r_[0,h[:,8]].astype(np.longdouble)*scale
    netK=-portK+np.sum(delta[:,2],axis=1)+tailK
    work=np.r_[0,h[:,5]].astype(np.longdouble)*scale
    energy_relative=float(np.max(abs(netK-work))/max(np.max(abs(work)),1e-100))
    assert relative<1e-12 and energy_relative<1e-5
    np.savez_compressed(OUT/f'mass-source-{n}.npz',source_seconds=t,radius_cm=r,J_forced_cm=np.asarray(J,float),
        scalar_source_coefficient_cm3=coefficient,u_seconds=utimes,charge_cm=q,normalized_charge=-q/(m.bg.M*m.R),
        unresolved_baryon_g=np.asarray(missing,float),net_Killing_energy_erg=np.asarray(netK*C*C,float),
        scattering_work_erg=np.asarray(work*C*C,float),kappa_log_variation=logk)
    row=dict(classification='Counterexample candidate',cells=n,
        endpoint_layer_mass_source_normalized=float(-q[-1]/(m.bg.M*m.R)),
        integrated_constraint_relative=relative,global_energy_identity_over_scattering_work=energy_relative,
        maximum_kappa_log_variation=float(max(abs(logk))),
        maximum_forced_J_cm=float(np.max(abs(J))),
        scope='The actual inner baryon and Killing-energy debit fixes J at the layer cut, then the distributed saved lab-frame gas energy is integrated through the layer. This is the matter forcing component. Interior profile, unresolved tail position, causal photon stress, f-dependent mass feedback and wave potential remain outside this component.')
    prior.write(OUT/f'mass-source-{n}.json',row)
    return -q/(m.bg.M*m.R),row


def prepare():
    assert not OUT.exists();OUT.mkdir();(OUT/'metric-stress').mkdir()
    paths=[Path(__file__),Path(evolution.__file__),Path(prior.__file__),evolution.old.OUT/'runtime-columns.npz']
    paths += [evolution.OUT/f'cells-{n}.npz' for n in [896,1792]]
    prior.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply the repaired conservative trajectories to the direct trace, boost-invariant metric stress E-Pr, and the actual lab-energy forcing of the GR mass constraint. Measure the separate outgoing source components on both saved grids.',
        equations='S_stress=-4*pi*G/c^4*N^2*A^4*r^2*Phi*d(E-Pr); J_prime+r*Phi^2*J=4*pi*G/c^4*r^2*A^4*dE for f=0. The full symbolic constraint is saved by the parent phase.',
        seconds=45,CPU_threads=1,memory_GB=2,native_calls=0,new_fluid_steps=0,grid_gate=.02,
        measured_basis='The prior two direct readouts cost8.72seconds. Four equivalent readouts plus two source integrations and setup are budgeted at45seconds; stop instead of expanding if this is insufficient.',
        stop='No EOS extension, new gas path, new horizon or automatic wave iteration. Do not classify the sum of these components as the full charge.',
        bindings={str(p.relative_to(prior.old.ROOT)):prior.old.photons.digest(p) for p in paths}))


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(45)
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert prior.old.photons.digest(prior.old.ROOT/p)==h,p
    times=np.load(prior.OUT/'cells-1792-linear-g12.npz')['u_seconds'];raw=np.load(prior.OUT/'native-bulk.npz')['raw']
    direct,stress=readouts();values={};rows={}
    for n in [896,1792]:
        da,dr=direct(n,'linear',times,raw);sa,sr=stress(n,'linear',times,raw);ja,jr=mass_source(n,times)
        values[n]=np.array([da,sa,ja]);rows[n]=dict(direct=dr,metric_stress=sr,layer_mass=jr)
    errors=np.max(abs(values[896]-values[1792]),axis=1)/np.max(abs(values[1792]),axis=1)
    previous=np.load(evolution.old.OUT/'charge/cells-1792-linear-g12.npz')['normalized_charge']
    change=float(max(abs(values[1792][0]-previous))/max(abs(values[1792][0])))
    total=np.sum(values[1792],axis=0)
    np.savez_compressed(OUT/'components.npz',u_seconds=times,direct_trace=values[1792][0],metric_stress=values[1792][1],layer_mass=values[1792][2],component_sum=total)
    result=dict(classification='Counterexample candidate',passed=bool(max(errors)<.02),
        component_grid_relative=errors.tolist(),endpoint_components_normalized=values[1792][:,-1].tolist(),
        endpoint_component_sum=float(total[-1]),boundary_repair_wave_change_relative=change,
        endpoint_inner_energy_mismatch_erg=rows[1792]['direct']['endpoint_bulk_energy_mismatch_erg'],
        seconds=time.monotonic()-start,source_mass_constraint_applied=True,distributed_metric_stress_applied=True,
        scalar_potential_and_feedback_iterated=False,causal_scattered_photon_metric_source_applied=False,
        full_inner_material_response=False,final_charge_solved=False,full_goal_complete=False)
    prior.write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
