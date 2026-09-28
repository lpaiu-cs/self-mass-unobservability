"""Return the updated free material motion to actual photon/thermal/H stages."""
from pathlib import Path
from types import FunctionType
import json,resource,sys,time
import numpy as np
import def_native_updated_material_return as matter
import def_native_matter_photon_feedback as base
import def_native_stage_feedback as stage
import def_native_parallel_return as parallel

prior=matter.prior;OUT=matter.OUT/'joint-feedback';write=matter.write;sha=matter.sha;AMP=matter.AMP
PATHS=matter.PATHS;replace=prior.replace


class Response(base.Response):
    def __init__(self,reference):
        base.mono.Response.__init__(self,reference);self.material=matter.Material(reference)
        d=np.load(matter.OUT/f'steps-128-reference-{reference}.npz');ids=[int(np.argmin(abs(d['t']-t))) for t in self.t]
        assert np.max(abs(d['t'][ids]-self.t))<1e-18
        self.motion=d['history_scaled'][ids].copy()
        self.energy_offset=(self.a-self.model.m.a0)*self.model.cx*base.C**2*self.motion[:,0]
        self.motion[:,2]-=self.energy_offset
        self.mechanical=self.motion[:,[2,3]]-self.material.transfer[:,[2,3]]
        self.xi=np.array([self.material.model.mech.xi@np.r_[0.,-np.cumsum(z[0,:self.nb])] for z in self.motion])
        self.velocity_jet_error=0.;self.mapping_error=0.;self.map_checks=[];self.angular=[];self.angular_times=[]

    def boundary_ports(self,t,x):
        k=max(0,min(np.searchsorted(self.t,t,side='left')-1,15));f=(t-self.t[k])/(self.t[k+1]-self.t[k])
        I=(1-f)*self.I[k]+f*self.I[k+1]
        speed=((1-f)*self.g['delta_log_speed'][k]+f*self.g['delta_log_speed'][k+1])/AMP
        actual=x[-1]*self.scale+I[-1]*speed[-1]
        self.angular_times.append(t)
        self.angular.append(AMP*2*np.pi*base.C*self.area[-1]*(actual[self.mu>0]@(self.num*self.E)))
        return super().boundary_ports(t,x)


runner=stage.runner.replace('rtol=1e-10','rtol=1e-12').replace('err<1e-10','err<1e-12').replace('self.max_residual<1e-10','self.max_residual<1e-12')
runner=replace(runner,"times=z['t'].tolist();","self.angular_times=list(z['accepted_angular_times']);self.angular=list(z['accepted_angular_luminosity']);self.max_residual=previous['linear_relative'];self.max_iterations=previous['max_Krylov_iterations'];times=z['t'].tolist();")
runner=replace(runner,'material_history=gas_history)',
    'material_history=gas_history,accepted_angular_times=self.angular_times,accepted_angular_luminosity=self.angular)')
runner=replace(runner,'gas_history=gas_history)',
    'gas_history=gas_history,accepted_angular_times=self.angular_times,accepted_angular_luminosity=self.angular)')
namespace=dict(vars(base.mono),OUT=OUT);exec(compile(runner,__file__,'exec'),namespace);Response.run=namespace['run']
dispatch=FunctionType(parallel.dispatch.__code__,dict(vars(parallel),OUT=OUT,__file__=__file__))


def worker(steps,reference,label,limit,restart):
    cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));matter.configure();start=time.monotonic();cpu=time.process_time()
    row=Response(reference).run(steps,label,limit,restart)
    row.update(worker_wall_seconds=time.monotonic()-start,worker_CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    write(OUT/f'{label}.json',row);assert row['passed'],row


def prepare():
    assert not OUT.exists();OUT.mkdir();assert json.loads((matter.OUT/'audit.json').read_text())['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',before_checkpoint='8dc7db609',
        claim='Resolve the newly measured material/photon energy-H mismatch by applying actual updated baryon/momentum/inventory/noncollisional transport to simultaneous moving photon/energy/H stages; then return new paired transfers to material and GR.',
        trigger='The first updated material sweep differs from the prior radiation-only material history by1.02237 in energy and0.439882 in H. Small compact GR correction alone does not close this coupling.',
        reuse='Existing exact nonrest energy coordinate and conservative primitive maps, signed velocity/inventory jets, actual SDIRK stage times and1e-12 full residual. Same three original paths,531 cells,8 angles,152 frequencies,17 background knots and3.434ms.',
        angular='Save both actual accepted SDIRK-stage outgoing angular luminosity increments and times, alongside full angular snapshots, thermal/H state, collision-only transfers and radial ports. Avoid another replay to reconstruct the next exterior lapse input.',
        limits='Prescribed previous material baryon/momentum and mechanical transport are returned once. A subsequent free-material sweep and GR/lapse reapplication are required; this is not yet a fixed point.',
        budgets=dict(pilot_wall_seconds=65,production_wall_seconds=650,CPU_processes=3,threads_each=1,total_virtual_memory_GiB=9,new_native_states=0),
        forecast='Concurrent prefixes4/8/4 steps. Per path17 measured operator points plus remaining measured stepping cost and20s overhead;2x maximum must fit650s. New6-variable primitive jets require3GiB virtual memory per process, explicitly increased from the prior2GiB cap. Later Krylov cost remains extrapolated.',
        gates=dict(conservation=1e-8,linear=1e-12,velocity_jet=1e-4,mapping=1e-10,time=.02,background=.02,angular_port=1e-12),
        stop='Stop all workers on any failed path or650s. Preserve prefixes; no automatic extra path, finer clock/grid, longer horizon, native support expansion or gate relaxation.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(matter.__file__),Path(base.__file__),Path(stage.__file__),matter.OUT/'audit.json',prior.OUT/'audit.json']}))
    (OUT/'expanded-run.py').write_text(runner)


def pilot():
    assert not (OUT/'pilot.json').exists()
    rows,seconds=dispatch([(64,128,'pilot-64-128',4,None),(128,128,'pilot-128-128',8,None),(128,64,'pilot-128-64',4,None)],65)
    estimates=[r['operator_point_seconds']/r['operator_points']*17+r['stepping_seconds']/r['new_steps']*(r['steps']-r['completed_steps'])+20 for r in rows]
    a=np.load(OUT/'pilot-64-128.npz')['moments'][-1,[0,1,2,3,5,6]];b=np.load(OUT/'pilot-128-128.npz')['moments'][-1,[0,1,2,3,5,6]]
    errors=(np.sum(abs(a-b),axis=1)/np.maximum(np.sum(abs(b),axis=1),1.)).tolist()
    p=dict(classification='Counterexample candidate',rows=rows,equal_horizon=errors,forecast_each_seconds=estimates,upper_seconds=2*max(estimates),
        eligible=all(r['passed'] for r in rows) and max(errors)<.02 and 2*max(estimates)<650,seconds=seconds)
    write(OUT/'pilot.json',p);print(json.dumps(p),flush=True)
    if p['eligible']:write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_cap_seconds=650,
        paths=[[n,r,f'pilot-{n}-{r}'] for n,r in PATHS],
        bindings={str(f):sha(f) for f in [Path(__file__),Path(matter.__file__),Path(base.__file__),OUT/'plan.json',OUT/'pilot.json',matter.OUT/'audit.json']}))


production=FunctionType(parallel.production.__code__,dict(vars(parallel),OUT=OUT,dispatch=dispatch))


def angular_repair():
    original=json.loads((OUT/'pilot.json').read_text());assert original['eligible'];assert not (OUT/'corrected-pilot.json').exists()
    remaining=65-original['seconds'];assert remaining>0
    write(OUT/'angular-repair-plan.json',dict(classification='Counterexample candidate',
        failure='Short-pilot angular readout does not reproduce the shared physical port because the new writer used the final frequency scale for every frequency. The physical equation and port are unaffected.',
        repair='Apply the full frequency scale. Replay only4/8/4 pilot steps with observational capture; require every non-angular result array bitwise identical and angular energy equal to the actual port below1e-12. Preserve rejected angular files.',
        remaining_original_pilot_wall_seconds=remaining,total_pilot_wall_seconds=65,production_wall_seconds=650,
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'first-producer.py',OUT/'pilot.json',OUT/'angular-prefix-failure.json']}))
    specs=[(64,128,'corrected-pilot-64-128',4,None),(128,128,'corrected-pilot-128-128',8,None),(128,64,'corrected-pilot-128-64',4,None)]
    rows,seconds=dispatch(specs,remaining);errors=[];gamma=1-1/np.sqrt(2)
    for n,r,label,_,_ in specs:
        old=np.load(OUT/f'pilot-{n}-{r}.npz');new=np.load(OUT/f'{label}.npz')
        assert all(np.array_equal(old[k],new[k]) for k in old.files if not k.startswith('accepted_angular_'))
        h=prior.run.flow.old.END/n;flux=new['accepted_angular_luminosity']@((np.arange(4)*2+1)/32);pairs=flux.reshape(-1,2)
        err=float(abs(h*np.sum(pairs*[1-gamma,gamma])-new['radial_ports'][-1,1,1])/max(h*np.sum(abs(pairs)),1.));assert err<1e-12;errors.append(err)
    estimates=[r['operator_point_seconds']/r['operator_points']*17+r['stepping_seconds']/r['new_steps']*(r['steps']-r['completed_steps'])+20 for r in rows]
    p=dict(classification='Counterexample candidate',rows=rows,all_nonangular_arrays_bitwise=True,angular_port_relative=errors,forecast_each_seconds=estimates,
        upper_seconds=2*max(estimates),eligible=all(r['passed'] for r in rows) and 2*max(estimates)<650,seconds=seconds,total_pilot_seconds=seconds+original['seconds'])
    write(OUT/'corrected-pilot.json',p);print(json.dumps(p),flush=True)
    if p['eligible']:write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_cap_seconds=650,
        paths=[[n,r,f'corrected-pilot-{n}-{r}'] for n,r in PATHS],
        bindings={str(f):sha(f) for f in [Path(__file__),Path(matter.__file__),Path(base.__file__),OUT/'angular-repair-plan.json',OUT/'corrected-pilot.json',matter.OUT/'audit.json']}))


def audit():
    p=json.loads((OUT/'result.json').read_text());assert p['passed'];errors=[]
    gamma=1-1/np.sqrt(2)
    for n,r in PATHS:
        d=np.load(OUT/f'steps-{n}-reference-{r}.npz');t=d['accepted_angular_times'];L=d['accepted_angular_luminosity'];assert len(t)==len(L)==2*n
        h=prior.run.flow.old.END/n;expected=np.ravel(h*(np.arange(n)[:,None]+np.array([gamma,1.])))
        assert np.max(abs(t-expected))<1e-18
        # Uniform8-angle bins have weights1/8 and positive centers1/8,3/8,5/8,7/8.
        flux=L@((np.arange(4)*2+1)/32);pairs=flux.reshape(n,2)
        integral=h*np.sum(pairs*np.array([1-gamma,gamma]))
        err=float(abs(integral-d['radial_ports'][-1,1,1])/max(h*np.sum(abs(pairs)),1.));assert err<1e-12;errors.append(err)
    write(OUT/'audit.json',dict(classification='Counterexample candidate',passed=True,actual_updated_material_motion_returned_to_photons=True,
        accepted_angular_port_relative=errors,full_EOS_history_error_enclosed=False,coupled_fixed_point_verified=False,final_charge_solved=False))


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None if sys.argv[5]=='None' else int(sys.argv[5]),None if sys.argv[6]=='None' else sys.argv[6])
    else:globals()[sys.argv[1]]()
