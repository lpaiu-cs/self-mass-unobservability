"""Return corrected-motion photons to the actual material and compact GR."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,resource,sys,time
import numpy as np
import def_native_corrected_motion_feedback as photons
import def_native_corrected_joint_return as previous
import verify_native_corrected_joint_return as execution

OUT=photons.OUT/'return';GR=OUT/'gr';write=photons.write;sha=photons.sha
physical=previous.matter;prior=previous.prior;PATHS=physical.PATHS


def photon_path(n,r):return photons.OUT/f'steps-{n}-reference-{r}.npz'


class Material(previous.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference,steps);p=np.load(photon_path(steps,reference))
        ids=[int(np.argmin(abs(p['t']-t))) for t in self.t];assert np.max(abs(p['t'][ids]-self.t))<1e-18
        c=p['collision_transfer'][ids]
        self.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/self.a,c[:,:,0],c[:,:,1]],axis=1)/physical.AMP
    run=FunctionType(previous.Material.run.__code__,dict(previous.Material.run.__globals__,OUT=OUT),
                     argdefs=previous.Material.run.__defaults__)


# Calibrate actual late RHS cost on saved states, without advancing or altering
# the accepted prefix. Model/point setup is excluded from the warm raw-call cost.
worker_source=inspect.getsource(physical.worker)
anchor="    write(OUT/f'{label}.json',row);print(json.dumps(row),flush=True);assert row['passed'],row"
calibration="""    row['nonstepping_setup_seconds']=row['worker_wall_seconds']-row['seconds']
    if limit is not None:
        sample=np.load(previous.OUT/f'steps-{steps}-reference-{reference}.npz');calibration=[];mark=time.monotonic()
        for k in [8,16]:
            j=min(k,15);m.point(j);m.point(j+1);t=float(m.t[k]);i=int(np.argmin(abs(sample['t']-t)))
            calls=m.raw_calls;start=time.monotonic();m.rhs(t,sample['history_scaled'][i])
            calibration.append(dict(time=t,seconds=time.monotonic()-start,raw_calls=m.raw_calls-calls))
        row['late_warm_calibration']=calibration;row['calibration_seconds']=time.monotonic()-mark
        row['calibration_branch_ratio']=m.physical_branch_ratio
        row['passed']=row['passed'] and m.physical_branch_ratio<.01
        row.update(worker_wall_seconds=time.monotonic()-started,worker_CPU_seconds=time.process_time()-cpu,
                   peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
"""
worker_source=prior.replace(worker_source,anchor,calibration+anchor)
worker_scope=dict(vars(physical),OUT=OUT,Material=Material,previous=previous)
exec(compile(worker_source,__file__,'exec'),worker_scope);worker=worker_scope['worker']
dispatch=FunctionType(execution.dispatch.__code__,dict(vars(execution),OUT=OUT,worker=worker))
source_scope=dict(previous.source_scope,OUT=OUT,GR=GR,Material=Material,
                  photons=SimpleNamespace(path=photon_path),material_path=lambda n,r:OUT/f'steps-{n}-reference-{r}.npz')


def source_worker(n,r,label,limit,restart):
    cap=2*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));physical.configure()
    scope=dict(source_scope,steps_arg=n,ref_arg=r)
    exec(compile(execution.source,__file__,'exec'),scope);scope['sources']()
    p=OUT/f'{label}.json';row=json.loads(p.read_text());assert row['passed']
    row['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024;write(p,row)


source_dispatch=FunctionType(execution.dispatch.__code__,dict(vars(execution),OUT=OUT,worker=source_worker))


class GRResponse(previous.GRResponse):
    run=FunctionType(previous.GRResponse.run.__code__,dict(previous.GRResponse.run.__globals__,OUT=GR))


gr_source=previous.old.source_owner.gr_source.replace('budget_seconds=120','budget_seconds=240').replace('signal.alarm(120)','signal.alarm(240)')
gr_scope=dict(previous.fields_owner.__globals__,OUT=OUT,GR=GR,GRResponse=GRResponse,__file__=__file__)
exec(compile(gr_source,__file__,'exec'),gr_scope);fields_owner=gr_scope['fields']


def prepare():
    assert not OUT.exists();assert json.loads((photons.OUT/'audit.json').read_text())['passed'];OUT.mkdir();GR.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply actual new paired photon transfers to the corrected source-knot material integrator and apply its real pressure/trace to GR; measure the remaining waveform residual and its scale relative to temporal comparisons.',
        reuse='Same corrected physical background/EOS,531 cells,64/128 paths,SSP2,CFL and3.434431ms. Reuse the completed photons. Retain all source fields and the2*T deep radial stress. No physical background replay.',
        budgets=dict(pilot_seconds=120,production_seconds=650,CPU_processes=3,threads_each=1,total_material_virtual_GiB=6,
                     source_dispatch_seconds=180,source_per_worker_seconds=120,GR_seconds=240,production_attempts=1),
        resource_basis='Phase136 material used293.08s under contention;650s allows a2x margin only if fresh prefixes and warm late native-RHS samples support it. Its source/GR setup hit serial caps; use independent forked source readouts and one240s four-path GR readout from the start. These are new declared budgets, not retroactive passes of earlier caps.',
        forecast='Use completed Phase136 raw-call counts times the maximum of fresh early cost and warm midpoint/endpoint raw cost, plus measured model setup and15s. Require2x max below650s. Later state/host changes remain unmeasured.',
        gates=dict(owner=1e-8,conservation=1e-8,directional=.002,branch=.01,small_state=1e-6,time=.02,background=.02,pressure=.002,GR_quadrature=.002,GR_independent=1e-9),
        interpretation='Compare absolute energy/H mismatch to the larger absolute photon/material time difference on the same17 times. This diagnoses numerical scale; it is not an error enclosure, uniform contraction or full fixed-point certificate.',
        stop='One production only; stop on any numerical gate or declared cap. No automatic budget expansion, additional waveform loop, finer mesh/clock/horizon or physical replay.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(photons.__file__),Path(previous.__file__),Path(execution.__file__),photons.OUT/'audit.json',previous.OUT/'audit.json']}))


def pilot():
    assert not (OUT/'pilot.json').exists()
    rows,seconds=dispatch([(n,r,f'pilot-{n}-{r}',2,None) for n,r in PATHS],120)
    old=json.loads((previous.OUT/'production.json').read_text())['paths'];estimates=[]
    for row in rows:
        ref=next(p for p in old if (p['steps'],p['reference'])==(row['steps'],row['reference']))
        rate=max([row['seconds']/row['raw_owner_calls']]+[p['seconds']/p['raw_calls'] for p in row['late_warm_calibration']])
        estimates.append(ref['raw_owner_calls']*rate+row['nonstepping_setup_seconds']+15)
    p=dict(classification='Counterexample candidate',rows=rows,forecast_each_seconds=estimates,upper_seconds=2*max(estimates),
           seconds=seconds,eligible=all(r['passed'] for r in rows) and 2*max(estimates)<650)
    write(OUT/'pilot.json',p);print(json.dumps(p),flush=True)
    if p['eligible']:write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_seconds=650,
        bindings={str(p):sha(p) for p in [Path(__file__),Path(photons.__file__),Path(previous.__file__),Path(execution.__file__),OUT/'plan.json',OUT/'pilot.json']}))


def production():
    assert not (OUT/'production.json').exists();p=json.loads((OUT/'execution-plan.json').read_text());assert p['eligible']
    for f,h in p['bindings'].items():assert sha(f)==h,f
    started=time.monotonic()
    try:rows,seconds=dispatch([(n,r,f'steps-{n}-reference-{r}',None,f'pilot-{n}-{r}') for n,r in PATHS],650)
    except Exception as exc:write(OUT/'production-failure.json',dict(error=repr(exc),seconds=time.monotonic()-started));raise
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),paths=rows,seconds=seconds,
                CPU_seconds=sum(r['worker_CPU_seconds'] for r in rows),sum_peak_RSS_bytes=sum(r['peak_RSS_bytes'] for r in rows),final_charge_solved=False)
    write(OUT/'production.json',result);print(json.dumps(result),flush=True);assert result['passed']


def sources():
    assert json.loads((OUT/'production.json').read_text())['passed'];assert not (OUT/'sources.json').exists()
    write(OUT/'source-plan.json',dict(classification='Counterexample candidate',budget_seconds=180,per_worker_seconds=120,
        claim='Use the unchanged complete conservative source map and actual new photons; checkpoint each pressure verdict before aggregation.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(execution.__file__),Path(previous.old.source_owner.__file__),OUT/'production.json',photons.OUT/'audit.json']}))
    rows,seconds=source_dispatch([(n,r,f'source-path-{n}-{r}',None,None) for n,r in PATHS],180)
    histories=[];stresses=[]
    for n,r in PATHS:
        d=np.load(OUT/f'steps-{n}-reference-{r}.npz');s=np.load(OUT/f'stress-{n}-reference-{r}.npz')
        ids=[int(np.argmin(abs(d['t']-v))) for v in s['t']];assert np.max(abs(d['t'][ids]-s['t']))<1e-18
        histories.append(d['history_scaled'][ids]);stresses.append(np.concatenate([s['material'],np.stack([s['photon_energy'],s['photon_radial_pressure']],axis=1)],axis=1))
    def compare(a,b):return (np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-300)).tolist()
    comparisons=dict(time=compare(histories[0],histories[1]),background=compare(histories[2],histories[1]),stress_time=compare(stresses[0],stresses[1]),stress_background=compare(stresses[2],stresses[1]))
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(v for c in comparisons.values() for v in c)<.02,
        comparisons=comparisons,paths=[r['path'] for r in rows],seconds=seconds,sum_peak_RSS_bytes=sum(r['peak_RSS_bytes'] for r in rows),final_charge_solved=False)
    write(OUT/'sources.json',result);print(json.dumps(result),flush=True);assert result['passed']


def fields():physical.configure();fields_owner()


def audit():
    a=json.loads((OUT/'production.json').read_text());b=json.loads((OUT/'sources.json').read_text());g=json.loads((GR/'result.json').read_text())
    assert a['passed'] and b['passed'] and g['passed'];checked=0;source_error=0.;changes=[];material=[];light=[]
    for plan in [photons.OUT/'plan.json',photons.OUT/'resource-plan.json',photons.OUT/'execution-plan.json',OUT/'plan.json',OUT/'execution-plan.json',OUT/'source-plan.json',GR/'plan.json']:
        for p,h in json.loads(plan.read_text())['bindings'].items():
            target=photons.OUT/'first-producer.py' if plan==photons.OUT/'plan.json' and Path(p)==Path(photons.__file__) else Path(p)
            assert sha(target)==h,p;checked+=1
    for n,r in PATHS:
        d=np.load(GR/f'source-{n}-reference-{r}.npz');s=np.load(OUT/f'stress-{n}-reference-{r}.npz');stress=s['material']
        rest=d['baryon_g'].astype(np.longdouble)*np.longdouble(d['cx'])*np.longdouble(prior.C)**2;total=rest+d['gas_nonrest_energy_erg']
        errors=[total-stress[:,0],d['nonrest_trace_erg']+rest-(stress[:,0]-stress[:,1]-2*stress[:,3]),
                d['nonrest_stress_erg']+rest-(stress[:,0]-stress[:,1]),d['metric_stress_erg']-(total+d['photon_energy_erg']-stress[:,1]-d['photon_radial_pressure_erg'])]
        source_error=max(source_error,float(max(np.max(abs(v)) for v in errors)/max(np.max(abs(stress)),1e-300)))
        old=np.load(previous.OUT/f'steps-{n}-reference-{r}.npz');new=np.load(OUT/f'steps-{n}-reference-{r}.npz');assert np.array_equal(old['t'],new['t'])
        x=old['history_scaled'];y=new['history_scaled'];norm=np.maximum(np.max(np.sum(abs(y),axis=2),axis=0),1e-300)
        changes.append(dict(steps=n,reference=r,times=len(new['t']),relative=(np.max(np.sum(abs(y-x),axis=2),axis=0)/norm).tolist()))
        ids=[int(np.argmin(abs(new['t']-t))) for t in s['t']];material.append(y[ids][:,[2,3]]*physical.AMP)
        p=np.load(photon_path(n,r));ids=[int(np.argmin(abs(p['t']-t))) for t in s['t']];assert np.max(abs(p['t'][ids]-s['t']))<1e-18
        light.append(p['moments'][ids][:,[1,2]])
    assert source_error<1e-12
    norm=lambda x:np.max(np.sum(abs(x),axis=2),axis=0)
    absolute=norm(material[1]-light[1]);material_time=norm(material[0]-material[1]);photon_time=norm(light[0]-light[1])
    scale=np.maximum(material_time,photon_time);before=json.loads((previous.OUT/'audit.json').read_text());charge=g['paths'][0]['endpoint_compact_with_metric']
    primary=float(np.load(prior.run.physical.GR/'wave-128-g8.npz')['free_scalar'][-1])
    result=dict(classification='Counterexample candidate',passed=True,bindings_checked=checked,source_identity_relative=source_error,
        corrected_motion_applied_to_actual_photons=True,new_paired_transfers_applied_to_actual_free_material=True,new_sources_applied_to_compact_GR=True,
        material_history_change=changes,canonical_common_times=17,previous_energy_H_waveform_residual=before['energy_H_waveform_residual'],
        energy_H_waveform_residual=b['paths'][1]['energy_H_waveform_residual'],absolute_energy_H_mismatch=absolute.tolist(),
        absolute_material_time_difference=material_time.tolist(),absolute_photon_time_difference=photon_time.tolist(),
        mismatch_over_larger_temporal_difference=(absolute/np.maximum(scale,1e-300)).tolist(),
        residual_below_temporal_comparison_scale=bool(np.all(absolute<scale)),
        additional_compact_charge=charge,change_from_corrected_prior_charge=charge-before['additional_compact_charge'],additional_over_primary_charge=charge/primary,
        uniform_error_bound=False,uniform_contraction_bound=False,coupled_fixed_point_verified=False,new_GR_reapplied_to_transport=False,
        full_EOS_history_error_enclosed=False,full_exterior_scalar=False,discarded_material_transport_closed=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
