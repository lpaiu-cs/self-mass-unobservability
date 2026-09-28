"""Finish the capped source readout from completed corrected trajectories."""
from pathlib import Path
from types import FunctionType
import json,resource,shutil,signal,sys,time
import numpy as np
import def_native_corrected_joint_return as run

OUT=run.OUT;GR=run.GR;REC=OUT/'source-recovery';write=run.write;sha=run.sha

# The original transformation is unchanged. Checkpoint each independent path
# before aggregation so a timeout cannot lose completed pressure-probe verdicts.
source=run.old.source_owner.source
source=run.prior.replace(source,'for steps,ref in [[64,128],[128,128],[128,64]]:',
                         'for steps,ref in [(steps_arg,ref_arg)]:')
source=source[:source.index('    def compare(a,b):')]
source+="    write(OUT/f'source-path-{steps_arg}-{ref_arg}.json',dict(classification='Counterexample candidate',passed=all(r['conservation']<1e-8 and r['pressure_probe']<.002 for r in rows),path=rows[0],seconds=time.monotonic()-start))\n    signal.alarm(0)\n"
dispatch=FunctionType(run.dispatch.__code__,dict(run.dispatch.__globals__,OUT=OUT,__file__=__file__))


def prepare():
    assert not REC.exists() and not (OUT/'sources.json').exists();REC.mkdir()
    assert json.loads((OUT/'production.json').read_text())['passed']
    originals=[GR/'source-64-reference-128.npz',OUT/'stress-64-reference-128.npz',OUT/'source-plan.json']
    for p in originals:shutil.copyfile(p,REC/p.name)
    write(REC/'original-failure.json',dict(classification='Counterexample candidate',passed=False,
        error='TimeoutError after the original75s source wall cap',
        interrupted_location='Constructing the second Material model; numpy zip input read',
        completed_exports=['source-64-reference-128.npz','stress-64-reference-128.npz'],
        missing='Per-path pressure-probe summary and final cross-path comparisons were not checkpointed.',
        physical_material_production_replayed=False))
    write(REC/'plan.json',dict(classification='Counterexample candidate',
        decision='The actual corrected material evolution is complete, but its new pressure/trace must reach GR to judge the charge change. Stop and preserve the75s serial failure; change only readout scheduling and checkpoint each completed path.',
        reuse='Use all completed material/photon/background arrays and the exact original source/pressure arithmetic. No trajectory, EOS bank, grid, clock or physical production rerun. Preserve the first coarse export and require exact array equality on its independently completed readout.',
        reassessment='One coarse export used almost the entire original75s, so another unchanged serial75s run is ineligible. Three source paths are independent. One parallel readout has a180s dispatcher cap and unchanged75s computation alarm per path,3 CPU processes with2GiB virtual memory each. This is an explicit readout-budget revision after failure, not a successful original75s verdict.',
        forecast='Observed coarse source export near75s includes GR setup. Up to75s concurrent path compute plus measured import/setup overhead suggests about90s; allow180s dispatcher including setup. This is an extrapolation under shared file I/O, not a guaranteed runtime. No further automatic extension.',
        budgets=dict(original_failed_source_seconds=75,parallel_dispatch_seconds=180,worker_source_seconds=75,CPU_processes=3,threads_each=1,total_virtual_GiB=6,readout_attempts=1,GR_seconds=120),
        stop='Stop on any worker/source gate,75s per-worker or180s dispatcher cap. Do not apply failed sources to GR. No extra physical production.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(run.old.source_owner.__file__),OUT/'plan.json',OUT/'production.json',REC/'original-failure.json',*list(REC.glob('*.npz'))]}))


def worker(steps,reference,label,limit,restart):
    cap=2*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));run.matter.configure()
    scope=dict(run.source_scope,steps_arg=steps,ref_arg=reference)
    exec(compile(source,__file__,'exec'),scope);scope['sources']()
    p=OUT/f'{label}.json';r=json.loads(p.read_text());assert r['passed']
    r['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024;write(p,r)


def finish():
    assert not (OUT/'sources.json').exists()
    plan=json.loads((REC/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    started=time.monotonic()
    try:rows,seconds=dispatch([(n,r,f'source-path-{n}-{r}',None,None) for n,r in run.matter.PATHS],180)
    except Exception as exc:
        write(REC/'parallel-failure.json',dict(classification='Counterexample candidate',passed=False,error=repr(exc),seconds=time.monotonic()-started));raise
    for name,folder in [('source-64-reference-128.npz',GR),('stress-64-reference-128.npz',OUT)]:
        a=np.load(REC/name);b=np.load(folder/name)
        assert a.files==b.files and all(np.array_equal(a[k],b[k]) for k in a.files),name
    histories=[];stresses=[];path_rows=[]
    for (n,r),row in zip(run.matter.PATHS,rows):
        d=np.load(OUT/f'steps-{n}-reference-{r}.npz');s=np.load(OUT/f'stress-{n}-reference-{r}.npz');t=s['t']
        ids=[int(np.argmin(abs(d['t']-v))) for v in t];assert np.max(abs(d['t'][ids]-t))<1e-18
        histories.append(d['history_scaled'][ids]);stresses.append(np.concatenate([s['material'],np.stack([s['photon_energy'],s['photon_radial_pressure']],axis=1)],axis=1));path_rows.append(row['path'])
    def compare(a,b):return (np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-300)).tolist()
    comparisons=dict(time=compare(histories[0],histories[1]),background=compare(histories[2],histories[1]),stress_time=compare(stresses[0],stresses[1]),stress_background=compare(stresses[2],stresses[1]))
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(v for c in comparisons.values() for v in c)<.02,
        comparisons=comparisons,paths=path_rows,seconds=seconds,prior_failed_source_seconds=75,
        sum_peak_RSS_bytes=sum(r['peak_RSS_bytes'] for r in rows),original_coarse_export_arrays_identical=True,
        returned_photon_transfer_applied_to_material=True,coupled_fixed_point_verified=False,final_charge_solved=False)
    write(OUT/'sources.json',result);print(json.dumps(result),flush=True);assert result['passed']


def audit():
    plan=json.loads((REC/'plan.json').read_text());checked=0
    for p,h in plan['bindings'].items():assert sha(p)==h,p;checked+=1
    run.audit();p=OUT/'audit.json';r=json.loads(p.read_text())
    r.update(readout_recovery_bindings_checked=checked,original_serial_source_budget_verdict=False,
             completed_material_trajectories_replayed=0,original_coarse_source_arrays_identical=True)
    write(p,r)


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None,None)
    else:globals()[sys.argv[1]]()
