"""Finish the capped source readout from completed corrected trajectories."""
from pathlib import Path
from types import FunctionType
import json,multiprocessing,resource,shutil,signal,sys,time
import numpy as np
import def_native_corrected_joint_return as run

OUT=run.OUT;GR=run.GR;REC=OUT/'source-recovery';write=run.write;sha=run.sha

# The original transformation is unchanged. Checkpoint each independent path
# before aggregation so a timeout cannot lose completed pressure-probe verdicts.
source=run.old.source_owner.source
source=run.prior.replace(source,'signal.alarm(75)','signal.alarm(120)')
source=run.prior.replace(source,'for steps,ref in [[64,128],[128,128],[128,64]]:',
                         'for steps,ref in [(steps_arg,ref_arg)]:')
source=source[:source.index('    def compare(a,b):')]
source+="    write(OUT/f'source-path-{steps_arg}-{ref_arg}.json',dict(classification='Counterexample candidate',passed=all(r['conservation']<1e-8 and r['pressure_probe']<.002 for r in rows),path=rows[0],seconds=time.monotonic()-start))\n    signal.alarm(0)\n"


def dispatch(specs,cap):
    # Linux fork shares the already imported, frozen owners. New interpreters
    # otherwise reread hundreds of files from the contended Windows mount.
    started=time.monotonic();children=[]
    try:
        for spec in specs:
            p=multiprocessing.get_context('fork').Process(target=worker,args=spec)
            p.start();children.append((p,spec[2]))
        while any(p.is_alive() for p,_ in children):
            if time.monotonic()-started>cap:raise TimeoutError('Registered source dispatch cap')
            for p,label in children:
                if p.exitcode not in [None,0]:raise RuntimeError((label,p.exitcode))
            time.sleep(.25)
        assert all(p.exitcode==0 for p,_ in children)
        return [json.loads((OUT/f'{label}.json').read_text()) for _,label in children],time.monotonic()-started
    finally:
        for p,_ in children:
            if p.is_alive():p.terminate()
        for p,_ in children:p.join()


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


def register_fork():
    assert not any(OUT.glob('source-path-*.json'))
    plan=json.loads((REC/'pre-dispatch-plan.json').read_text())
    plan['pre_dispatch_resource_observation']='Before readout dispatch, the import-only preparer took over120s wall but about4.5s CPU. Live process inspection found a separate VLLM engine using about457percent CPU and another active workload. Those unrelated jobs are not modified. The early quiet-host timing forecast is not reliable under this contention.'
    plan['dispatch_repair']='Use Linux fork after owner imports, preserving each independent process model and120s source computation alarm;180s dispatcher includes child setup and all computation. This removes three redundant full interpreter/import paths. No source worker has run under the pre-dispatch plan.'
    plan['budgets']['worker_source_seconds']=120
    plan['stop']='One forked dispatch only;120s per source computation and180s total. Preserve any completed per-path checks on failure. No automatic subsequent extension or physical trajectory replay.'
    for p in [Path(__file__),REC/'pre-dispatch-plan.json',REC/'pre-dispatch-producer.py']:plan['bindings'][str(p)]=sha(p)
    write(REC/'plan.json',plan)


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
    for p,h in plan['bindings'].items():
        target=REC/'source-producer.py' if Path(p)==Path(__file__) else Path(p)
        assert sha(target)==h,p;checked+=1
    extra=GR/'resume-plan.json'
    if extra.exists():
        for p,h in json.loads(extra.read_text())['bindings'].items():assert sha(p)==h,p;checked+=1
    run.audit();p=OUT/'audit.json';r=json.loads(p.read_text())
    r.update(readout_recovery_bindings_checked=checked,original_serial_source_budget_verdict=False,
             completed_material_trajectories_replayed=0,original_coarse_source_arrays_identical=True)
    write(p,r)


def resume_fields():
    assert not (GR/'result.json').exists() and not (GR/'resume-plan.json').exists()
    names=['128-reference-128-g8','128-reference-128-g4']
    saved=[GR/f'fields-{name}{suffix}' for name in names for suffix in ['.json','.npz']]
    assert all(p.exists() for p in saved)
    write(GR/'original-failure.json',dict(classification='Counterexample candidate',passed=False,
        error='Original120s GR alarm interrupted the coarse-time potential propagation',
        completed_paths=names,remaining_paths=['64-reference-128-g8','128-reference-64-g8'],
        independent_direct_pending=True))
    write(GR/'resume-plan.json',dict(classification='Counterexample candidate',
        decision='Complete only the two missing registered GR controls and independent direct comparison. The fine8/fine4 fields and all physical histories are reused byte-for-byte. No extra path or resolution.',
        reassessment='Completed fine rules took17.87s and14.88s, while the original120s timer also included a slow model initialization. Two remaining propagations plus direct comparison are estimated40-60s after setup; allow one additional120s readout cap including setup. This is an explicit budget revision under measured host contention; the original120s verdict remains failed.',
        budgets=dict(original_failed_seconds=120,remaining_readout_seconds=120,CPU_processes=1,threads=1,new_physical_trajectory_steps=0),
        stop='One resume only. Stop on120s or any original scientific comparison gate; retain all completed fields. No further automatic readout or refinement.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),GR/'plan.json',GR/'original-failure.json',OUT/'sources.json',*saved]}))
    run.matter.configure();started=time.monotonic()
    signal.signal(signal.SIGALRM,run.matter.branch.base.flow.old.optical.timeout);signal.alarm(120)
    try:
        m=run.GRResponse();rows=[json.loads((GR/f'fields-{name}.json').read_text()) for name in names]
        for label in ['64-reference-128','128-reference-64']:rows.append(m.run(label,8))
        fine=np.load(GR/'fields-128-reference-128-g8.npz');norm=max(np.max(abs(fine['U'])),1e-300)
        comparisons={key:float(np.max(abs(np.load(GR/f'fields-{name}.npz')['U']-fine['U']))/norm)
            for key,name in [('quadrature','128-reference-128-g4'),('time','64-reference-128-g8'),('background','128-reference-64-g8')]}
        data=dict(np.load(GR/'source-128-reference-128.npz'))
        direct,error=run.matter.previous.independent.direct(m,data,8);agreement=abs(direct/rows[0]['endpoint_direct']-1)
        r=dict(classification='Counterexample candidate',passed=comparisons['quadrature']<.002 and max(comparisons['time'],comparisons['background'])<.02 and agreement<1e-9,
               comparisons=comparisons,independent_direct_relative=agreement,paths=rows,seconds=time.monotonic()-started,prior_failed_seconds=120,
               returned_source_applied_to_represented_GR=True,coupled_fixed_point_verified=False,final_charge_solved=False)
        write(GR/'result.json',r);print(json.dumps(r),flush=True);assert r['passed']
    except Exception as exc:
        write(GR/'resume-failure.json',dict(classification='Counterexample candidate',passed=False,error=repr(exc),seconds=time.monotonic()-started));raise
    finally:signal.alarm(0)


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None,None)
    else:globals()[sys.argv[1]]()
