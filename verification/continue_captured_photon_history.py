"""Resume same-history recovery from an exact original coupled-stage capture."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import resume_coordinate_exact_photons as prior
import capture_original_joint_photons as capture

OUT=Path('native-captured-photon223-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,coarse=7200,fine=10800,audit=300)


def seed_path(n):return OUT/f'input-{n}.npz'


def prepare():
    assert not OUT.exists();OUT.mkdir();files=[];paths=[]
    result=read(capture.OUT/'result.json');assert result['passed'] and result['all_physical_array_values_exact']
    assert read(capture.OUT/'physical-replay.json')['passed'] and result['snapshot']['passed']
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-64','clock-128']:(OUT/folder).mkdir(parents=True)
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        os.link(src,OUT/src.relative_to(OLD));files.append(src)
    for n,src in [(64,OLD/'accepted-64.npz'),(128,capture.OUT/'accepted-128.npz')]:
        os.link(src,seed_path(n));p=np.load(seed_path(n));step=int(p['step']);assert step=={64:11,128:16}[n]
        z=np.load(prior.prior.saved(n));assert p['t']==z['actual_step_edges'][step]
        assert len(json.loads(str(p['logs'])))==step
        for key in ['moments','collisions','ports','packets']:assert len(p[key])==2*step
        paths.append(dict(clock=n,reused_steps=step,remaining=len(z['actual_step_edges'])-1-step))
        files += [src,seed_path(n),prior.prior.saved(n)]
    files += [capture.OUT/n for n in ['result.json','physical-replay.json','fine-receipt.json','symbolic.json']]
    files += [OLD/n for n in ['restart-check.json','resume-controller-status.json','clock-128/snapshot-02.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the actual original coupled-stage capture to unblock the failed fine photon history, then complete the already accepted common15/16trajectory under all original gates.',
        repair='Use222fine16-step checkpoint whose material,photon,native,collision,port and internal endpoint arrays exactly reproduce184. Reuse214coarse11accepted steps. Every new conditional pair keeps214exact coordinate inversion, original model lifetime, whole equation/native identity, endpoint, ledger and angular/radial gates. Original214/220failures remain.',
        no_replay='No accepted photon prefix or native-fluid step is replayed here.222replayed only8old fine steps to obtain the original missing photon data; this is not a new physics path or longer physical interval.',
        persistence='Keep immutable input checkpoints. Write accepted photon/history states before and after every conditional step. Additionally preserve a failed snapshot pair immediately, avoiding the lost-pair recomputation needed in219.',
        gates=read(prior.prior.OUT/'plan.json')['gates'],paths=paths,budgets=CAPS,CPU_threads_per_path=1,virtual_GiB_per_path=6,maximum_parallel_paths=2,
        forecast='220conditional8steps250.96s including setup. At20..45s per missing step, coarse100steps about33..75minutes and fine199about66..149minutes. Later costs unmeasured; allow2/3hours. Existing main218source and plan are unchanged.',
        stop='First original gate or wall failure cancels the other owned recovery and retains both accepted checkpoints and failed candidate. No automatic longer interval, tighter tolerance, new grid or full-native replay.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'restart-check.json',dict(classification='Counterexample candidate',passed=True,paths=paths,
        captured_interval_original_physical_arrays_exact=True,prior652_native_coordinate_check=read(OLD/'restart-check.json')['passed'],inherited_prefix_not_reaudited=True,new_physical_steps=0))


def run(n):
    base_snapshot=FunctionType(prior.prior.snapshot.__code__,dict(prior.prior.snapshot.__globals__,OUT=OUT))
    def capture_snapshot(n,m,z,step,photons,collisions,ports,packets,*,initial_x,photon_pair,gas,t,h):
        try:return base_snapshot(n,m,z,step,photons,collisions,ports,packets)
        except BaseException:
            np.savez_compressed(OUT/f'rejected-snapshot-{n}.npz',step=step,time=t,step_size=h,x_initial=initial_x,
                photon_stage_solution=photon_pair,gas=gas,stage_times=z['joint_stage_times'][2*step:2*step+2],
                collisions=collisions,ports=ports,packets=packets)
            raise
    source=inspect.getsource(prior.run)
    old="snapshot=FunctionType(prior.snapshot.__code__,dict(prior.snapshot.__globals__,OUT=OUT))"
    assert source.count(old)==1;source=source.replace(old,'snapshot=capture_snapshot')
    mark="    exec(compile(source,__file__,'exec'),ns);"
    a="row['snapshot']=snapshot(n,m,z,step,pairs[-1],collisions,ports,packets)"
    b="row['snapshot']=snapshot(n,m,z,step,pairs[-1],collisions,ports,packets,initial_x=x,photon_pair=pairs,gas=gas,t=t,h=h)"
    assert source.count(mark)==1;source=source.replace(mark,f'    assert source.count({a!r})==1\n    source=source.replace({a!r},{b!r})\n'+mark)
    seed=FunctionType(prior.seed.__code__,dict(prior.seed.__globals__,seed_path=seed_path))
    ns=dict(prior.run.__globals__,OUT=OUT,RESUME=True,seed=seed,seed_path=seed_path,capture_snapshot=capture_snapshot)
    exec(compile(source,__file__,'exec'),ns);ns['run'](n)


def audit():
    FunctionType(prior.audit.__code__,dict(prior.audit.__globals__,OUT=OUT))()
    result=read(OUT/'result.json');result.update(original_coupled_capture_applied=True,source_capture_sha256=sha(capture.OUT/'result.json'))
    write(OUT/'result.json',result)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.base.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in ['coarse','fine']:run(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
