"""Continue the already planned fine clock while preserving coarse rejection."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import continue_precise_conserved_native as prior
import continue_precise_momentum as dynamics
import continue_postfloor_native as legacy

OUT=Path('native-independent-fine232-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,check=300,fine=10800)


def prepare():
    assert read(OLD/'controller-status.json')['state']=='failed'
    assert 'True native joint Radau equation' in read(OLD/'coarse-receipt.json')['error']
    assert not OUT.exists();OUT.mkdir();files=[]
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [OLD/f'sweep-1/photons/interval-15-{n}{ext}' for n in [64,128] for ext in ['.npz','.json']]
    for src in inputs:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    (OUT/'sweep-1/material').mkdir(exist_ok=True)
    files += [OLD/n for n in ['plan.json','coarse-receipt.json','controller-status.json','rejected-joint-stage.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Finish the already planned independent128fine clock from215accepted steps, without treating coarse failure as a dependency of its differential equations.',
        evidence='231coarse failed actual118after263.89s despite8linear passes. Fine has16remaining original physical substeps and has not run. Continuing this existing path costs less than serially deferring it through further coarse arithmetic repairs.',
        change='Remove ONLY the old scheduling assertion requiring a passed coarse path before fine dispatch. Keep231equations, precise conserved native inputs, integer solver, accepted checkpoint, original fine grid/clock/period and every per-path physical/numerical gate. No fabricated coarse pass and no paired audit or full-history admission.',
        controls='Zero-step fine restart must reproduce every saved physical/history array exactly. Capture actual new accepted photons/ports. Coarse rejection, original paired2percent criterion and required prefix/native revalidation remain open even if this independent path passes.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,virtual_GiB=8,
        forecast='231coarse8proposals263.89s; fine stiffness is unmeasured. Existing16fine substeps only,3hour cap and8Newton/12linear limits. No accepted-prefix replay, finer grid, longer duration or additional physical parameter path.',
        stop='Any exact restart,original linear/nonlinear/constitutive/balance gate or3hour cap. Save last accepted state before every step. A fine success is not paired convergence or final charge.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},coarse_rejection_preserved=True,paired_time_admitted=False,final_charge_conclusion='unadjudicated',full_goal_complete=False))


def initialize(seed=True):FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(seed)
def factory(log,n):return FunctionType(prior.factory.__code__,dict(prior.factory.__globals__,OUT=OUT))(log,n)


def check():
    initialize(False);m=prior.prior.owner.Model(128);row=m.run(128,'resume-check-128',120,restart='interval-15-128')
    a=np.load(OUT/'sweep-1/photons/interval-15-128.npz');b=np.load(OUT/'sweep-1/photons/resume-check-128.npz')
    assert set(a.files)==set(b.files)
    for key in a.files:assert np.array_equal(a[key],b[key]),key
    assert row['actual_new_steps']==row['new_steps']==0 and row['passed']
    write(OUT/'restart-check.json',dict(classification='Counterexample candidate',passed=True,all_saved_arrays_exact=True,accepted_steps=215,new_physical_steps=0))


def fine():
    assert read(OUT/'restart-check.json')['passed'];dynamics.prior.stable.stable_operator=dynamics.stable_operator
    base=inspect.getsource(legacy.before.evolve);guard="    if n==128:assert read(OUT/'path-64.json')['passed']\n"
    assert base.count(guard)==1;base=base.replace(guard,'')
    source=inspect.getsource(legacy.evolve)
    old='factory=FunctionType(prior.factory.__code__,dict(prior.factory.__globals__,OUT=OUT))'
    assert source.count(old)==1;source=source.replace(old,'factory=current_factory')
    source=source.replace('source=inspect.getsource(before.evolve)','source=base_source')
    ns=dict(legacy.evolve.__globals__,OUT=OUT,initialize=initialize,current_factory=factory,base_source=base)
    exec(compile(source,__file__,'exec'),ns);ns['evolve'](128)
    row=read(OUT/'path-128.json');row.update(coarse_rejection_preserved=True,paired_time_admitted=False,prefix_revalidation_pending=True,full_goal_complete=False)
    write(OUT/'result.json',row)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
