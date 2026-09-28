"""Complete the original fine clock with the coarse path's verified arithmetic."""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import continue_true_momentum as prior
import continue_independent_fine_native as independent

OUT=Path('native-common-arithmetic239-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha;joint=prior.joint
CAPS=dict(prepare=300,check=300,fine=21600,audit=600)


def bind(fn,**values):
    return FunctionType(fn.__code__,dict(fn.__globals__,**values),argdefs=fn.__defaults__)


def prepare():
    coarse=read(OLD/'path-64.json');assert coarse['passed'] and coarse['actual_completed_steps']==119
    assert read(OLD/'coarse-receipt.json')['error'] is None;assert not OUT.exists();OUT.mkdir();files=[]
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','path-64.json']]
    inputs += [OLD/f'sweep-1/photons/{name}{ext}' for name in ['interval-15-128','complete-64'] for ext in ['.npz','.json']]
    for src in inputs:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    (OUT/'sweep-1/material').mkdir(exist_ok=True)
    files += [OLD/n for n in ['plan.json','coarse-receipt.json','stage-progress-64.json','restart-check.json']]
    files += [prior.repair.OUT/n for n in ['result.json','check-receipt.json','symbolic.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Complete the SAME128fine clock under the B/thermal/S arithmetic that actually completed coarse119, then use both actual paths in the unchanged10channel time comparison.',
        evidence='238last actual coarse stage passed9.16708e-15in116.38s including setup and full-path readout.232is still computing under the previous arithmetic and cannot certify this corrected pair without further work. Starting the required corrected remaining16fine stages avoids serially waiting for that older arithmetic.',
        reuse='Reuse the exact canonical215accepted-stage fine checkpoint and all its histories; zero-step restart must reproduce EVERY NPZarray. Do not replay that prefix or coarse119. The few232post-prefix stages are under different arithmetic and are retained as separate evidence, not spliced into this corrected solve.',
        method='Reuse238true B/thermal/S native arithmetic, faces/gravity/ledgers and original integer linear solver. Apply to all remaining actual fine stages and capture their actual photons/ports. Use the real238coarse completion, never a fabricated scheduling pass.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,CPU_affinity=10,virtual_GiB=12,max_Newton=8,max_linear=12,
        forecast='238one actual coarse step116.38s including setup;16fine stages would take about31minutes if that cost held. Fine/later branch costs are unmeasured;32..64minutes is a planning estimate only. Allow6hours so necessary work is not repeatedly stopped. One existing physical clock/period, no grid or accuracy relaxation.',
        comparison='Wait for238same-native saved-prefix material-ledger checks before paired admission. These checks do not certify all earlier vector equations uniformly. Preserve the original early-time/global failures, selfGR, EOS/derivative/spatial/boundary/nonlinear/observational/infinity scope.',
        stop='Exact-restart,original linear/nonlinear/physical/native/constitutive/conservation or paired-time failure,8Newton/12linear or cap. Preserve complete accepted state and rejected proposals.232,236and238sources/plans remain frozen.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(prior.repair.OUT/'symbolic.json'))


def initialize(seed=True):bind(prior.initialize,OUT=OUT)(seed)
def factory(log,n):return bind(prior.factory,OUT=OUT)(log,n)
def check():bind(independent.check,OUT=OUT,initialize=initialize)()


def fine():
    assert read(OUT/'path-64.json')['passed']
    bind(independent.fine,OUT=OUT,initialize=initialize,factory=factory)()
    row=read(OUT/'result.json');row.update(actual_corrected_coarse_path_passed=True,same_corrected_native_arithmetic=True)
    write(OUT/'result.json',row)


def audit():
    assert read(OLD/'prefix-result.json')['passed'] and read(OLD/'prefix-receipt.json')['error'] is None
    assert read(OLD/'prefix-receipt.json')['source_sha256']==sha(prior.__file__)
    src=OLD/'prefix-result.json';os.link(src,OUT/'prefix-result.json')
    write(OUT/'prefix-admission.json',dict(classification='Counterexample candidate',passed=True,
        same_coarse118_and_fine215_saved_prefix=True,source=str(src),sha256=sha(src),uniform_prefix_vector_certificate=False))
    bind(prior.prior.prior.audit,OUT=OUT)()
    r=read(OUT/'result.json');r.update(same_corrected_native_arithmetic=True,common_saved_prefix_material_ledger_revalidated=True,
        uniform_prefix_vector_certificate=False,full_goal_complete=False,final_charge_conclusion='unadjudicated')
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_AS,(12*1024**3,12*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
