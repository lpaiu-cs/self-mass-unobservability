"""Apply the verified true-native thermal repair to actual coupled stages."""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
import repair_native_thermal_precision as repair

prior=repair.prior.prior;OUT=Path('native-thermal-continuation235-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha;joint=prior.joint
CAPS=dict(prepare=180,check=300,coarse=7200,prefix=1800)


def prepare():
    result=read(repair.OUT/'result.json');assert result['affine_native_B_mismatch']<1e-18
    assert max(r['norm_stage_effect'] for r in result['increment_linearity'])<1e-25
    assert not OUT.exists();OUT.mkdir();files=[]
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [OLD/f'sweep-1/photons/interval-15-{n}{ext}' for n in [64,128] for ext in ['.npz','.json']]
    for src in inputs:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    (OUT/'sweep-1/material').mkdir(exist_ok=True)
    files += [OLD/n for n in ['last-accepted-64.npz','rejected-joint-stage.npz','rejected-joint-stage.json','coarse-receipt.json','failure-64.json']]
    files += [repair.OUT/n for n in ['plan.json','result.json','comparison.npz','check-receipt.json','symbolic.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Resolve actual118and finish the two remaining original coarse substeps using the demonstrated same-equation native arithmetic repair.',
        evidence='233exactly reproduces231linear9.97529e-16and actual1.21893e-6; actual native-minus-affine B accounts for it.234removes internal long-double thermal conversion while preserving its fixed coefficient, reducing that mismatch to2.01209e-22 and half/double local-increment nonlinearity below2.1e-70.40/80digits and original constitutive probes pass. The old proposal still fails the NEW RHS, so no acceptance is inherited.',
        method='Restore117state/history exactly. Use231last rejected pair only as a Newton guide. Build the repaired high-precision tangent for every actual model; its original branch rules evaluate the true nonlinear flux. Recompute original231linear RHS with this owner and use the existing integer solver. Keep all original stage,physical,native,constitutive,balance and time gates.',
        history='Run the required exact-inverse saved-prefix material ledger check only after actual coarse success. Neither that check nor coarse completion admits paired time convergence, uniform earlier vector accuracy, self-GR or final charge. The separately running232fine uses frozen earlier arithmetic and needs explicit cross-evaluation before any paired admission.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,virtual_GiB=8,max_Newton=8,max_linear=12,
        forecast='231failed8proposals263.89s.234checks one reconstructed pair in its bound receipt. Allow2hours for the2coarse original substeps plus30minutes for saved-prefix verification; no accepted-state replay, new clock/grid/period or physical parameter path.',
        stop='Any original numerical/physical/restart/prefix gate or8Newton/12linear/wall cap. Preserve failed proposals and full accepted state.224and232source/plans are unchanged.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(repair.OUT/'symbolic.json'))


def initialize(seed=True):
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(seed)
    constructor=prior.prior.owner.Model.__init__
    def construct(m,n):
        constructor(m,n)
        with repair.precision.mp.workdps(60):m.precise_tangent=repair.build()
        if seed and n==64:
            z=np.load(OLD/'rejected-joint-stage.npz')
            m.resume_seed=dict(time=z['time'][()],solution=z['solution'].copy(),equations=[])
    prior.prior.owner.Model.__init__=construct


def factory(log,n):return FunctionType(prior.factory.__code__,dict(prior.factory.__globals__,OUT=OUT))(log,n)
def check():FunctionType(prior.check.__code__,dict(prior.check.__globals__,OUT=OUT,initialize=initialize))()
def coarse():FunctionType(prior.evolve.__code__,dict(prior.evolve.__globals__,OUT=OUT,initialize=initialize,factory=factory))(64)
def prefix():FunctionType(prior.prefix.__code__,dict(prior.prefix.__globals__,OUT=OUT,initialize=initialize))()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action=='coarse':assert read(OUT/'restart-check.json')['passed']
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
