"""Apply the checked conserved-map arithmetic in actual native evolution."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import repair_native_conserved_precision as repair
import read_complete_radau_history as bridge

prior=repair.prior;flux=repair.flux;OUT=Path('native-conserved-continuation231-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha;joint=prior.joint
CAPS=dict(prepare=180,check=300,coarse=7200,fine=10800,prefix=1800,audit=180)


def prepare():
    result=read(repair.OUT/'result.json');assert result['controls_passed'] and result['arithmetic_effect']>1e-12
    assert not OUT.exists();OUT.mkdir();files=[]
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [OLD/f'sweep-1/photons/interval-15-{n}{ext}' for n in [64,128] for ext in ['.npz','.json']]
    for src in inputs:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    (OUT/'sweep-1/material').mkdir(exist_ok=True)
    files += [OLD/n for n in ['last-accepted-64.npz','rejected-joint-stage.npz','rejected-joint-stage.json','failure-64.json','coarse-receipt.json']]
    files += [repair.OUT/n for n in ['plan.json','result.json','comparison.npz','check-receipt.json','symbolic.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply independently checked high-precision conserved inputs to the true native RHS AND actual nonlinear residual of the remaining same coupled period.',
        evidence='229eight linear passes still fail nonlinear118.230reproduces that full defect exactly. Promoting the identical gas-to-conserved map changes the actual defect by9.67826e-7;40/80digits agree and half/double constitutive probes are below2.24e-27. The old proposal remains rejected at9.68446e-7.',
        method='Restore the same117accepted checkpoint exactly. Use229last rejected pair as a NEW Newton guide only. Recompute the RHS with the promoted native owner; never reuse228linear solution as though its RHS were unchanged. Keep229joint integer solver,8Newton/12linear limits and every original gate. Capture actual accepted photon moments and original boundary returns.',
        history='Coarse/fine calculations are provisional until the changed native arithmetic is checked against the saved prefix material ledger using the existing exact conserved inverse. Run this required prefix check after successful actual continuation, before paired audit/admission. Prefix balance is not a uniform recertification of every earlier vector equation.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,virtual_GiB=8,
        forecast='229eight proposals201.8s;230native controls34.9s including setup. Allow2/3hours for2coarse/16fine steps,30minutes prefix control, with unchanged physical grid and time span. Generous caps, not predicted completion times.',
        stop='Any original exact-restart,linear,nonlinear,constitutive,balance,prefix or paired-time gate;8Newton/12linear or wall cap. No automatic resolution expansion.227and224 remain frozen.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(repair.OUT/'symbolic.json'))


def initialize(seed=True):
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(seed)
    repair.precision.native_B=repair.native_function()
    if seed:
        constructor=prior.owner.Model.__init__
        def construct(m,n):
            constructor(m,n)
            if n==64:
                z=np.load(OLD/'rejected-joint-stage.npz')
                m.resume_seed=dict(time=z['time'][()],solution=z['solution'].copy(),equations=[])
        prior.owner.Model.__init__=construct


def factory(log,n):
    source=inspect.getsource(prior.factory);assert source.count('return seeded')==1;source=source.replace('return seeded','return solve')
    ns=dict(prior.factory.__globals__,OUT=OUT);exec(compile(source,__file__,'exec'),ns);return ns['factory'](log,n)


def check():FunctionType(prior.check.__code__,dict(prior.check.__globals__,OUT=OUT,initialize=initialize))()
def evolve(n):FunctionType(prior.evolve.__code__,dict(prior.evolve.__globals__,OUT=OUT,initialize=initialize,factory=factory))(n)


def prefix():
    source=inspect.getsource(flux.prefix)
    old="g=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])"
    assert source.count(old)==1;source=source.replace(old,'g=restored_gas(m,q)')
    ns=dict(flux.prefix.__globals__,OUT=OUT,OLD=OLD,initialize=lambda seed=False:initialize(False),restored_gas=bridge.recovery.prior.restored_gas)
    exec(compile(source,__file__,'exec'),ns);ns['prefix']()


def audit():
    assert read(OUT/'prefix-result.json')['passed']
    FunctionType(prior.audit.__code__,dict(prior.audit.__globals__,OUT=OUT))()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
