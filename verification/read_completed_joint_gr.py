"""Counterexample candidate: read GR from a completed, accepted joint solution.

No new photon/material evolution. The input must pass its whole-period paired
time and same-solution ledger checks before this consumer can be prepared.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json,os,resource,sys,time
import read_full_incident_joint_gr as base

OUT=Path('native-complete-joint-gr-work')
read,write,sha=base.read,base.write,base.sha
CAPS=dict(prepare=60,source=600,fields=300)


def saved(n):return INPUT/f'sweep-1/photons/complete-{n}.npz'


def prepare():
    result=read(INPUT/'result.json')
    assert result['passed'] and result['full_horizon_photon_material_completed']
    assert read(INPUT/'controller-status.json')['state']=='completed'
    assert not OUT.exists();OUT.mkdir()
    for folder in ['sweep-0','sweep-1/photons','sweep-1/material','gr']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    files=list((INPUT/'sweep-0').rglob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    reused={}
    for p in files:
        dst=OUT/p.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [INPUT/n for n in ['result.json','plan.json','path-64.json','path-128.json','coarse-receipt.json','fine-receipt.json','audit-receipt.json']]
    files += [saved(n) for n in [64,128]]+[INPUT/f'sweep-1/photons/interval-15-{n}.npz' for n in [64,128]]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Use the full accepted same joint solution, all17canonical states and its own energy/boundary histories to construct actual material/photon GR sources and the retarded compact field.',
        decision='Passing source identity, original2percent time,0.2percent quadrature and independent GR controls admits using this same field in the next physical GR-return solve. Compact GR is not the final infinity-normalized scalar charge.',
        reuse='Extend the already checked186source/GR mapping over the saved complete period. No new physical steps, EOS/background solve, grid, source waveform or earlier charge is added.',
        precision='Recover Etilde directly from stored normalized gas. Preserve long-double native primitive variations and radial stress in the readout; retain original EOS bank and current inventory displacement. The original full-solution pressure and conserved-moment identities are checked.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='The186three-time source check used about20..30s per attempt; prior full17-time GR field work took12.43s. This reader uses2models and17saved states. Allow10minutes source and5minutes fields, with no physical reintegration.',
        gates=dict(time=.02,quadrature=.002,pressure=.002,identity=1e-12,independent_GR=1e-9,conservation=1e-8),
        input_directory=str(INPUT),full_horizon_seconds=result['full_horizon_seconds'],
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as sp
    a,a0,c,b,e=sp.symbols('a a0 c b e')
    assert sp.simplify((e+(a-a0)*c*b+a0*c*b)/a-c*b-e/a)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Stable nonrest-energy identity only; not an EOS, self-GR or final-charge certificate.'))


def pressure_function():
    fn=FunctionType(base.pressure_function.__code__,dict(base.pressure_function.__globals__,OUT=OUT))()
    source=(OUT/'expanded-pressure.py').read_text()
    source=source.replace('np.zeros(m.n)','np.zeros(m.n,dtype=LD)').replace(',float)',',LD)')
    ns=dict(fn.__globals__);exec(compile(source,__file__,'exec'),ns)
    (OUT/'expanded-pressure.py').write_text(source)
    return ns['primitive']


def source():
    def verify(file,unused):
        n=int(file.stem.rsplit('-',1)[1])
        return base.run.verify(file,INPUT/f'sweep-1/photons/interval-15-{n}.npz')
    run=SimpleNamespace(**dict(vars(base.run),verify=verify))
    fn=FunctionType(base.source.__code__,dict(base.source.__globals__,OUT=OUT,INPUT=INPUT,INTERVAL=16,
        saved=saved,run=run,pressure_function=pressure_function));fn()


def fields():
    fn=FunctionType(base.fields.__code__,dict(base.fields.__globals__,OUT=OUT));fn()
    result=read(OUT/'fields.json')
    assert abs(result['physical_horizon_seconds']-read(OUT/'plan.json')['full_horizon_seconds'])<1e-18
    result.update(full_horizon_completed=True,full_declared_input_period=True,
        self_GR_return_closed=False,physical_final_charge_solved=False)
    write(OUT/'fields.json',result)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    INPUT=Path(sys.argv[2]) if action=='prepare' else Path(read(OUT/'plan.json')['input_directory'])
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));base.run.owner.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
