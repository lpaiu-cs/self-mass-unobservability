"""Return the accepted common15/16 GR history at the original actual stages."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import apply_actual_stage_gr as prior

OUT=Path('native-complete-return236-work');INPUT=prior.bridge.OUT
bridge,geometry=prior.bridge,prior.geometry
read,write,sha=prior.read,prior.write,prior.sha
saved=bridge.saved
CAPS=dict(prepare=300,field1288=10800,field648=10800,field1284=10800,
          collect=300,metric=3600,coarse=10800,fine=21600,audit=600)
FIELDS=dict(field1288=(128,8),field648=(64,8),field1284=(128,4))


def bind(fn,**values):
    return FunctionType(fn.__code__,dict(fn.__globals__,**values),argdefs=fn.__defaults__)


representation=bind(prior.representation,INPUT=INPUT)
model=SimpleNamespace(**dict(vars(prior.prior),saved=saved))
model.initialize=bind(prior.prior.initialize,OUT=OUT,saved=saved)


class Response(bridge.Response):
    setup=bind(prior.Response.setup,INPUT=INPUT)


class Lapse(prior.Lapse):
    boundary=bind(prior.Lapse.boundary,prior=model)


def prepare():
    full=read(INPUT/'result.json');assert full['GR_return_admitted'] and full['accepted_common_history_applied']
    assert read(prior.OUT/'result.json')['passed'];assert not OUT.exists();OUT.mkdir();files=[]
    for part in ['sweep-1/photons','sweep-1/material','gr','metric']:(OUT/part).mkdir(parents=True)
    for src in list((INPUT/'sweep-0').rglob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    candidates=[0.];steps=[]
    for n in [64,128]:
        p=np.load(saved(n));steps.append(len(p['actual_step_edges'])-1)
        candidates+=list(p['joint_stage_times'])+list(p['actual_step_edges']);files.append(saved(n))
    clock=[]
    for t in sorted(candidates):
        if not clock or t-clock[-1]>1e-18:clock.append(t)
    clock=np.array(clock);assert steps==[111,215]
    assert abs(clock[-1]-full['physical_horizon_seconds'])<1e-18
    checks=[]
    for n in [64,128]:
        raw,poly=representation(n);d=dict(raw,t=clock,original_clock=np.array(n));d.update({k:f(clock) for k,f in poly.items()})
        ids=np.array([np.argmin(abs(clock-t)) for t in raw['t']]);assert np.max(abs(clock[ids]-raw['t']))<1e-18
        jumps={k:float(np.max(abs(d[k][ids]-raw[k]))/max(np.max(abs(raw[k])),prior.LD('1e-290'))) for k in poly}
        for k in poly:d[k][ids]=raw[k]
        assert all(np.array_equal(d[k][ids],np.asarray(raw[k],float)) for k in poly)
        np.savez_compressed(OUT/f'gr/source-{n}.npz',**d)
        checks.append(dict(clock=n,original_endpoint_values_exact=True,pre_floor_to_saved_endpoint_relative=jumps))
        files += [INPUT/f'gr/source-{n}.npz']
    # Each independent field owns its initializer's generated files as well as
    # its output, so concurrent initializers cannot overwrite another worker.
    for action in FIELDS:
        worker=OUT/action
        for src in list((OUT/'sweep-0').rglob('*.npz'))+list((OUT/'gr').glob('source-*.npz'))+[OUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
            dst=worker/src.relative_to(OUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst)
        for part in ['sweep-1/photons','sweep-1/material']:(worker/part).mkdir(parents=True)
    files += [INPUT/n for n in ['plan.json','result.json','sources.json','driver-polynomial-audit.json','fields-receipt.json']]
    files += [INPUT/f'gr/fields-{n}-g{q}{ext}' for n,q in FIELDS.values() for ext in ['.npz','.json']]
    files += [prior.OUT/n for n in ['result.json','metric-result.json','symbolic.json','coarse-receipt.json','fine-receipt.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the completed224common15/16dense GR source to the SAME actual matter/photon stage equations, using the227actual-stage derivative repair, and judge the original paired2percent gate.',
        evidence='224common15/16source and field time gates pass; maximum field time3.46043e-5, quadrature7.35e-16, independentGR1.30e-15.227actual-stage return passed the shorterT/8global gate but early local differences remain. This extends that actual application, not another source-only diagnosis.',
        method='Reuse224dense source and original111/215stage histories. Evaluate fields at the union of the ORIGINAL stages and endpoints, preserving post-floor endpoint values and analytic same-polynomial derivatives. Reuse227same Radau angular emission, compensated high/low equations, exact conserved inverse and canonical model lifetimes. Capture actual accepted photons/ports. No independent diagnostic charge addition.',
        gates=read(prior.OUT/'plan.json')['decision'],budgets=CAPS,CPU_threads_per_process=1,virtual_GiB_per_process=16,
        parallel='Three independent GR fields in disjoint worker directories, CPUs4/6/8; then one metric and sequential original coarse/fine coupled returns. Existing232CPU0and235CPU2remain unchanged. Available memory58GB before launch; previous224peakRSS6.08GB.',
        forecast='224threefields1280.95/616.24/875.03seconds with original output times. Extra existing stage output times may cost1.5..3times; unmeasured estimate32..65minutes for parallel longest field.22715/29steps334/640seconds suggests41/79minutes for111/215steps if later marginal cost holds; later branches unmeasured. Generous3hours perfield,1hour metric,3/6hours coarse/fine; not completion guarantees.',
        scope='Same531cells,64/128clocks, physical3.219779ms and one GR feedback iterate. Last1/16, uniform prefix native accuracy, selfGR fixed point, full nonlinear evolution, EOS/derivative/spatial/boundary/observational/final-charge closure remain open.',
        stop='Original provenance, endpoint reproduction, metric, stage, physical moment, native anchor, branch, constitutive, conservation or paired-time failure; or wall/RSS cap. Preserve each canonical accepted checkpoint. Do not refine grid/clock or relax gates.',
        physical_steps=steps,output_times=len(clock),physical_horizon_seconds=float(clock[-1]),
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'source-endpoint-check.json',dict(classification='Counterexample candidate',passed=True,rows=checks))
    write(OUT/'symbolic.json',read(prior.OUT/'symbolic.json'))
    print(json.dumps(dict(steps=steps,output_times=len(clock),physical_horizon_seconds=float(clock[-1]))),flush=True)


def field(action):
    n,q=FIELDS[action];worker=OUT/action
    bind(prior.field,OUT=worker,INPUT=INPUT,Response=Response)(n,q)


def collect():
    for action,(n,q) in FIELDS.items():
        assert read(OUT/f'{action}-receipt.json')['error'] is None
        worker=OUT/action;assert read(worker/f'field-{n}-g{q}-check.json')['passed']
        for ext in ['.npz','.json']:os.link(worker/f'gr/fields-{n}-g{q}{ext}',OUT/f'gr/fields-{n}-g{q}{ext}')
    write(OUT/'field-collection.json',dict(classification='Counterexample candidate',passed=True,independent_workers=True))


def metric():
    assert read(OUT/'field-collection.json')['passed'];geometry.prior.saved=saved
    bind(prior.metric,OUT=OUT,Lapse=Lapse,representation=representation)()


def evolve(n):
    source=inspect.getsource(prior.evolve);assert source.count('for j in [1,2]:')==1
    source=source.replace('for j in [1,2]:','for j in range(1,16):')
    ns=dict(prior.evolve.__globals__,OUT=OUT,prior=model)
    exec(compile(source,__file__,'exec'),ns)
    (OUT/f'expanded-evolve-{n}.py').write_text(source);ns['evolve'](n)


def audit():bind(prior.audit,OUT=OUT)()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));prior.prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in FIELDS:field(action)
        elif action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
