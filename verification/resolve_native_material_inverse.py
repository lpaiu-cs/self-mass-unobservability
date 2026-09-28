"""Resolve native atmospheric primitive accuracy before differencing it.

Counterexample candidate. Same physical input, probe16, fluxes and gates.
"""
from pathlib import Path
from types import FunctionType,MethodType
import json,resource,sys,time
import numpy as np
import return_native_pressure_matter as run

OUT=run.OUT;read,write,sha=run.read,run.write,run.sha
BASE_INITIALIZE=run.initialize;BASE_PATHS=run.paths
CAPS=dict(check=30,retry=20,pilot=70,production=450,residual=40)


def paths(sweep):
    p,m=BASE_PATHS(sweep)
    return p,m.with_name('material-precise-inverse') if sweep else m


def precise_inverse(m):
    f=m.model.flow;fn=f.primitive.__func__;constants=fn.__code__.co_consts
    assert constants.count(2e-14)==1 and constants.count(1e-12)==1
    assert not m.cache,'Precision must be installed before reference fluxes are cached'
    code=fn.__code__.replace(co_consts=tuple(1e-14 if v in (2e-14,1e-12) else v for v in constants))
    f.primitive=MethodType(FunctionType(code,fn.__globals__,argdefs=fn.__defaults__,closure=fn.__closure__),f)


def initialize():
    run.original.paths=paths;run.paths=paths;BASE_INITIALIZE();Parent=run.c.Material
    class Material(Parent):
        def __init__(self,reference=128,steps=128):
            super().__init__(reference,steps);precise_inverse(self)
    run.c.Material=Material


def check(retry=False):
    failed=BASE_PATHS(1)[1]/'pilot-64.npz';assert not read(failed.with_suffix('.json'))['passed']
    if retry:
        assert read(OUT/'primitive-check-receipt.json')['error']=='AssertionError()'
        assert read(OUT/'primitive-check-receipt.json')['seconds']+CAPS['retry']<30
        write(OUT/'primitive-dispatch-repair.json',dict(classification='Counterexample candidate',
            correction='Live primitive code already uses2e-14 energy tolerance from the earlier tight_primitive owner. The first guard inspected a historical2e-11 constant and failed before any numerical comparison. Preserve that failed plan and producer; correct the live constant lookup. The originally proposed1e-14 target and existing probe range are unchanged.',
            source_sha256=sha(__file__),preserved_source_sha256=sha(OUT/'primitive-introspection-producer.py'),budget_seconds=CAPS['retry']))
    else:
        assert not (OUT/'primitive-precision-plan.json').exists()
        prepare_check(failed)
    initialize();m=run.c.Material(128,64);d=np.load(failed);values=[];rows=[]
    for factor in [.5,1.,2.]:
        m.min_probe=np.inf;r,_,dt=m.rhs(float(d['time']),d['delta_scaled'],factor);values.append(r)
        rows.append(dict(factor=factor,epsilon=float(m.min_probe),cfl=float(dt)))
    errors=[(np.sum(abs(v-values[1]),axis=1)/np.maximum(np.sum(abs(values[1]),axis=1),1.)).astype(float).tolist() for v in [values[0],values[2]]]
    np.savez_compressed(OUT/'primitive-precision-rates.npz',t=d['time'],state=d['delta_scaled'],rates=values)
    result=dict(classification='Counterexample candidate',passed=bool(np.max(errors)<.002),half_nominal_double=errors,
        rows=rows,probe_unchanged=16.,inverse_tolerance=1e-14,actual_trajectory_accepted=False,final_charge_conclusion='unadjudicated')
    write(OUT/'primitive-precision-check.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def prepare_check(failed):
    files=[Path(__file__),Path(run.__file__),failed,failed.with_suffix('.json'),OUT/'probe-localization.json',OUT/'pilot-receipt.json']
    write(OUT/'primitive-precision-plan.json',dict(classification='Counterexample candidate',
        failure='The actual material pilot fails the original0.2percent derivative gate. Repeating the existing8/16/32arithmetic probes localizes the largest B discrepancy to atmospheric cells146/147, with zero discrepancy on analytic deep faces.',
        hypothesis='The atmospheric native primitive inverse stops at2e-11 relative energy residual and1e-12 pressure-settling residual. Differencing small local changes can amplify this inverse error; the existing global probe is constrained by other components/cells.',
        change='Tighten both inverse stopping constants to1e-14 BEFORE constructing cached reference fluxes. Keep native EOS, reconstruction, HLL/donor branches, SSP2, physical eta, probe16and original acceptance gates. One fixed precision trial, no probe or precision ladder.',
        admission='All half/nominal/double RHS differences on the saved actual endpoint must be below0.2percent before rerunning the same4/8macro-step material prefix in a separate folder.',
        budgets=CAPS,total_original_action_budget_seconds=640,
        stop='Failure stops this route. No automatic smaller step, larger probe, new EOS or full material run.',
        final_charge_conclusion='unadjudicated',bindings={str(p):sha(p) for p in files}))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'primitive-{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        spent=sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+read(OUT/'probe-localization.json')['seconds']
        assert spent+CAPS[action]<=640,(spent,CAPS[action])
        if action in ['check','retry']:check(action=='retry')
        else:
            assert read(OUT/'primitive-precision-check.json')['passed']
            for p,h in read(OUT/'primitive-precision-plan.json')['bindings'].items():
                frozen=OUT/'primitive-introspection-producer.py' if Path(p)==Path(__file__) else p
                assert sha(frozen)==h,p
            assert sha(__file__)==read(OUT/'primitive-dispatch-repair.json')['source_sha256']
            paths(1)[1].mkdir(exist_ok=True);run.initialize=initialize;run.paths=paths
            if action in ['pilot','production']:run.material(action=='pilot')
            else:run.residual()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
