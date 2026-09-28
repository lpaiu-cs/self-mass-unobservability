"""Read the actual compensated return with its own applied GR geometry.

Counterexample candidate. The lower component belongs to the same stage solve;
it is retained separately from its high anchor to avoid rounding away feedback.
This endpoint source is not a dense-time or final-charge certificate.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json,os,resource,sys,time
import numpy as np
import read_complete_radau_history as prior
import apply_actual_stage_gr as returned

CHECK=len(sys.argv)>1 and sys.argv[1].startswith('check_')
ROOT=Path('native-returned-source245-work');OUT=ROOT/'check' if CHECK else ROOT/'full'
INPUT=returned.OUT if CHECK else Path('native-complete-return236-work')
read,write,sha,LD,AMP=prior.read,prior.write,prior.sha,prior.LD,prior.AMP
CAPS=dict(prepare=600,source=2700)
saved=lambda n:INPUT/f'sweep-1/photons/return-{n}.npz'


def bind(fn,**values):
    return FunctionType(fn.__code__,dict(fn.__globals__,**values),argdefs=fn.__defaults__)


def prepare():
    result=read(INPUT/'result.json')
    assert result['passed'] and result['actual_return_time_evolved']
    assert read(INPUT/'controller-status.json')['state']=='completed'
    assert read(INPUT/'metric-result.json')['passed']
    if not CHECK:assert read(ROOT/'check/result.json')['representation_controls_passed']
    assert not OUT.exists();OUT.mkdir(parents=True);files=[]
    for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr','metric']:(OUT/part).mkdir(parents=True)
    for src in list((INPUT/'sweep-0').rglob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','metric/metric-128-g8.npz']]:
        dst=OUT/src.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    horizons=[];steps=[]
    for n in [64,128]:
        p=np.load(saved(n));q=np.load(INPUT/f'recovered-{n}.npz')
        assert np.array_equal(p['joint_stage_times'],q['times']) and np.array_equal(p['joint_stage_weights'],q['weights'])
        assert np.array_equal(p['joint_collision_rates_scaled'],q['collision_rates'])
        assert read(INPUT/f'run-{n}.json')['passed']
        horizons.append(float(p['actual_step_edges'][-1]));steps.append(len(p['actual_step_edges'])-1)
        files += [saved(n),INPUT/f'recovered-{n}.npz',INPUT/f'run-{n}.json']
    assert abs(horizons[0]-horizons[1])<1e-18
    assert steps==([15,29] if CHECK else [111,215])
    files += [INPUT/n for n in ['plan.json','result.json','metric-result.json','coarse-receipt.json','fine-receipt.json','audit-receipt.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Read all accepted endpoint sources of the actual paired GR-return solve, with the identical returned metric, photon moments, material floor and boundary history.',
        method='Reuse224stable physical source map and exact conserved inverse. Replace only the readout driver by227StageDriver on the very same saved fine GR metric consumed by both actual return clocks. Keep high/low components separate; no zero-state nonlinear evolution or post-hoc unrelated charge.',
        decision='Original pressure0.2percent, mapping1e-12, material-ledger1e-8 and source-time2percent gates. Passing endpoints enables a later dense same-solution retarded charge, not final closure. The primary degree8 incident-mode geometry representation cannot represent this returned metric.',
        reuse='Use completed227for a real short-history check, then completed236full common15/16. No accepted physical stage, photon equation, GR field or EOS background is recomputed.',
        original_return_horizon_seconds=horizons[0],actual_steps=steps,budgets=CAPS,CPU_threads=1,CPU_affinity=6,virtual_GiB=16,
        forecast='224endpoint map measured185.92seconds for326steps. Short44step check roughly1..3minutes and full326roughly4..8minutes if marginal cost holds; tiny-component cancellation is unmeasured. Allow45minutes perreadout to avoid artificial interruption.',
        stop='Provenance, actual applied geometry, exact source times, pressure/mapping/ledger failure stops. Record any source-time failure without admitting a later physical feedback. No automatic new grid, physical path or tolerance relaxation.',
        scope='One actual compensated GR feedback iterate; endpoint source only. Continuous geometry, dense source, selfGR fixedpoint, EOS/derivative/spatial/boundary/nonlinear/observational/infinity gates remain.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as s
    a,x,d,b,y,e=s.symbols('a x d b y e');assert s.expand(a*(x+d)+b*(y+e)-(a*x+b*y)-(a*d+b*e))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Linear retained-source component identity only, not a nonlinear or physical error certificate.'))


def initialize():
    bind(prior.endpoint.initialize,OUT=OUT)()
    returned.geometry.OUT=OUT;Parent=prior.base.run.owner.Model
    class ReadModel(Parent):
        def __init__(self,n):
            super().__init__(n);primary=self.driver;driver=returned.StageDriver(primary)
            self.driver=self.redshift_driver=self.material.driver=driver;self.set_stage(0.)
            p=np.load(saved(n));checks=[]
            for t in p['actual_step_edges']:
                j=int(np.argmin(abs(driver.clock-t)));assert abs(driver.clock[j]-t)<1e-18
                values=driver.at(float(t));field,rates=self.geometry(float(t))
                assert np.array_equal(values['delta_u'],driver.g['delta_u'][j])
                assert np.array_equal(values['delta_lambda'],driver.g['delta_lambda'][j])
                assert np.array_equal(values['delta_lambda_rate'],driver.g['actual_delta_lambda_rate'][j])
                assert np.array_equal(field[0],values['delta_u']*(1./AMP))
                assert np.array_equal(field[2],values['delta_lambda']*(1./AMP))
                checks.append(float(t))
            write(OUT/f'applied-geometry-{n}.json',dict(classification='Counterexample candidate',passed=True,
                actual_endpoint_times=checks,returned_u_lambda_and_actual_rate_exact=True,
                same_metric_sha256=sha(OUT/'metric/metric-128-g8.npz')))
    prior.base.run.owner.Model=ReadModel


def recovered(n):
    p=np.load(INPUT/f'recovered-{n}.npz');return len(p['times'])//2,p['photon_moments'],p['radial_ports']


def source():
    endpoint=SimpleNamespace(**dict(vars(prior.endpoint),initialize=initialize))
    bind(prior.endpoints,OUT=OUT,saved=saved,recovered=recovered,endpoint=endpoint)()
    r=read(OUT/'endpoint-sources.json');r.update(
        representation_controls_passed=True,actual_returned_metric_applied=True,
        source_time_passed=r['passed'],component='low component of the same actual compensated solve',
        input_directory=str(INPUT),dense_time_representation_admitted=False,
        new_material_steps=0,new_photon_solves=0,new_GR_fields=0,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1].removeprefix('check_');assert action in CAPS
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    os.sched_setaffinity(0,{6});resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    prior.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
