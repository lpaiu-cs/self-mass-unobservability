"""Return the actual native collision response to free matter and GR."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,signal,sys,time
import numpy as np
import repair_native_collision_defect as repaired
import verify_native_corrected_motion_feedback as previous

photons=repaired.response;physical=previous.physical;AMP=physical.AMP
OUT=photons.OUT/'material-charge';GR=OUT/'gr';write=repaired.write;sha=repaired.sha


def configure():physical.configure();repaired.configure()
def photon_path(n,r):return photons.OUT/f'steps-{n}-reference-{r}.npz'


class Material(previous.previous.Material):
    def __init__(self,reference,steps=128):
        assert reference==128;super().__init__(reference,steps)
        for key in self.metric:
            if key.startswith('delta_'):self.metric[key]=np.zeros_like(self.metric[key])
        p=dict(np.load(photon_path(steps,reference)));ids=[int(np.argmin(abs(p['t']-t))) for t in self.t]
        assert np.max(abs(p['t'][ids]-self.t))<1e-18;c=p['collision_transfer'][ids]
        self.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/self.a,c[:,:,0],c[:,:,1]],axis=1)/AMP
    run=FunctionType(previous.Material.run.__code__,dict(previous.Material.run.__globals__,OUT=OUT),argdefs=previous.Material.run.__defaults__)


worker=FunctionType(physical.worker.__code__,dict(physical.worker.__globals__,OUT=OUT,Material=Material,configure=configure))


def prepare():
    assert not OUT.exists() and json.loads((photons.OUT/'result.json').read_text())['passed'];OUT.mkdir();GR.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply native-minus-table collision photon/H/energy/momentum histories to actual free material and the existing retarded GR charge; quantify whether the earlier charge survives.',
        reuse='Two64/128 paths on the existing fine background,531 cells,3.434ms,corrected SSP stage convention, original conservative primitive pressure and retarded one-observer readout. No new EOS roots, mesh, background evolution or horizon.',
        model='Isolated constitutive correction with zero new metric forcing. Preserve its actual amplitude; no rescaling to pass donor/small-state gates. Free-material motion returning to photons and total coupling amplification remain open.',
        budgets=dict(pilot_seconds=60,production_seconds=900,source_seconds=90,GR_seconds=90,CPU_threads=1,virtual_GiB=3),
        forecast='Measure two actual prefix steps per path; extrapolate using prior completed late-CFL raw-call counts and previous full worker wall time as a floor. Require twice sum of forecasts below900s. Resume the prefixes; no repeat of completed physics.',
        gates=dict(owner=1e-8,conservation=1e-8,directional=.002,branch=.01,small_state=1e-6,time=.02,pressure=.002,quadrature=.002,independent_GR=1e-9),
        stop='Stop at a physical regime or time-comparison failure. No automatic finer paths, amplitude reduction, relaxed gates or extra response iterations.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(previous.__file__),Path(physical.__file__),Path(repaired.__file__),photons.OUT/'result.json',repaired.OUT/'result.json',previous.OUT/'production.json']}))


def pilot():
    assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.alarm(60);rows=[]
    for n in [64,128]:
        worker(n,128,f'pilot-{n}',2,None);rows.append(json.loads((OUT/f'pilot-{n}.json').read_text()))
    old=json.loads((previous.OUT/'production.json').read_text())['paths'];estimates=[]
    for r in rows:
        p=next(v for v in old if (v['steps'],v['reference'])==(r['steps'],128))
        estimates.append(max(p['raw_owner_calls']*r['seconds']/r['raw_owner_calls']+r['worker_wall_seconds']-r['seconds']+15,p['worker_wall_seconds']))
    result=dict(classification='Counterexample candidate',rows=rows,forecast_each_seconds=estimates,
        upper_seconds=2*sum(estimates),eligible=all(r['passed'] for r in rows) and 2*sum(estimates)<900,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def production():
    assert json.loads((OUT/'pilot.json').read_text())['eligible'];assert not (OUT/'production.json').exists()
    start=time.monotonic();signal.alarm(900);rows=[]
    for n in [64,128]:
        worker(n,128,f'steps-{n}-reference-128',None,f'pilot-{n}')
        rows.append(json.loads((OUT/f'steps-{n}-reference-128.json').read_text()))
    write(OUT/'production.json',dict(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start));signal.alarm(0)


def sources():
    assert json.loads((OUT/'production.json').read_text())['passed'];configure()
    source=previous.previous.old.source_owner.source
    source=source.replace('signal.alarm(75)','signal.alarm(90)')
    source=source.replace('[[64,128],[128,128],[128,64]]','[[64,128],[128,128]]')
    source=source.replace('background=compare(histories[2],histories[1]),','')
    source=source.replace(',stress_background=compare(allstress[2],allstress[1])','')
    scope=dict(previous.source_scope,OUT=OUT,GR=GR,Material=Material,
               photons=SimpleNamespace(path=photon_path),material_path=lambda n,r:OUT/f'steps-{n}-reference-{r}.npz')
    exec(compile(source,__file__,'exec'),scope);(OUT/'expanded-source.py').write_text(source);scope['sources']()
    assert json.loads((OUT/'sources.json').read_text())['passed']


def charge():
    import verify_native_stage_energy_charge as readout
    assert json.loads((OUT/'sources.json').read_text())['passed'];configure();start=time.monotonic();signal.alarm(90)
    m=previous.GRResponse();paths={}
    for n,q in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(GR/f'source-{n}-reference-128.npz'));paths[n,q]=readout.read(m,d,q)
        np.savez_compressed(GR/f'wave-{n}-g{q}.npz',**paths[n,q])
    fine=paths[128,8]['free_scalar'];norm=max(np.max(abs(fine)),1e-300)
    terr=float(np.max(abs(fine-paths[64,8]['free_scalar']))/norm);qerr=float(np.max(abs(fine-paths[128,4]['free_scalar']))/norm)
    d=dict(np.load(GR/'source-128-reference-128.npz'));direct,coordinate=readout.independent.direct(m,d,8)
    ierr=abs(direct/paths[128,8]['direct_scalar'][-1]-1)
    prior=json.loads((repaired.forcing.history.OUT/'result.json').read_text())
    result=dict(classification='Counterexample candidate',passed=terr<.02 and qerr<.002 and ierr<1e-9,
        endpoint_collision_charge=float(fine[-1]),time_relative=terr,quadrature_relative=qerr,independent_direct_relative=ierr,
        retarded_coordinate_residual=coordinate,prior_native_source_result=prior,seconds=time.monotonic()-start,
        actual_native_collision_in_photon_thermal_H_evolution=True,actual_collision_transfer_in_free_material=True,
        actual_collision_response_in_GR=True,motion_returned_to_photons=False,uniform_derivative_bound=False,
        coupled_fixed_point_verified=False,final_charge_solved=False,full_goal_complete=False)
    write(GR/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':
    signal.signal(signal.SIGALRM,repaired.forcing.history.flow.old.optical.timeout);action=sys.argv[1];start=time.monotonic()
    try:globals()[action]()
    except Exception as exc:write(OUT/f'{action}-failure.json',dict(error=repr(exc),seconds=time.monotonic()-start));raise
