"""Consume the measured finite collision remainder in the actual coupled return."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
from functools import lru_cache
import io,json,resource,signal,sys,time
import numpy as np
file_load=np.load


@lru_cache(maxsize=8)
def frozen_bytes(path):return Path(path).read_bytes()


def memory_reads():
    def load(path,*args,**kwargs):
        if isinstance(path,(str,Path)) and str(path).endswith('.npz'):
            p=Path(path)
            # Cache immutable INPUT bytes, never arrays or this phase's outputs.
            data=p.read_bytes() if 'finite-collision' in p.parts or 'checkpoint' in p.name else frozen_bytes(str(p))
            path=io.BytesIO(data)
        return file_load(path,*args,**kwargs)
    np.load=load


memory_reads()
import repair_native_finite_collision_remainder as repair
run=repair.run;SOURCE=repair.OUT;OUT=SOURCE/'applied';TOTAL=OUT/'total-photons';RETURN=OUT/'return';GR=RETURN/'gr'
write=run.write;sha=run.sha;AMP=run.AMP;configure=run.configure


def prepare():
    assert not OUT.exists()
    for p in [OUT,TOTAL,RETURN,GR]:p.mkdir()
    pilot=json.loads((SOURCE/'completed-pilot.json').read_text());assert all(r['passed'] for r in pilot['rows']) and not pilot['eligible']
    assert pilot['forecast_upper_seconds']<330
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Propagate the actual finite radiation collision remainder through the accepted simultaneous photon/E/H operator, add that increment to the saved response, then apply the new TOTAL photon transfers to finite material and actual GR.',
        decision='Quantify the effect of the previously omitted finite collision term on the positive conditional charge. Preserve its sampled/arithmetic limits; this is a deferred correction, not full nonlinear or uniform microscopic certification.',
        resource_reassessment='The original245s readout allowance was insufficient: after81.523589098s the measured conservative projection312.241693236s exceeds163.476410902s remaining. Preserve that ineligibility. Authorize one330s remaining-state readout cap, adding166.523589098s to the earlier cap; unchanged17 states, all accepted points reused. Longdouble is necessary because the binary64 remainder was unresolved, not because a physical acceptance gate was relaxed.',
        cheaper_alternatives='No background or native roots are replayed. Reuse zero/midpoint/endpoint and all old photons. Linear superposition solves only the new photon defect; full saved photons are not re-integrated. The subsequent finite material equation must consume the new total force and is not interchangeable with an old path.',
        budgets=dict(sample_seconds=330,photon_pilot_seconds=90,photon_production_seconds=900,material_pilot_seconds=60,material_production_seconds=780,source_seconds=90,GR_seconds=90,CPU_threads=1,virtual_GiB=3),
        forecast='Remainder14 remaining knots at twice the measured slowest point plus10s. Photon4/8 prefixes use the completed Phase142 step costs as late-cost floors; require2x summed forecast below900s. Material2-step prefixes use completed full raw-call counts; require2x summed forecast below780s.',
        gates=dict(number=1e-10,linear=1e-12,balance=1e-8,time=.02,finite_resolution=.002,pressure=.002,quadrature=.002,independent_GR=1e-9,angular=1e-12),
        limitations='Actual radiation primitive/rate remainder on17 saved states of the interpolation model. Native cross sections remain binary64. The16eps arithmetic floor, half-amplitude comparison and finite waveform changes are not rigorous error bounds. Native derivative, further coupling, floor/exterior, nonlinear GR and observation remain open.',
        stop='One source readout and one two-path correction/return only. Stop on a gate or cap; no automatic further waveform loop, time/grid/horizon increase or relaxed thresholds.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(repair.__file__),Path(repair.old.__file__),Path(run.__file__),SOURCE/'completed-pilot.json',run.OUT/'audit.json',run.GR/'result.json']}))


def sample():
    assert not (OUT/'sample-failure.json').exists()
    assert not (SOURCE/'production.json').exists();start=time.monotonic();signal.alarm(330);s=repair.State();rows=[]
    # Zero increments are identically zero at either precision.
    for suffix in ['npz','json']:(SOURCE/f'point-0.{suffix}').write_bytes((repair.old.OUT/f'point-0.{suffix}').read_bytes())
    for k in range(17):
        row=json.loads((SOURCE/f'point-{k}.json').read_text()) if (SOURCE/f'point-{k}.json').exists() else s.point(k)
        assert row['passed'],row;rows.append(row)
    p=dict(classification='Counterexample candidate',passed=True,rows=rows,new_states=14,saved_points_reused=[0,8,16],seconds=time.monotonic()-start,
        maximum_remainder_over_linear=max(r['rows'][0]['remainder_over_linear'] for r in rows),uniform_remainder_bound=False)
    write(SOURCE/'production.json',p);signal.alarm(0);print(json.dumps({k:v for k,v in p.items() if k!='rows'}),flush=True)


def resume_sample():
    assert not (OUT/'IO-reallocation-plan.json').exists()
    assert not (SOURCE/'production.json').exists();failure=json.loads((OUT/'sample-failure.json').read_text())
    assert 'TimeoutError' in failure['error']
    saved=[k for k in range(17) if (SOURCE/f'point-{k}.json').exists()];missing=[k for k in range(17) if k not in saved]
    write(OUT/'IO-reallocation-plan.json',dict(classification='Counterexample candidate',
        failure=failure,observed='During the interrupted readout the process spent60 CPU seconds in279 wall seconds and was observed waiting in p9_client_rpc. Saved-state initialization dominated before the first new collision point.',
        repair='Load each immutable NPZ byte stream sequentially into BytesIO before NumPy reads its members. No array, physical owner, acceptance gate, state count or interval changes.',
        budgets=dict(resume_sample_seconds=330,material_production_seconds=450,previous_material_production_seconds=780,aggregate_budget_increase_seconds=0),
        reuse=saved,remaining=missing,
        forecast='One missing-state prefix after initialization; require twice its measured cost times the remaining count plus10 seconds below the330s remaining allowance. Stop if ineligible. Material pilot must now fit450s.',
        stop='One IO repair and resume only; no automatic second resume or budget extension.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'first-producer.py',OUT/'plan.json',OUT/'sample-failure.json']}))
    start=time.monotonic();signal.alarm(330);s=repair.State();row=s.point(missing[0]);assert row['passed']
    forecast=2*row['seconds']*(len(missing)-1)+10;remaining=330-(time.monotonic()-start)
    pilot=dict(seconds=time.monotonic()-start,point=row,forecast_seconds=forecast,remaining_seconds=remaining,eligible=forecast<remaining)
    write(OUT/'IO-resume-pilot.json',pilot);print(json.dumps(pilot),flush=True);assert pilot['eligible']
    for k in missing[1:]:assert s.point(k)['passed']
    rows=[json.loads((SOURCE/f'point-{k}.json').read_text()) for k in range(17)];assert all(r['passed'] for r in rows)
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,new_states=len(missing),saved_points_reused=saved,
        seconds=time.monotonic()-start,interrupted_readout_seconds=failure['seconds'],
        maximum_remainder_over_linear=max(r['rows'][0]['remainder_over_linear'] for r in rows),uniform_remainder_bound=False)
    write(SOURCE/'production.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


def cached_sample():
    assert not (OUT/'cached-input-plan.json').exists()
    assert not (SOURCE/'production.json').exists();prior=json.loads((OUT/'IO-resume-pilot.json').read_text());assert not prior['eligible']
    prior_spent=json.loads((OUT/'resume_sample-failure.json').read_text())['seconds']
    budget=int(330-prior_spent);saved=[k for k in range(17) if (SOURCE/f'point-{k}.json').exists()]
    missing=[k for k in range(17) if k not in saved]
    write(OUT/'cached-input-plan.json',dict(classification='Counterexample candidate',
        reassessment='The previous resume stays ineligible. Code inspection found repeated full reads of the same95MB background and60MB historical photon archives in inherited constructors. BytesIO per open did not remove those repeated reads.',
        repair='Cache at most8 immutable input byte streams before importing owners; still decode independent arrays. Exclude all phase outputs and checkpoints. This is an input-IO change, not a second unchanged resume or a physical gate change.',
        budget_seconds=budget,aggregate_budget_increase_seconds=0,material_production_cap_seconds=450,prior_spent_seconds=prior_spent,
        reuse=saved,remaining=missing,forecast='Use one missing point and require2x its measured cost times the uncomputed count plus10s within the unused original330s resume allowance.',
        stop='If this measured cache repair does not fit, close the source attempt without a response or another dispatch.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'second-producer.py',OUT/'IO-resume-pilot.json',OUT/'resume_sample-failure.json']}))
    start=time.monotonic();signal.alarm(budget);s=repair.State();row=s.point(missing[0]);assert row['passed']
    forecast=2*row['seconds']*(len(missing)-1)+10;remaining=budget-(time.monotonic()-start)
    pilot=dict(seconds=time.monotonic()-start,point=row,forecast_seconds=forecast,remaining_seconds=remaining,eligible=forecast<remaining,
        input_cache=frozen_bytes.cache_info()._asdict())
    write(OUT/'cached-input-pilot.json',pilot);print(json.dumps(pilot),flush=True);assert pilot['eligible']
    for k in missing[1:]:assert s.point(k)['passed']
    rows=[json.loads((SOURCE/f'point-{k}.json').read_text()) for k in range(17)];assert all(r['passed'] for r in rows)
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,new_states=len(missing),saved_points_reused=saved,
        seconds=time.monotonic()-start,earlier_readout_seconds=json.loads((OUT/'sample-failure.json').read_text())['seconds']+prior['seconds'],
        maximum_remainder_over_linear=max(r['rows'][0]['remainder_over_linear'] for r in rows),uniform_remainder_bound=False,
        input_cache=frozen_bytes.cache_info()._asdict())
    write(SOURCE/'production.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


class Response(run.Response):
    def __init__(self,reference=128,steps=128):
        super().__init__(reference,steps)
        for key in ['motion','energy_offset','mechanical','xi']:setattr(self,key,np.zeros_like(getattr(self,key)))
    def defect(self,k):
        if k not in self.defects:
            d=dict(np.load(SOURCE/f'point-{k}.npz'));self.defects[k]={key:np.asarray(d[key],float) for key in ['photon','bound','escape']}
            for j in list(self.defects):
                if j not in [k-1,k,k+1]:del self.defects[j]
        return self.defects[k]
    run=FunctionType(run.Response.run.__code__,dict(run.Response.run.__globals__,OUT=OUT),argdefs=run.Response.run.__defaults__)


photon_worker=FunctionType(run.photon_worker.__code__,dict(run.photon_worker.__globals__,OUT=OUT,Response=Response),argdefs=run.photon_worker.__defaults__)


def photon_pilot():
    assert json.loads((SOURCE/'production.json').read_text())['passed'];assert not (OUT/'pilot.json').exists()
    configure();start=time.monotonic();signal.alarm(90)
    rows=[photon_worker(n,f'pilot-{n}',k) for n,k in [(64,4),(128,8)]]
    prior=json.loads((run.OUT/'result.json').read_text())['paths'];estimates=[]
    for r in rows:
        p=next(x for x in prior if x['steps']==r['steps'])
        cost=max(r['stepping_seconds']/r['new_steps'],p['stepping_seconds']/p['new_steps'])
        estimates.append(17*r['operator_point_seconds']/r['operator_points']+(r['steps']-r['completed_steps'])*cost+20)
    a=np.load(OUT/'pilot-64.npz')['moments'][-1,[0,1,2,3,5,6]];b=np.load(OUT/'pilot-128.npz')['moments'][-1,[0,1,2,3,5,6]]
    errors=(np.sum(abs(a-b),axis=1)/np.maximum(np.sum(abs(b),axis=1),1.)).tolist()
    p=dict(classification='Counterexample candidate',rows=rows,equal_horizon=errors,forecast_each_seconds=estimates,upper_seconds=2*sum(estimates),eligible=max(errors)<.02 and 2*sum(estimates)<900,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',p);signal.alarm(0);print(json.dumps(p),flush=True)


photon_production=FunctionType(run.photon_production.__code__,dict(run.photon_production.__globals__,OUT=OUT,photon_worker=photon_worker))


def photon_path(n,r):return TOTAL/f'steps-{n}-reference-{r}.npz'


def combine():
    assert json.loads((OUT/'result.json').read_text())['passed'];configure();m=Response();rows=[]
    keys=['moments','delta_packet_scaled_occupation','delta_material','ledger','escape','collision_transfer','radial_ports',
          'photon_history_scaled_occupation','material_history','accepted_angular_luminosity']
    for n in [64,128]:
        original=dict(np.load(run.photon_path(n,128)));correction=dict(np.load(OUT/f'steps-{n}-reference-128.npz'))
        assert np.array_equal(original['t'],correction['t']) and np.array_equal(original['accepted_angular_times'],correction['accepted_angular_times'])
        total=dict(original)
        for key in keys:total[key]=original[key]+correction[key]
        total['moments'][:,4]=np.sum(abs(total['photon_history_scaled_occupation'])*(m.Eweight/m.scale),axis=(2,3))
        np.savez_compressed(photon_path(n,128),**total)
        relative=(np.max(np.sum(abs(correction['moments'][:,[0,1,2,3,5,6]]),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(original['moments'][:,[0,1,2,3,5,6]]),axis=2),axis=0),1.)).tolist()
        rows.append(dict(steps=n,correction_over_original=relative))
    write(TOTAL/'result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,original_photon_paths_replayed=False,
        interpretation='Exact linear superposition for the same stored tangent operator and frozen prescribed mechanical input. New finite collision remainder is included once; its new-state Jacobian and remainder are not certified.'))


class Material(run.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference,steps);p=dict(np.load(photon_path(steps,reference)));ids=[int(np.argmin(abs(p['t']-t))) for t in self.t]
        assert np.max(abs(p['t'][ids]-self.t))<1e-18;c=p['collision_transfer'][ids]
        self.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/self.a,c[:,:,0],c[:,:,1]],axis=1)/AMP
    run=FunctionType(run.Material.run.__code__,dict(run.Material.run.__globals__,OUT=RETURN),argdefs=run.Material.run.__defaults__)


material_worker=FunctionType(run.material_worker.__code__,dict(run.material_worker.__globals__,OUT=RETURN,Material=Material),argdefs=run.material_worker.__defaults__)
material_pilot_owner=FunctionType(run.material_pilot.__code__,dict(run.material_pilot.__globals__,OUT=TOTAL,RETURN=RETURN,material_worker=material_worker,matter=SimpleNamespace(OUT=run.RETURN)))
def material_pilot():
    material_pilot_owner();p=json.loads((RETURN/'pilot.json').read_text())
    p.update(eligible=p['upper_seconds']<450,production_cap_seconds=450)
    write(RETURN/'pilot.json',p);print(json.dumps(p),flush=True)
material_production=FunctionType(run.material_production.__code__,dict(run.material_production.__globals__,OUT=RETURN,worker=material_worker,
    signal=SimpleNamespace(alarm=lambda seconds:signal.alarm(min(seconds,450)))))
sources=FunctionType(run.sources.__code__,dict(run.sources.__globals__,OUT=RETURN,GR=GR,Material=Material,old=SimpleNamespace(previous=run.matter.old.old.previous,photon_path=photon_path)))


def charge():
    fn=FunctionType(run.charge.__code__,dict(run.charge.__globals__,RETURN=RETURN,GR=GR));fn()
    r=json.loads((GR/'result.json').read_text());before=json.loads((run.GR/'result.json').read_text())
    r.update(finite_collision_remainder_applied=True,previous_native_finite_motion_endpoint=before['endpoint_collision_charge'],
             endpoint_finite_collision_change=r['endpoint_collision_charge']-before['endpoint_collision_charge'],
             full_nonlinear_radiation=False,uniform_nonlinear_remainder_bound=False)
    write(GR/'result.json',r);print(json.dumps(r),flush=True)


if __name__=='__main__':
    cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));signal.signal(signal.SIGALRM,run.native.forcing.history.flow.old.optical.timeout)
    memory_reads()
    action=sys.argv[1];start=time.monotonic()
    try:globals()[action]()
    except Exception as exc:
        failure=OUT/f'{action}-failure.json'
        if not failure.exists():write(failure,dict(error=repr(exc),seconds=time.monotonic()-start))
        raise
