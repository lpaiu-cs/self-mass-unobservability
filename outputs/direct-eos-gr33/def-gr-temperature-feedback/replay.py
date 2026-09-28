"""Independent saved-field replay of the declared temperature-feedback decision."""
from pathlib import Path
import hashlib
import json
import signal
import time
import numpy as np

OUT=Path('outputs/direct-eos-gr33/def-gr-temperature-feedback')
read=lambda p:json.loads(p.read_text());start=time.monotonic();signal.alarm(60)
result=read(OUT/'result.json');bindings=0
for name in ['plan.json','fast-plan.json']:
    for path,h in read(OUT/name)['bindings'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h,path
        bindings+=1
old=np.load(OUT.parent/'def-gr-interface-patch/p4-2048.npz');d=np.load(OUT/'response.npz')
assert np.array_equal(d['radius'],old['native_radius']) and np.array_equal(d['times'],np.arange(65)/64)
for key in ['weights','masks']:assert np.array_equal(d[key],old[key])
nr=len(d['radius']);w=d['weights'];masks=d['masks']
assert np.isfinite(d['correction']).all() and np.all(d['correction'][:,0]==0)
def norm(a):
    v=a[...,:nr];f=a[...,nr:]
    return np.stack([np.sqrt(np.sum(w*v*v,axis=-1)),np.sqrt(np.sum(w*f*f,axis=-1)),
        *[np.sqrt(np.sum(w[m]*v[...,m]**2,axis=-1)/w[m].sum()) for m in masks]],axis=-1)
base=np.concatenate([old['native_velocity'],old['native_scalar']],axis=1);scale=norm(base).max(0)
c=d['correction'];effects=norm(c[-1]).max(0)/scale
errors=norm(c[-1]-c[-2]).max(0)/scale;contour=norm(c[-1]-d['coarse_correction']).max(0)/scale
for j,(key,row) in enumerate(result['comparisons'].items()):
    for metric,value in [('correction_relative',effects[j]),('time_absolute_to_baseline',errors[j]),('contour_absolute_to_baseline',contour[j])]:
        assert abs(value-row[metric])<1e-13*max(abs(value),1e-100),(key,metric)
assert np.array_equal(d['total_native_velocity'],old['native_velocity']+c[-1,:,:nr])
assert np.array_equal(d['total_native_scalar'],old['native_scalar']+c[-1,:,nr:])
energy=d['correction_heat_energy'];assert energy[0]==energy[-1]==0
balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
assert balance<2e-13 and abs(balance-result['heat_balance'])<1e-15
endpoint=float(abs(d['base_q']-old['q']).max()/max(abs(old['q']).max(),1e-100))
assert endpoint<1e-10 and abs(endpoint-result['baseline_endpoint_replay'])<1e-15
assert read(OUT/'symbolic.json')['passed'] and read(OUT/'temperature-check.json')['passed']
assert read(OUT/'fast-pilot.json')['correction_native_max_relative']<2e-8
passed=bool(max(errors)<1e-6 and max(contour)<1e-6 and result['linear_residual']<1e-9 and result['feedback_residual']<1e-10)
assert passed==result['numerical_target_passed']
small=bool(passed and max(effects)<1e-4);assert small==result['temperature_only_feedback_below_target']
cost=result['total_compute_seconds']+read(OUT/'temperature-check.json')['seconds']+time.monotonic()-start
assert cost<720
row=dict(classification='Counterexample candidate',saved_replay_passed=True,numerical_target_passed=passed,
    decision_below_declared_target=small,source_bindings=bindings,native_material_samples=nr,
    largest_response_correction_relative=float(max(effects)),largest_time_absolute_to_baseline=float(max(errors)),
    largest_contour_absolute_to_baseline=float(max(contour)),heat_balance=balance,
    accounted_compute_seconds=cost,seconds=time.monotonic()-start,new_EOS_calls=0,new_time_paths=0,
    scope='Replay of the fixed p4 temperature-only loop and its declared absolute decision tolerance; no new spatial or full nonlinear physical certification.')
(OUT/'replay.json').write_text(json.dumps(row,indent=2)+'\n');signal.alarm(0);print(json.dumps(row))
