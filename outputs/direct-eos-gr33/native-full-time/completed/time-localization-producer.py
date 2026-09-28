"""Saved-only decomposition of the failed full-input time pair."""
import time
from pathlib import Path
started=time.monotonic()
import numpy as np
import compare_full_incident_fluid_time as run

run.prior.joint.previous.original.inf.incident.native.deadline(10)
assert run.read(run.OUT/'audit-receipt.json')['seconds']+10<15
paths=[]
for n in [64,128]:
    p=dict(np.load(run.OUT/f'sweep-1/photons/pilot-{n}.npz'))
    native=np.sum(p['joint_stage_weights'][:,None,None]*p['joint_native_rates_scaled'],axis=0,dtype=run.LD)
    collision=np.sum(p['joint_stage_weights'][:,None,None]*p['joint_collision_rates_scaled'],axis=0,dtype=run.LD)
    floor=p['material_floor_discard_scaled'];actual=native+collision-floor
    paths.append(dict(data=p,native=native,collision=collision,floor=floor,actual=actual))
a,b=paths;rows=[]
for j,name in enumerate(['Etilde','H','B','S']):
    diff=a['actual'][:,j]-b['actual'][:,j];norm=np.sum(abs(diff));idx=np.argsort(abs(diff))[::-1]
    parts={k:a[k][:,j]-b[k][:,j] for k in ['native','collision','floor']}
    rows.append(dict(channel=name,relative=float(norm/max(np.sum(abs(b['actual'][:,j])),run.LD('1e-290'))),
        component_L1_over_total={k:float(np.sum(abs(v))/norm) for k,v in parts.items()},
        cells=[dict(cell=int(i),radius_cm=float(b['data']['radius_E'][i]),share=float(abs(diff[i])/norm),
                    difference=float(diff[i]),**{k:float(v[i]) for k,v in parts.items()}) for i in idx[:12]]))
gas_mid=[]
for path in paths:
    p=path['data'];ids=np.array([int(np.argmin(abs(p['joint_stage_times']-t))) for t in a['data']['actual_step_edges'][1:]])
    assert np.max(abs(p['joint_stage_times'][ids]-a['data']['actual_step_edges'][1:]))<1e-18
    gas_mid.append(p['joint_stage_conserved_scaled'][ids])
closing=(np.sum(abs(gas_mid[0]-gas_mid[1]),axis=2)/np.maximum(np.sum(abs(gas_mid[1]),axis=2),run.LD('1e-290'))).astype(float).tolist()
ph=[]
for j,name in [(0,'photon_energy'),(1,'material_reference_energy'),(6,'material_pressure')]:
    diff=a['data']['moments'][-1,j]-b['data']['moments'][-1,j];norm=np.sum(abs(diff));idx=np.argsort(abs(diff))[::-1]
    ph.append(dict(channel=name,cells=[dict(cell=int(i),share=float(abs(diff[i])/norm),difference=float(diff[i])) for i in idx[:8]]))
result=dict(classification='Counterexample candidate',trajectory_replayed=False,rows=rows,photon_rows=ph,
    same_time_closing_material_relative=closing,closing_material_order=['B','S','Eref','H'],
    seconds=time.monotonic()-started,source_sha256=run.sha(__file__),
    scope='Exact saved-integral and cell localization. Closing stages are before the endpoint floor projection. Component L1 values may cancel and are not disjoint error shares. No unique temporal-error cause established.')
run.write(run.OUT/'time-localization.json',result)
print(result)
