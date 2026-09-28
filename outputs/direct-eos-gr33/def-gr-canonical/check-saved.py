"""Replay the stored acceptance decision; do not run another GR path."""
from pathlib import Path
import json
import hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
read=lambda p:json.loads(p.read_text())
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert read(OUT/'symbolic.json')['passed']
assert read(OUT/'regular/symbolic.json')['passed']
assert not read(OUT/'evolution/result.json')['time_passed']
assert read(OUT/'modes/result.json')['time_passed']
assert not read(OUT/'modes/result.json')['passed']
assert not read(OUT/'nested/result.json')['passed']
result=read(OUT/'regular/result.json');assert result['projection_passed'] and not result['passed']
fields=['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']
maximum_replay_error=0.;paths=0
for directory in ['evolution','modes','nested','regular']:
    for path in (OUT/directory).glob('*.json'):
        row=read(path)
        if 'history' not in row:continue
        state=dict(np.load(path.with_suffix('.npz')));v=state['native_velocity'];psi=state['native_scalar'];w=state['weights'];masks=state['masks']
        expected=np.column_stack([np.sqrt((v*v)@w),np.sqrt((psi*psi)@w),*[np.sqrt((v[:,m]**2)@w[m]/w[m].sum()) for m in masks]])
        saved=np.array([[h[f] for f in fields] for h in row['history']])
        err=float(np.max(abs(saved-expected)/np.maximum(abs(expected).max(axis=0),1e-100)))
        assert err<1e-12,(path,err);maximum_replay_error=max(maximum_replay_error,err);paths+=1
        assert np.all(v[0]==0) and np.all(psi[0]==0)
        energy=state['heat_energy'];balance=abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100)
        assert float(balance)<2e-13
for field in fields:
    a,b,c=[np.array([h[field] for h in read(OUT/f'regular/fine-{n}.json')['history']]) for n in [128,256,512]]
    norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
    measured=dict(previous=d1,last=d2,order=np.log2(d1/d2))
    for label in ['coefficient','outer','spatial']:
        other=np.array([h[field] for h in read(OUT/f'regular/{label}-512.json')['history']])
        measured[label]=np.max(abs(c-other))/norm
    for key,value in measured.items():assert abs(result['comparisons'][field][key]-value)<1e-12
assert all(v['last']<.02 and v['order']>1.5 and v['coefficient']<.02 and v['outer']<.002 for v in result['comparisons'].values())
assert all(result['comparisons'][f]['spatial']>.02 for f in [fields[0],fields[2],fields[3]])
assert result['comparisons'][fields[1]]['spatial']<.02
assert not result['original_failure_resolved'] and not result['full_dynamic_charge_solved']
reconstruction=read(OUT/'regular/reconstruction.json')
assert reconstruction['physical_pressure_metric_reconstruction_checked']
assert reconstruction['independent_baryon_pressure_backward_error']<1e-12
assert reconstruction['heat_faces_not_aligned_to_spatial_nodes']==0
assert read(OUT/'reconstruction-failure.json')['passed'] is False
bindings=0
for plan in OUT.rglob('*plan.json'):
    row=read(plan)
    for rel,sha in row.get('bindings',{}).items():
        p=Path(rel)
        if not p.is_absolute():p=ROOT/p
        assert digest(p)==sha,(plan,rel);bindings+=1
    for key,name in [('source','def_gr_canonical.py'),('canonical_source','def_gr_canonical.py')]:
        if key not in row:continue
        if plan.name=='assembly-plan.json' and key=='source':name='def_gr_energy_fem.py'
        if plan.name=='reconstruction-plan.json':
            p=OUT/('reconstruct-regular.py' if plan.parent.name=='regular' else 'reconstruct.py')
        else:p=ROOT/'verification'/name
        assert digest(p)==row[key],(plan,key);bindings+=1
check=dict(classification='Counterexample candidate',passed=True,stored_paths_replayed=paths,
    maximum_native_readout_replay_error=maximum_replay_error,source_bindings_checked=bindings,
    fixed_grid_propagation_passed=True,spatial_velocity_passed=False,original_failure_resolved=False,
    full_dynamic_charge_solved=False)
(OUT/'saved-check.json').write_text(json.dumps(check,indent=2)+'\n');print(json.dumps(check))
