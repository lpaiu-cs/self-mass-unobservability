"""Replay every saved native history and its unchanged decision, without evolution."""
from pathlib import Path
import hashlib
import json
import numpy as np

OUT=Path(__file__).resolve().parent
FIELDS=['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']
write=lambda p,d:p.write_text(json.dumps(d,indent=2)+'\n')
cases={};bindings=0;max_error=0.
for p in sorted(OUT.rglob('*.json')):
    data=json.loads(p.read_text())
    for source,sha in data.get('bindings',{}).items():
        assert hashlib.sha256(Path(source).read_bytes()).hexdigest()==sha,(p,source)
        bindings+=1
    if 'history' not in data:continue
    arrays=np.load(p.with_suffix('.npz'));v=arrays['native_velocity'];f=arrays['native_scalar'];w=arrays['weights'];masks=arrays['masks']
    assert v.shape==f.shape and len(v)==65 and np.all(np.isfinite(v)) and np.all(np.isfinite(f))
    assert np.all(v[0]==0) and np.all(f[0]==0)
    replay=np.column_stack([np.sqrt(np.sum(w*v*v,axis=1)),np.sqrt(np.sum(w*f*f,axis=1)),
        *[np.sqrt(np.sum(w[m]*v[:,m]**2,axis=1)/w[m].sum()) for m in masks]])
    saved=np.array([[row[key] for key in FIELDS] for row in data['history']])
    errors=np.max(abs(replay-saved),axis=0)/np.maximum(np.max(abs(saved),axis=0),1e-100)
    error=float(max(errors));assert error<1e-12,(p,error);max_error=max(max_error,error)
    energy=arrays['heat_energy'];heat=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    assert heat<2e-13 and abs(heat-data['heat_telescoping'])<1e-15
    assert list(np.array([row['tau'] for row in data['history']]))==list(np.arange(65)/64)
    cases[str(p.relative_to(OUT))]=dict(replay_relative=errors.tolist(),heat_telescoping=heat)


def comparisons(paths):
    rows=[json.loads(p.read_text())['history'] for p in paths];result={}
    for key in FIELDS:
        a,b,c=[np.array([row[key] for row in history]) for history in rows]
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
        result[key]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
    return result


decisions={}
for name in ['factor','balanced','hierarchical']:
    saved=json.loads((OUT/name/'result.json').read_text())
    cmp=comparisons([OUT/name/f'p4-{n}.json' for n in [128,256,512]])
    for key in FIELDS:
        for item in ['previous','last','order']:assert abs(cmp[key][item]-saved['propagation']['4'][key][item])<1e-12
    passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    assert not passed and not saved['passed'];decisions[name]=dict(passed=passed,comparisons=cmp)
cmp=comparisons([OUT/'direct-time'/f'p4-{n}.json' for n in [512,1024,2048]])
saved=json.loads((OUT/'direct-time/stage-result.json').read_text())
for key in FIELDS:
    for item in ['previous','last','order']:assert abs(cmp[key][item]-saved['temporal']['4'][key][item])<1e-12
assert not all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
decisions['direct-time']=dict(passed=False,comparisons=cmp)
# Independent full-system time reference is itself not fully accepted. This
# contrast is diagnostic only and cannot promote any reduced path to success.
reference=json.loads((OUT/'direct-time/p4-2048.json').read_text())['history']
contrast={}
for n in [128,256,512]:
    reduced=json.loads((OUT/'hierarchical'/f'p4-{n}.json').read_text())['history'];contrast[str(n)]={}
    for key in FIELDS:
        a,b=[np.array([row[key] for row in rows]) for rows in [reduced,reference]]
        contrast[str(n)][key]=float(np.max(abs(a-b))/max(abs(b).max(),1e-100))
result=dict(classification='Counterexample candidate',saved_replay_passed=True,scientific_acceptance=False,
    histories=len(cases),source_bindings=bindings,max_native_RMS_replay_relative=max_error,cases=cases,
    decisions=decisions,reduced_vs_direct_unaccepted_reference=contrast,
    original_failure_resolved=False,new_evolutions=0,new_EOS_calls=0)
write(OUT/'replay.json',result)
print(json.dumps({k:v for k,v in result.items() if k not in ['cases','decisions']},indent=2))
