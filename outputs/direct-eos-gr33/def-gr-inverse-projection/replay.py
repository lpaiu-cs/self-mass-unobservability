"""Independent saved-history and unchanged-gate replay; no new evolution."""
from pathlib import Path
import json
import hashlib
import numpy as np

OUT=Path(__file__).resolve().parent
FIELDS=['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']
read=lambda p:json.loads(p.read_text())
count=0;bindings=0;largest=0.;checks={}
for p in sorted(OUT.rglob('*.json')):
    row=read(p)
    for name,sha in row.get('bindings',{}).items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha,(p,name)
        bindings+=1
    if 'history' not in row:continue
    d=np.load(p.with_suffix('.npz'));v=d['native_velocity'];f=d['native_scalar'];w=d['weights'];masks=d['masks']
    assert v.shape==f.shape and len(v)==65 and np.all(np.isfinite(v)) and np.all(np.isfinite(f))
    assert np.all(v[0]==0) and np.all(f[0]==0)
    measured=np.column_stack([np.sqrt(np.sum(w*v*v,axis=1)),np.sqrt(np.sum(w*f*f,axis=1)),
        *[np.sqrt(np.sum(w[m]*v[:,m]**2,axis=1)/w[m].sum()) for m in masks]])
    expected=np.array([[r[k] for k in FIELDS] for r in row['history']])
    error=float(np.max(abs(measured-expected)/np.maximum(abs(expected).max(0),1e-100)))
    assert error<1e-12,(p,error);largest=max(largest,error);count+=1
    assert np.array_equal([r['tau'] for r in row['history']],np.arange(65)/64)
    energy=d['heat_energy'];balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    assert balance<2e-13 and abs(balance-row['heat_telescoping'])<1e-15
    if 'pole_expansion_absolute_bounds' in row:
        for k,bound in row['pole_expansion_absolute_bounds'].items():assert bound<1e-12*max(abs(r[k]) for r in row['history'])
    checks[str(p.relative_to(OUT))]=dict(native_replay_relative=error,heat_telescoping=balance)
decisions={}
for directory,levels in [('',[128,256,512]),('fields',[128,256,512]),('gauss',[512,1024,2048]),('multishift',[128,256,512]),('inverse-gram',[128,256,512])]:
    result=read(OUT/directory/'stage-result.json');rows=[read(OUT/directory/f'p4-{n}.json')['history'] for n in levels];comparisons={}
    for key in FIELDS:
        a,b,c=[np.array([r[key] for r in h]) for h in rows];norm=max(abs(c).max(),1e-100)
        first=np.max(abs(a-b))/norm;last=np.max(abs(b-c))/norm
        comparisons[key]=dict(previous=float(first),last=float(last),order=float(np.log2(first/last)))
        for k,v in comparisons[key].items():assert abs(v-result['comparisons'][key][k])<1e-12
    passed=all(v['last']<.02 and v['order']>1.5 for v in comparisons.values())
    assert passed==result['propagation_passed'] and not passed
    decisions[directory or 'inverse']=dict(passed=passed,comparisons=comparisons)
midpoint=read(OUT.parent/'def-gr-spatial-repair/direct-time/result.json')['temporal']['4']
gauss=decisions['gauss']['comparisons']
gains={key:midpoint[key]['last']/gauss[key]['last'] for key in [FIELDS[0],FIELDS[2]]}
result=dict(classification='Counterexample candidate',saved_replay_passed=True,scientific_acceptance=False,
    histories=count,source_bindings=bindings,maximum_native_replay_relative=largest,checks=checks,decisions=decisions,
    direct_velocity_difference_reduction=gains,original_failure_resolved=False,new_time_paths=0,new_EOS_calls=0)
(OUT/'replay.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['checks','decisions']},indent=2))
