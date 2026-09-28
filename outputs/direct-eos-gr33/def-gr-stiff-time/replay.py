"""Independent native readout, provenance and unchanged-gate replay."""
from pathlib import Path
import json
import hashlib
import numpy as np

ROOT=Path('outputs/direct-eos-gr33');names=['def-gr-stiff-time','def-gr-gauss-refined','def-gr-allocated-projection']
FIELDS=['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']
read=lambda p:json.loads(p.read_text())
checks={};bindings=0
for name in names:
    for p in sorted((ROOT/name).glob('*.json')):
        row=read(p)
        for path,sha in row.get('bindings',{}).items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha,(p,path)
            bindings+=1
        if 'history' not in row:continue
        d=np.load(p.with_suffix('.npz'));v=d['native_velocity'];f=d['native_scalar'];w=d['weights'];masks=d['masks']
        assert len(v)==65 and v.shape==f.shape and np.isfinite(v).all() and np.isfinite(f).all()
        assert np.all(v[0]==0) and np.all(f[0]==0)
        assert np.array_equal([r['tau'] for r in row['history']],np.arange(65)/64)
        actual=np.column_stack([np.sqrt(np.sum(w*v*v,axis=1)),np.sqrt(np.sum(w*f*f,axis=1)),
            *[np.sqrt(np.sum(w[m]*v[:,m]**2,axis=1)/w[m].sum()) for m in masks]])
        expected=np.array([[r[k] for k in FIELDS] for r in row['history']])
        error=float(np.max(abs(actual-expected)/np.maximum(abs(expected).max(0),1e-100)));assert error<1e-12
        energy=d['heat_energy'];balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
        assert balance<2e-13 and abs(balance-row['heat_telescoping'])<1e-15
        assert row.get('linear_residual',0)<1e-9
        for key,bound in row.get('pole_expansion_absolute_bounds',{}).items():assert bound<1e-12*max(abs(r[key]) for r in row['history'])
        checks[str(p.relative_to(ROOT))]=dict(readout_relative=error,heat_balance=balance)
decisions={}
for name,levels in [('def-gr-stiff-time',[512,1024,2048]),('def-gr-allocated-projection',[128,256,512])]:
    result=read(ROOT/name/'stage-result.json');rows=[read(ROOT/name/f'p4-{n}.json')['history'] for n in levels];cmp={}
    for key in FIELDS:
        a,b,c=[np.array([r[key] for r in row]) for row in rows];norm=max(abs(c).max(),1e-100)
        first=float(max(abs(a-b))/norm);last=float(max(abs(b-c))/norm);order=float(np.log2(first/last))
        cmp[key]=dict(previous=first,last=last,order=order)
        for k,v in cmp[key].items():assert abs(v-result['comparisons'][key][k])<1e-12
    passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    if 'saved_Gauss_differences' in result:passed=passed and max(result['saved_Gauss_differences'].values())<.02
    assert passed==result['propagation_passed']
    decisions[name]=dict(passed=passed,comparisons=cmp)
p=ROOT/'def-gr-gauss-refined';result=read(p/'arithmetic-result.json');old=read(ROOT/'def-gr-inverse-projection/gauss/p4-2048.json');new=read(p/'p4-2048.json')
for key in FIELDS:
    a=np.array([r[key] for r in old['history']]);b=np.array([r[key] for r in new['history']])
    change=float(max(abs(a-b))/max(abs(b).max(),1e-100));assert abs(change-result['comparisons'][key]['same_step_arithmetic_difference'])<1e-15
assert not result['convergence_triple_executed'] and not result['propagation_passed']
output=dict(classification='Counterexample candidate',saved_replay_passed=True,histories=len(checks),source_bindings=bindings,
    maximum_readout_relative=max(r['readout_relative'] for r in checks.values()),checks=checks,decisions=decisions,
    arithmetic_comparison_reproduced=True,original_failure_resolved=False,new_evolutions=0)
(ROOT/'def-gr-stiff-time/replay.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({k:v for k,v in output.items() if k not in ['checks','decisions']},indent=2))
