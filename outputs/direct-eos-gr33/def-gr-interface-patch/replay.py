"""Independent replay of every four-field time, space and conditional gate."""
from pathlib import Path
import hashlib
import json
import time
import signal
import numpy as np

OUT=Path('outputs/direct-eos-gr33/def-gr-interface-patch')
FIELDS=['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']
read=lambda p:json.loads(p.read_text());signal.alarm(60);start=time.monotonic();checks={};bindings=0
reference=np.load(OUT/'p4-2048.npz')
for p in sorted(OUT.glob('*.json')):
    row=read(p)
    for path,h in row.get('bindings',{}).items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h,(p,path);bindings+=1
    if 'history' not in row:continue
    d=np.load(p.with_suffix('.npz'));v=d['native_velocity'];f=d['native_scalar'];w=d['weights'];masks=d['masks']
    assert len(v)==65 and v.shape==f.shape and np.isfinite(v).all() and np.isfinite(f).all()
    assert np.all(v[0]==0) and np.all(f[0]==0)
    assert np.array_equal([r['tau'] for r in row['history']],np.arange(65)/64)
    for key in ['native_radius','weights','masks']:assert np.array_equal(d[key],reference[key]),(p,key)
    actual=np.column_stack([np.sqrt(np.sum(w*v*v,axis=1)),np.sqrt(np.sum(w*f*f,axis=1)),
        *[np.sqrt(np.sum(w[m]*v[:,m]**2,axis=1)/w[m].sum()) for m in masks]])
    expected=np.array([[r[k] for k in FIELDS] for r in row['history']]);assert np.all(abs(expected).max(0)>1e-100)
    error=float(np.max(abs(actual-expected)/abs(expected).max(0)));assert error<1e-12
    E=d['heat_energy'];balance=float(abs(np.sum(-np.diff(E),dtype=np.longdouble))/max(abs(E).max(),1e-100))
    assert balance<2e-13 and abs(balance-row['heat_telescoping'])<1e-15
    assert row['linear_residual']<1e-9
    checks[p.name]=dict(readout_relative=error,heat_balance=balance)

def series(label):return np.array([[r[k] for k in FIELDS] for r in read(OUT/(label+'.json'))['history']])

decisions={}
labels=['p1','p2','p4','coefficient','outer','quadrature','abscissa']
for label in labels:
    a,b,c=[series(f'{label}-{n}') for n in [512,1024,2048]];scale=abs(c).max(0)
    first=np.max(abs(a-b),axis=0)/scale;last=np.max(abs(b-c),axis=0)/scale;orders=np.log2(first/last)
    saved=read(OUT/f'{label}-result.json')
    for j,k in enumerate(FIELDS):
        for name,value in [('previous',first[j]),('last',last[j]),('order',orders[j])]:
            assert abs(value-saved['comparisons'][k][name])<1e-12,(label,k,name)
    passed=bool(np.all(last<.02)&np.all(orders>1.5)) and all(v['contour']<.0002 for v in saved['comparisons'].values())
    assert passed==saved['propagation_passed'];decisions[label]=passed
a,b,c=[series(f'p{p}-2048') for p in [1,2,4]];scale=abs(c).max(0)
first=np.max(abs(a-b),axis=0)/scale;last=np.max(abs(b-c),axis=0)/scale
spatial=bool(np.all(last<.02)&np.all(last<first));saved=read(OUT/'spatial-result.json')
assert spatial==saved['spatial_passed']
for j,k in enumerate(FIELDS):
    assert abs(first[j]-saved['comparisons'][k]['previous'])<1e-15
    assert abs(last[j]-saved['comparisons'][k]['last'])<1e-15
contrasts={};result=read(OUT/'contrasts-result.json')
for label,gate in [('coefficient',.02),('outer',.002),('quadrature',.002),('abscissa',.0002)]:
    delta=np.max(abs(series(label+'-2048')-c),axis=0)/scale;saved=read(OUT/f'{label}-contrast.json')
    for j,k in enumerate(FIELDS):assert abs(delta[j]-saved['relative'][k])<1e-15
    passed=bool(np.all(delta<gate)) and decisions[label];assert passed==saved['passed']
    contrasts[label]=passed
passed=all(decisions.values()) and spatial and all(contrasts.values()) and read(OUT/'operator-check.json')['passed']
assert result['all_original_contrasts_passed']==all(contrasts.values())
row=dict(classification='Counterexample candidate',saved_replay_passed=True,histories=len(checks),source_bindings=bindings,
    maximum_readout_relative=max(v['readout_relative'] for v in checks.values()),checks=checks,
    propagation=decisions,spatial_passed=spatial,contrasts=contrasts,all_four_component_gates_passed=passed,
    seconds=time.monotonic()-start,new_evolutions=0,
    scope='The original declared gates for this frozen-source linear GR model; no claim of a true-error bound, full nonlinear stellar evolution or observational closure.')
(OUT/'replay.json').write_text(json.dumps(row,indent=2)+'\n');signal.alarm(0)
print(json.dumps({k:v for k,v in row.items() if k!='checks'},indent=2))
