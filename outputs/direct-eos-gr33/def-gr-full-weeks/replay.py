"""Replay saved native fields, immutable gates and the remaining spatial error."""
from pathlib import Path
import hashlib
import json
import signal
import time
import numpy as np

ROOT=Path('outputs/direct-eos-gr33');OUT=ROOT/'def-gr-full-weeks';B=OUT/'beta1024'
FIELDS=['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']
read=lambda p:json.loads(p.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
signal.alarm(60);started=time.monotonic();checks={};bindings=0
for folder in ['def-gr-transfer-repair','def-gr-response-basis','def-gr-full-weeks']:
    for p in sorted((ROOT/folder).rglob('*.json')):
        row=read(p)
        for path,h in row.get('bindings',{}).items():
            assert sha(Path(path))==h,(p,path);bindings+=1
        if 'history' not in row:continue
        d=np.load(p.with_suffix('.npz'));v=d['native_velocity'];f=d['native_scalar'];w=d['weights'];masks=d['masks']
        assert len(v)==65 and v.shape==f.shape and np.isfinite(v).all() and np.isfinite(f).all()
        assert np.all(v[0]==0) and np.all(f[0]==0)
        assert np.array_equal([r['tau'] for r in row['history']],np.arange(65)/64)
        actual=np.column_stack([np.sqrt(np.sum(w*v*v,axis=1)),np.sqrt(np.sum(w*f*f,axis=1)),
            *[np.sqrt(np.sum(w[m]*v[:,m]**2,axis=1)/w[m].sum()) for m in masks]])
        expected=np.array([[r[k] for k in FIELDS] for r in row['history']])
        error=float(np.max(abs(actual-expected)/np.maximum(abs(expected).max(0),1e-100)));assert error<1e-12
        E=d['heat_energy'];balance=float(abs(np.sum(-np.diff(E),dtype=np.longdouble))/max(abs(E).max(),1e-100))
        assert balance<2e-13 and abs(balance-row['heat_telescoping'])<1e-15
        assert row.get('linear_residual',0)<1e-9
        checks[str(p.relative_to(ROOT))]=dict(readout_relative=error,heat_balance=balance)

def compare(rows):
    result={}
    for key in FIELDS:
        a,b,c=[np.array([r[key] for r in row['history']]) for row in rows];norm=max(abs(c).max(),1e-100)
        first=float(max(abs(a-b))/norm);last=float(max(abs(b-c))/norm)
        result[key]=dict(previous=first,last=last,order=float(np.log2(first/last)),decreased=last<first)
    return result

decisions={}
for folder,prefix,levels,result_name in [(ROOT/'def-gr-response-basis','p4',[128,256,512],'stage-result'),
        (OUT,'p4',[512,1024,2048],'p4-result'),*[(B,f'p{p}',[512,1024,2048],f'p{p}-result') for p in [1,2,4]]]:
    result=read(folder/(result_name+'.json'));cmp=compare([read(folder/f'{prefix}-{n}.json') for n in levels])
    for key,row in cmp.items():
        for k in ['previous','last','order']:assert abs(row[k]-result['comparisons'][key][k])<1e-12
    passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    passed=passed and all(v.get('contour',0)<.0002 for v in result['comparisons'].values())
    assert passed==result['propagation_passed']
    decisions[str((folder/prefix).relative_to(ROOT))]=dict(passed=passed,comparisons=cmp)
cmp=compare([read(B/f'p{p}-2048.json') for p in [1,2,4]]);result=read(B/'spatial-result.json')
for key,row in cmp.items():
    for k in ['previous','last','decreased']:assert row[k]==result['comparisons'][key][k]
assert result['spatial_passed']==all(v['last']<.02 and v['decreased'] for v in cmp.values())

# Reuse the actual failed histories; no new solve or altered readout.
a=np.load(B/'p2-2048.npz');b=np.load(B/'p4-2048.npz');mask=b['masks'][1];w=b['weights'];radius=b['native_radius']
for key in ['native_radius','weights','masks','cells','heat_energy','heat_flux']:assert np.array_equal(a[key],b[key]),key
history=[read(B/f'p{p}-2048.json')['history'] for p in [2,4]]
diff=np.array([x[FIELDS[3]]-y[FIELDS[3]] for x,y in zip(*history)]);j=int(np.argmax(abs(diff)))
delta=b['native_velocity'][j]-a['native_velocity'][j];squared=w*delta**2;total=squared[mask].sum()
ids=np.flatnonzero(mask);rank=ids[np.argsort(squared[ids])[::-1]];cells=b['cells'];cids=np.clip(np.searchsorted(cells,radius,side='right')-1,0,len(cells)-2)
top=[dict(native_index=int(i),radius_R=float(radius[i]),cell=int(cids[i]),cell_width_R=float(cells[cids[i]+1]-cells[cids[i]]),
          cell_coordinate=float((radius[i]-cells[cids[i]])/(cells[cids[i]+1]-cells[cids[i]])),
          p2_velocity_m_s=float(a['native_velocity'][j,i]),p4_velocity_m_s=float(b['native_velocity'][j,i]),
          share_weighted_squared_difference=float(squared[i]/total)) for i in rank[:12]]
location=dict(classification='Counterexample candidate',tau=j/64,native_samples_in_mask=int(mask.sum()),
    rms_peak_tau=float(np.argmax([r[FIELDS[3]] for r in history[1]])/64),
    p4_peak_rms_m_s=max(r[FIELDS[3]] for r in history[1]),
    top12_share_weighted_squared_difference=float(squared[rank[:12]].sum()/total),top=top,
    interpretation='Location of the saved p2/p4 discrepancy, not a causal explanation or a true-error bound. RMS difference gate is unchanged.',new_evolutions=0)
(OUT/'spatial-location.json').write_text(json.dumps(location,indent=2)+'\n')
output=dict(classification='Counterexample candidate',saved_replay_passed=True,histories=len(checks),source_bindings=bindings,
    maximum_readout_relative=max(r['readout_relative'] for r in checks.values()),checks=checks,decisions=decisions,
    spatial_comparison_reproduced=True,propagation_passed_all_three_spaces=True,spatial_passed=False,
    original_failure_resolved=False,new_evolutions=0,seconds=time.monotonic()-started)
(OUT/'replay.json').write_text(json.dumps(output,indent=2)+'\n');signal.alarm(0)
print(json.dumps({k:v for k,v in output.items() if k not in ['checks','decisions']},indent=2));print(json.dumps(location,indent=2))
