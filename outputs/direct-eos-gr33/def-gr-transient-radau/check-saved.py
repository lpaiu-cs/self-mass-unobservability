"""Replay both rejected GR candidates and the precision control from disk."""
from pathlib import Path
import json
import numpy as np
import def_gr_transient_radau as work

out=work.OUT
def read(p):return json.loads(p.read_text())
counts=[]
for directory in [out,out/'consistent-inertia']:
    for p,h in read(directory/'plan.json')['bindings'].items():assert work.digest(Path(p))==h,p
    result=read(directory/'result.json');hist={n:read(directory/f'heat-{n}.json')['history'] for n in [16,32,64]}
    for field,expected in result['comparisons'].items():
        a,b,c=[np.array([x[field] for x in hist[n]]) for n in [16,32,64]];norm=max(abs(c).max(),1e-100)
        e1=float(max(abs(a-b[::2]))/norm);e2=float(max(abs(b-c[::2]))/norm)
        assert e1==expected['time_previous'] and e2==expected['time_last']
        assert float(np.log2(e1/e2))==expected['order']
        for name,label in [('coarse','coefficient-64'),('outer','outer-64')]:
            v=np.array([x[field] for x in read(directory/(label+'.json'))['history']])
            key='coefficients' if name=='coarse' else name
            assert float(max(abs(c-v))/norm)==expected[key]
    assert not result.get('passed',result.get('numerical_gates_passed'))
    assert result['max_heat_balance']<2e-13 and result['max_linear_residual']<1e-9
    counts.append(5)
old=read(work.prior.OUT/'heat-64.json')['history'];new=read(out/'precision-64.json')['history']
for field,expected in read(out/'precision.json')['differences'].items():
    a=np.array([x[field] for x in old]);b=np.array([x[field] for x in new])
    assert float(max(abs(a-b))/max(abs(a).max(),1e-100))==expected
assert max(read(out/'precision.json')['differences'].values())<1e-4
assert read(out/'control.json')['symbolic_passed']
assert min(read(out/'control.json')['orders'])>4.7
assert read(out/'consistent-inertia/symbolic.json')['passed']
assert not read(out/'consistent-inertia/rejection.json')['candidate_physically_adopted']
work.write(out/'verification.json',dict(classification='Counterexample candidate',saved_replay_passed=True,
    actual_paths_replayed=sum(counts)+1,reproduced_failed_candidates=2,new_evolution_paths=0,
    scope='Confirms precision comparison and both unchanged failed verdicts, not successful GR closure.'))
print('PASS replay of11 actual paths; both candidate failures remain false')
