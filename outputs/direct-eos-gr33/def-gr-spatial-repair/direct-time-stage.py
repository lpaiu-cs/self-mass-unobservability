"""Budget-gated first actual degree4 path; never repeat it in later contrasts."""
from pathlib import Path
import time
import signal
import json
import numpy as np
import def_gr_direct_time as direct

out=direct.OUT;task=direct.task
pilot=json.loads((out/'pilot-budget.json').read_text())
forecast=1.4*(pilot['setup_seconds']+3584/128*pilot['steps128_seconds']+4)
assert forecast<120
task.write(out/'staged-budget-plan.json',dict(classification='Counterexample candidate',
    reason='The original simultaneous six-case forecast364.58s exceeds300s. Run only the decisive degree4 full coupled path first; do not launch the unmeasured contrasts. Reassess remaining cost from this actual path without repeating it.',
    first_case_forecast_seconds=forecast,first_case_cap_seconds=120,total_cap_seconds=300,
    gates='Unchanged512/1024/2048 steps, relative2%, order1.5 for all four original histories; stop immediately on failure.',
    bindings={str(p):task.digest(p) for p in [Path(__file__),Path(direct.__file__),out/'plan.json',out/'pilot-budget.json']}))
signal.alarm(120);start=time.monotonic()
model=direct.space.Model(4,task.BANK/'fine-bank.npz')
rows={n:direct.evolve(model,n,f'p4-{n}') for n in [512,1024,2048]}
comparison={}
for field in task.FIELDS:
    a,b,c=[np.array([r[field] for r in rows[n]['history']]) for n in [512,1024,2048]]
    norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
    comparison[field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
passed=all(v['last']<.02 and v['order']>1.5 for v in comparison.values())
result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,
    time_passed=passed,passed=False,spatial_passed=False,temporal={'4':comparison},
    paths=['4'],seconds=time.monotonic()-start,original_failure_resolved=False,full_dynamic_charge_solved=False)
task.write(out/'stage-result.json',result)
signal.alarm(0);print('STAGED DIRECT',json.dumps(result),flush=True)
