"""Reallocate the existing600s envelope before any actual Gauss paths start."""
from pathlib import Path
import time
import signal
import resource
import json
import def_gr_gauss_time as method

OUT=method.OUT;task=method.task
pilot=json.loads((OUT/'pilot-budget.json').read_text());assert pilot['first_case_forecast_seconds']<240
plan=json.loads((OUT/'plan.json').read_text())
for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
task.write(OUT/'budget-reallocation-plan.json',dict(classification='Counterexample candidate',
    reason='Measured complex direct steps predict215.81s including40 percent margin, above initial180s first-case allocation. Reallocate240s of the unchanged600s total envelope to the decisive degree4 path. No added paths, steps, degree, horizon or relaxed accuracy gate; remaining contrasts keep a600s cumulative cap and require another measured-cost check.',
    first_case_forecast_seconds=pilot['first_case_forecast_seconds'],first_case_cap_seconds=240,total_cap_seconds=600,
    decision='Run exactly512/1024/2048 on degree4. Stop on any of the four unchanged propagation gates. No repeat of this path.',
    bindings={str(p):task.digest(p) for p in [Path(__file__),Path(method.__file__),OUT/'plan.json',OUT/'pilot-budget.json']}))
signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic()
model=method.space.Model(4,task.BANK/'fine-bank.npz');rows=[method.evolve(model,n,f'p4-{n}') for n in [512,1024,2048]]
cmp=method.common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
    comparisons=cmp,seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
    original_failure_resolved=False,full_dynamic_charge_solved=False)
task.write(OUT/'stage-result.json',result);signal.alarm(0);print('GAUSS RESULT',json.dumps(result),flush=True)
