"""Reuse measured broad-band build cost; preserve the120s first-case cap."""
from pathlib import Path
import time
import signal
import resource
import json
import numpy as np
import def_gr_inverse_gram as method

OUT=method.OUT;task=method.task
pilot=json.loads((OUT/'pilot-budget.json').read_text());prior=json.loads((OUT.parent/'multishift/stage-result.json').read_text())
oldpilot=json.loads((OUT.parent/'multishift/pilot-budget.json').read_text())
forecast=1.4*(pilot['setup_seconds']+prior['seconds']+64*max(0,pilot['build64_seconds']-oldpilot['build64_seconds'])+8*pilot['response64_seconds'])
assert forecast<120
plan=json.loads((OUT/'plan.json').read_text())
for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
task.write(OUT/'measured-budget-plan.json',dict(classification='Counterexample candidate',
    reason='The naive quadratic pilot forecast131.89s charges fixed setup64 times. Reuse the measured full512 multishift path21.67s, charge the entire new setup again, and scale only the new Gram overhead quadratically, plus response and40 percent margin. No increase of120s first-case or300s total cap.',
    first_case_forecast_seconds=forecast,first_case_cap_seconds=120,total_cap_seconds=300,
    decision='Exactly degree4 at128/256/512; same propagation gates and stop on failure. Check signed-face pole amplitudes before using the absolute remainder aggregation.',
    bindings={str(p):task.digest(p) for p in [Path(__file__),Path(method.__file__),OUT/'plan.json',OUT/'pilot-budget.json',OUT.parent/'multishift/stage-result.json',OUT.parent/'multishift/pilot-budget.json']}))
signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
model=method.space.Model(4,task.BANK/'fine-bank.npz')
assert np.all(np.all(model.heat.amplitude>=0,axis=1)|np.all(model.heat.amplitude<=0,axis=1))
R,F,Q,S,tail,meta=method.basis(model,512)
rows=[method.series(model,R[:n,:n],F[:n],Q[:,:n],S[:n],(tail[0][:,:n],tail[1],tail[2]),f'p4-{n}') for n in [128,256,512]]
cmp=method.common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
    comparisons=cmp,seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
    pole_amplitude_sign_consistency=True,original_failure_resolved=False,full_dynamic_charge_solved=False,**meta)
task.write(OUT/'stage-result.json',result);signal.alarm(0);print('GRAM RESULT',json.dumps(result),flush=True)
