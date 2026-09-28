"""One same2048 path within original cap; no triple after failed forecast."""
import json
import signal
import time
import resource
from pathlib import Path
import numpy as np
import def_gr_gauss_refined as run

OUT=run.OUT
assert not (OUT/'arithmetic-plan.json').exists()
pilot=json.loads((OUT/'pilot-budget.json').read_text());plan=json.loads((OUT/'plan.json').read_text())
forecast=1.4*(pilot['setup_seconds']+2048/64*pilot['steps64_seconds']+2)
run.write(OUT/'arithmetic-plan.json',dict(classification='Counterexample candidate',
    reason='Original three-path forecast exceeded360s. Do not run the triple. One same2048 path separates arithmetic from temporal error at the existing finest path.',
    forecast_seconds=forecast,hard_seconds=300,
    decision='Compare all four full histories with saved same2048 and previous1024/2048 Gauss differences. No convergence acceptance from one path. Do not launch512/1024 or spatial/contrast paths in this branch.',
    bindings={str(p):run.task.digest(p) for p in [Path(__file__),Path(run.__file__),OUT/'plan.json',OUT/'pilot-budget.json',run.common.OUT/'gauss/p4-2048.json',run.common.OUT/'gauss/p4-2048.npz']}))
assert forecast<300
for p,h in plan['bindings'].items():assert run.task.digest(Path(p))==h,p
signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
model=run.space.Model(4,run.task.BANK/'fine-bank.npz');row=run.evolve(model,2048,'p4-2048')
old=json.loads((run.common.OUT/'gauss/p4-2048.json').read_text())
last=json.loads((run.common.OUT/'gauss/stage-result.json').read_text())['comparisons'];comparisons={}
for key in run.task.FIELDS:
    a=np.array([r[key] for r in old['history']]);b=np.array([r[key] for r in row['history']])
    change=float(max(abs(a-b))/max(abs(b).max(),1e-100))
    comparisons[key]=dict(same_step_arithmetic_difference=change,previous_temporal_last=last[key]['last'],arithmetic_to_previous_temporal_ratio=change/last[key]['last'])
result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,comparisons=comparisons,
    seconds=time.monotonic()-start,propagation_passed=False,convergence_triple_executed=False,
    original_failure_resolved=False,full_dynamic_charge_solved=False)
run.write(OUT/'arithmetic-result.json',result);signal.alarm(0);print('ARITHMETIC RESULT',json.dumps(result),flush=True)
