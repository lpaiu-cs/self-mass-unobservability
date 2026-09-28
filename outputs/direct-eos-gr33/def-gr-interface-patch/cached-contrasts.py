"""Run the unchanged four contrasts after validated invariant-data caching."""
from pathlib import Path
import json
import signal
import time
import resource
import numpy as np
import def_gr_interface_patch as task
import def_gr_cached_resolvent as cached

go=task.go;OUT=task.OUT;read=lambda p:json.loads(p.read_text());task.verify_plan();task.install()
pilot=read(OUT/'cached-pilot.json');assert pilot['budget_passed'] and pilot['transfer_equivalence_passed']
assert read(OUT/'spatial-result.json')['spatial_passed'] and read(OUT/'operator-check.json')['passed']
for row in pilot['checks']:
    values=np.array(list(row['relative'].values()));assert np.isfinite(values).all() and max(values)<1e-10
for name in ['contrasts-plan.json','cached-plan.json']:
    for p,h in read(OUT/name)['bindings'].items():assert go.task.digest(Path(p))==h,p
assert not (OUT/'cached-execution-plan.json').exists() and not (OUT/'contrasts-result.json').exists()
remaining=pilot['remaining_compute_seconds'];cases=read(OUT/'contrasts-plan.json')['cases']
task.write(OUT/'cached-execution-plan.json',dict(classification='Counterexample candidate',
    claim='Complete exactly the four previously planned contrasts with validated invariant-data caching and unchanged numerical corrections.',
    forecast_seconds=pilot['forecast_seconds'],remaining_compute_seconds=remaining,parent_total_compute_seconds=1800,
    decision='Stop the first original propagation or contrast failure; stop at the remaining-budget alarm. No new geometry,parameter search or threshold change.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(cached.__file__),Path(task.__file__),Path(go.__file__),OUT/'cached-pilot.json',OUT/'contrasts-plan.json']}))
go.Problem=cached.Problem
signal.alarm(int(remaining));resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));started=time.monotonic()
reference=read(OUT/'p4-2048.json')['history'];results={};passed=True
for case in cases:
    label=case['label'];gate=case['relative_gate']
    row=go.solve(bank=Path(case['bank']),outer=case['outer'],quadrature=case['quadrature'],sigma=case['sigma'],label=label)
    history=read(OUT/f'{label}-2048.json')['history'];cmp={}
    for key in go.task.FIELDS:
        a=np.array([r[key] for r in reference]);b=np.array([r[key] for r in history])
        cmp[key]=float(max(abs(a-b))/max(abs(a).max(),1e-100))
    accepted=row['propagation_passed'] and max(cmp.values())<gate
    results[label]=dict(passed=accepted,propagation_passed=row['propagation_passed'],relative=cmp,gate=gate,seconds=row['seconds'])
    task.write(OUT/f'{label}-contrast.json',dict(classification='Counterexample candidate',**results[label]))
    print('ORIGINAL CONTRAST',label,json.dumps(results[label]),flush=True)
    if not accepted:passed=False;break
seconds=time.monotonic()-started
result=dict(classification='Counterexample candidate',all_original_contrasts_passed=passed and len(results)==4,
    cases=results,seconds=seconds,parent_accounted_compute_seconds=1800-remaining+seconds,
    original_failure_resolved=False,full_dynamic_charge_solved=False)
task.write(OUT/'contrasts-result.json',result);signal.alarm(0);print('PATCH CONTRAST RESULT',json.dumps(result),flush=True)
