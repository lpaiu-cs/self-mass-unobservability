"""Conditional original input, exterior, quadrature and abscissa contrasts."""
from pathlib import Path
import json
import time
import signal
import resource
import numpy as np
import def_gr_interface_patch as task

go=task.go;OUT=task.OUT;read=lambda p:json.loads(p.read_text());task.verify_plan();task.install()
assert read(OUT/'spatial-result.json')['spatial_passed']
assert not (OUT/'contrasts-plan.json').exists()
assert read(OUT/'operator-check.json')['passed']
used=sum(read(OUT/p)['seconds'] for p in ['pilot-budget.json','p4-result.json','spatial-pilot.json','spatial-result.json','operator-check.json'])
cases=[('coefficient',go.task.BANK/'coarse-bank.npz',2,6,12,.02),
       ('outer',go.task.BANK/'fine-bank.npz',3,6,12,.002),
       ('quadrature',go.task.BANK/'fine-bank.npz',2,8,12,.002),
       ('abscissa',go.task.BANK/'fine-bank.npz',2,6,14,.0002)]
task.write(OUT/'contrasts-plan.json',dict(classification='Counterexample candidate',
    claim='Complete the four original conditional contrasts on the same accepted patched coupled GR response.',
    retained='One fixed candidate per contrast,all65 times and four fields,full native histories and unchanged original thresholds. Every contrast also retains its own propagation/contour gates. Coarse pole bank is a sensitivity contrast, not a physical error certificate.',
    cases=[dict(label=l,bank=str(b),outer=o,quadrature=q,sigma=s,relative_gate=g) for l,b,o,q,s,g in cases],
    budget=dict(pilot_seconds=90,parent_total_compute_seconds=1800,used_measured_seconds=used,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
    decision='Measure4 representative resolvents per candidate. Require forecast for all four within the remaining1800s parent compute budget before production. Stop on the first propagation or contrast failure; no post-result candidate replacement.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(task.__file__),Path(go.__file__),OUT/'spatial-result.json',OUT/'p4-result.json',go.task.BANK/'coarse-bank.npz']}))
signal.alarm(90);pilot_started=time.monotonic();pilots={}
inversion=read(task.PRIOR/'pilot-budget.json')['inversion_forecast_seconds']
for label,bank,outer,quad,sigma,gate in cases:
    t=time.monotonic();p=go.Problem(4,bank,outer,quad);setup=time.monotonic()-t;t=time.monotonic()
    for k in [1,683,1365,2048]:p.transform(go.contour(k,sigma)[0])
    seconds=time.monotonic()-t
    # Same native inversion size and almost identical banded solve; use the
    # actual per-case resolvents with a20% margin, not the initial40% pilot.
    forecast=1.2*(setup+seconds/4*2048+inversion+20)
    pilots[label]=dict(setup_seconds=setup,solve4_seconds=seconds,forecast_seconds=forecast,dofs=p.model.size)
    del p
pilot_seconds=time.monotonic()-pilot_started;remaining=1800-used-pilot_seconds
forecast=sum(p['forecast_seconds'] for p in pilots.values())
task.write(OUT/'contrasts-pilot.json',dict(classification='Counterexample candidate',cases=pilots,seconds=pilot_seconds,
    forecast_seconds=forecast,remaining_compute_seconds=remaining,budget_passed=forecast<remaining,
    assumption='Measured case-specific4 solves scaled to2048 plus prior measured native inversion,20s output and20 percent margin. Four sampled frequencies do not bound worst-case runtime; remaining-budget alarm is authoritative.'))
signal.alarm(0);print('PATCH CONTRAST FORECAST',forecast,'REMAINING',remaining,flush=True)
assert forecast<remaining
signal.alarm(int(remaining));resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));started=time.monotonic()
reference=read(OUT/'p4-2048.json')['history'];results={};passed=True
for label,bank,outer,quad,sigma,gate in cases:
    row=go.solve(bank=bank,outer=outer,quadrature=quad,sigma=sigma,label=label)
    history=read(OUT/f'{label}-2048.json')['history'];cmp={}
    for key in go.task.FIELDS:
        a=np.array([r[key] for r in reference]);b=np.array([r[key] for r in history])
        cmp[key]=float(max(abs(a-b))/max(abs(a).max(),1e-100))
    accepted=row['propagation_passed'] and max(cmp.values())<gate
    results[label]=dict(passed=accepted,propagation_passed=row['propagation_passed'],relative=cmp,gate=gate,seconds=row['seconds'])
    task.write(OUT/f'{label}-contrast.json',dict(classification='Counterexample candidate',**results[label]))
    if not accepted:passed=False;break
result=dict(classification='Counterexample candidate',all_original_contrasts_passed=passed and len(results)==4,
    cases=results,seconds=time.monotonic()-started,parent_measured_compute_seconds=used+pilot_seconds+time.monotonic()-started,
    original_failure_resolved=False,full_dynamic_charge_solved=False)
task.write(OUT/'contrasts-result.json',result);signal.alarm(0);print('PATCH CONTRAST RESULT',json.dumps(result),flush=True)
