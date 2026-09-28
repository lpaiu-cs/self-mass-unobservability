"""Conditional fixed degree2/1 paths after common degree4 acceptance."""
from pathlib import Path
import json
import time
import signal
import resource
import numpy as np
import def_gr_full_weeks as go

ROOT=go.OUT;OUT=ROOT/'beta1024';go.OUT=OUT;go.BETA=1024
assert json.loads((OUT/'p4-result.json').read_text())['propagation_passed']
assert json.loads((OUT/'scale-contrast.json').read_text())['scale_contrast_passed']
assert not (OUT/'spatial-plan.json').exists()
go.write(OUT/'spatial-plan.json',dict(classification='Counterexample candidate',
    claim='Complete the original1/2/4 polynomial-space comparison using the now jointly accepted full-matrix time inversion.',
    retained='Same native readouts,heat poles,zero initial state,horizon,beta1024,sigma12 and4096 contour nodes. No finer degree or grid. Both degree2 and1 require their own unchanged four-field propagation and contour gates.',
    gates=dict(propagation_relative=.02,propagation_order=1.5,contour=.0002,spatial_relative=.02,spatial_decrease=True),
    budget=dict(pilot_seconds=60,spatial_seconds=300,parent_total_seconds=1200,parent_reserved_for_prior_work_seconds=650,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
    decision='Measure4 representative resolvents on each space. Run only if combined forecast fits300s; stop a failed propagation. Compare degree1,2,4 only if both additional spaces pass. Do not refine space on failure.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(go.__file__),OUT/'plan.json',OUT/'p4-result.json',OUT/'scale-contrast.json',go.task.BANK/'fine-bank.npz']}))
lu=go.splu;go.splu=lambda A:lu(A,permc_spec='NATURAL')
signal.alarm(60);pilots={};old=json.loads((ROOT/'pilot-budget.json').read_text())
for degree in [2,1]:
    start=time.monotonic();p=go.Problem(degree);setup=time.monotonic()-start;start=time.monotonic()
    for k in [1,683,1365,2048]:p.transform(go.contour(k,12)[0])
    elapsed=time.monotonic()-start
    forecast=1.3*(setup+2048/4*elapsed+old['inversion_forecast_seconds']+15)
    pilots[degree]=dict(setup_seconds=setup,resolvent4_seconds=elapsed,forecast_seconds=forecast,dofs=p.model.size)
    del p
total=sum(r['forecast_seconds'] for r in pilots.values())
go.write(OUT/'spatial-pilot.json',dict(classification='Counterexample candidate',cases=pilots,forecast_seconds=total,
    assumption='Measured smaller matrices; same native inversion cost plus15s output per case and30 percent allowance. Full contours at these degrees unmeasured.'))
signal.alarm(0);print('SPATIAL FORECAST',total,flush=True)
assert total<300
signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic();accepted=True;paths=[]
for degree in [2,1]:
    row=go.solve(degree=degree,label=f'p{degree}');paths.append(degree)
    if not row['propagation_passed']:accepted=False;break
cmp={}
if accepted:
    rows=[json.loads((OUT/f'p{degree}-2048.json').read_text())['history'] for degree in [1,2,4]]
    for key in go.task.FIELDS:
        a,b,c=[np.array([r[key] for r in row]) for row in rows];norm=max(abs(c).max(),1e-100)
        first=float(max(abs(a-b))/norm);last=float(max(abs(b-c))/norm)
        cmp[key]=dict(previous=first,last=last,decreased=last<first)
    accepted=all(v['last']<.02 and v['decreased'] for v in cmp.values())
result=dict(classification='Counterexample candidate',spatial_passed=accepted,comparisons=cmp,paths=paths,
    seconds=time.monotonic()-start,original_failure_resolved=False,full_dynamic_charge_solved=False)
go.write(OUT/'spatial-result.json',result);signal.alarm(0);print('SPATIAL RESULT',json.dumps(result),flush=True)
