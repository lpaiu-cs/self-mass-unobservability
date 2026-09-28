"""Original degree1/2/4 spatial gates on the fixed, source-preserving patch."""
from pathlib import Path
import json
import signal
import time
import resource
import numpy as np
import def_gr_interface_patch as task

go=task.go;OUT=task.OUT;read=lambda p:json.loads(p.read_text())
assert read(OUT/'p4-result.json')['propagation_passed'];task.verify_plan();task.install()
assert not (OUT/'spatial-plan.json').exists()
task.write(OUT/'spatial-plan.json',dict(classification='Counterexample candidate',
    claim='Apply the original1/2/4 polynomial-space acceptance to the single frozen local geometry repair, preserving all four native readouts.',
    gates=dict(propagation_relative=.02,propagation_order=1.5,contour=.0002,spatial_relative=.02,spatial_decrease=True),
    budget=dict(pilot_seconds=60,production_seconds=300,parent_total_compute_seconds=1800,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
    decision='Measure4 resolvents on each remaining space; run only within300s forecast. Stop on a time failure or the final spatial failure. No further cell split or degree increase.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(task.__file__),Path(go.__file__),OUT/'plan.json',OUT/'p4-result.json']}))
signal.alarm(60);pilot_start=time.monotonic();pilots={};old=read(task.PRIOR/'pilot-budget.json')
for degree in [2,1]:
    start=time.monotonic();p=go.Problem(degree);setup=time.monotonic()-start;start=time.monotonic()
    for k in [1,683,1365,2048]:p.transform(go.contour(k,12)[0])
    seconds=time.monotonic()-start;forecast=1.3*(setup+seconds/4*2048+old['inversion_forecast_seconds']+15)
    pilots[degree]=dict(setup_seconds=setup,solve4_seconds=seconds,forecast_seconds=forecast,dofs=p.model.size)
    del p
forecast=sum(p['forecast_seconds'] for p in pilots.values())
task.write(OUT/'spatial-pilot.json',dict(classification='Counterexample candidate',cases=pilots,forecast_seconds=forecast,
    seconds=time.monotonic()-pilot_start,assumption='Measured setups and4 resolvents per space plus saved native inversion,15s output and30 percent margin.'))
signal.alarm(0);assert forecast<300;print('PATCH SPATIAL FORECAST',forecast,flush=True)
signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic();passed=True;paths=[]
for degree in [2,1]:
    row=go.solve(degree=degree,label=f'p{degree}');paths.append(degree)
    if not row['propagation_passed']:passed=False;break
cmp={}
if passed:
    rows=[read(OUT/f'p{p}-2048.json')['history'] for p in [1,2,4]]
    for key in go.task.FIELDS:
        a,b,c=[np.array([r[key] for r in row]) for row in rows];norm=max(abs(c).max(),1e-100)
        first=float(max(abs(a-b))/norm);last=float(max(abs(b-c))/norm)
        cmp[key]=dict(previous=first,last=last,decreased=last<first)
    passed=all(v['last']<.02 and v['decreased'] for v in cmp.values())
result=dict(classification='Counterexample candidate',spatial_passed=passed,comparisons=cmp,paths=paths,
    seconds=time.monotonic()-start,original_failure_resolved=False,full_dynamic_charge_solved=False)
task.write(OUT/'spatial-result.json',result);signal.alarm(0);print('PATCH SPATIAL RESULT',json.dumps(result),flush=True)
