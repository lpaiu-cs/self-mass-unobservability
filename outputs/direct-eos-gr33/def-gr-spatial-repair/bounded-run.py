"""Cost re-evaluation: execute/reuse the single worst planned spatial case first."""
from pathlib import Path
import json
import time
import signal
import resource
import numpy as np
import def_gr_spatial_repair as task

OUT=task.OUT;assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']>300
task.write(OUT/'budget-reassessment-plan.json',dict(classification='Counterexample candidate',
    reason='The64-vector pilot extrapolated every cost quadratically and forecast992.86s, exceeding300s; the original production launch was not executed.',
    decision='Measure the single most expensive already planned degree4 case under60s, reuse its results, then permit the remaining unchanged comparisons only if six times measured cost with25 percent allowance is below300s. No increased total cap, degree, time horizon, precision gate or source.',
    budget=dict(first_case_hard_seconds=60,total_hard_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
    bindings={str(p):task.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'plan.json',OUT/'pilot-budget.json']}))
signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic();cases={};propagation={};checks=[]
for degree in [4,2,1]:
    began=time.monotonic();model=task.Model(degree,task.BANK/'fine-bank.npz');P=task.Projection(model)
    _,Q,orth=P.basis(512);K,skew=model.energy_matrix(Q)
    rows={n:task.series(model,K[:n,:n],Q[:,:n],f'p{degree}-{n}') for n in [128,256,512]}
    cases[str(degree)]=rows[512];propagation[str(degree)]={}
    for field in task.FIELDS:
        a,b,c=[np.array([h[field] for h in rows[n]['history']]) for n in [128,256,512]]
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
        propagation[str(degree)][field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
    elapsed=time.monotonic()-began
    checks.append(dict(degree=degree,dofs=model.size,quadrature_points=len(model.weights),seconds=elapsed,orthogonality=orth,skew=skew,linear_residual=P.error))
    del model,P,Q,K
    if degree==4:
        forecast=6*elapsed*1.25
        task.write(OUT/'measured-budget.json',dict(classification='Counterexample candidate',worst_planned_case_seconds=elapsed,
            forecast_total_seconds=forecast,reused_case='p4-128/256/512',
            assumption='All six cases charged at the measured most expensive polynomial degree with25 percent allowance; larger exterior and8-point quadrature still unmeasured. Hard300s unchanged.'))
        print('MEASURED FORECAST',forecast,flush=True)
        assert forecast<300,'Budget exceeded; retain this case and do not launch remaining cases.'
        signal.alarm(max(1,int(300-(time.monotonic()-start))))
comparisons={}
for field in task.FIELDS:
    a,b,c=[np.array([h[field] for h in cases[str(n)]['history']]) for n in [1,2,4]]
    norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
    comparisons[field]=dict(previous=float(d1),last=float(d2),decreased=bool(d2<d1))
prop_pass=all(v['last']<.02 and v['order']>1.5 for p in propagation.values() for v in p.values())
spatial_pass=all(v['last']<.02 and v['decreased'] for v in comparisons.values())
if prop_pass and spatial_pass:
    for label,bank,outer,nquad in [('coefficient','coarse-bank.npz',2,6),('outer','fine-bank.npz',3,6),('quadrature','fine-bank.npz',2,8)]:
        began=time.monotonic();model=task.Model(4,task.BANK/bank,outer,nquad);P=task.Projection(model);_,Q,orth=P.basis(512);K,skew=model.energy_matrix(Q)
        cases[label]=task.series(model,K,Q,label+'-p4-512')
        for field in task.FIELDS:
            c=np.array([h[field] for h in cases['4']['history']]);d=np.array([h[field] for h in cases[label]['history']])
            comparisons[field][label]=float(np.max(abs(c-d))/max(abs(c).max(),1e-100))
        checks.append(dict(label=label,seconds=time.monotonic()-began,orthogonality=orth,skew=skew,linear_residual=P.error));del model,P,Q,K
passed=prop_pass and spatial_pass and all(v.get('coefficient',1)<.02 and v.get('outer',1)<.002 and v.get('quadrature',1)<.002 for v in comparisons.values())
result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,passed=passed,
    propagation_passed=prop_pass,spatial_passed=spatial_pass,propagation=propagation,comparisons=comparisons,checks=checks,
    paths=list(cases),seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
    original_failure_resolved=False,full_dynamic_charge_solved=False)
task.write(OUT/'result.json',result);signal.alarm(0);print('RESULT',json.dumps(result),flush=True)
