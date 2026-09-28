"""Test cached full transfers and remeasure the remaining fixed contrasts."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_gr_interface_patch as task
import def_gr_cached_resolvent as cached

go=task.go;OUT=task.OUT;task.install();signal.alarm(90);started=time.monotonic()
assert not (OUT/'cached-pilot.json').exists()
task.write(OUT/'cached-plan.json',dict(classification='Counterexample candidate',
    claim='Remove repeated absolute sparse copies,rate casts and invariant sparse-pattern assembly from the identical full resolvent; preserve exact source formula,three original residual corrections and all scientific gates.',
    evidence='16-call profiler:0.226s sparse absolute copies,0.109s sparse diagonal products and0.035s sparse addition among1.716s total. Original contrast forecast1495.92s exceeded1285.38s remaining; no contrast trajectory was launched.',
    check='Compare original and cached full solutions at4 fixed frequencies for each of the four original contrast cases; all four weighted field errors below1e-10. Keep one known16-point original pilot comparison on the baseline. Fixed3 repeated128-column reconstruction timings,median used in forecast.',
    budget=dict(hard_seconds=90,new_EOS_calls=0,new_time_paths=0,parent_total_compute_seconds=1800),
    decision='Only if equivalence passes and measured four-case forecast fits the remaining original budget launch the unchanged four contrasts. Do not change accuracy or physical input.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(cached.__file__),Path(go.__file__),Path(task.__file__),OUT/'contrasts-pilot.json']}))
cases=[('baseline',go.task.BANK/'fine-bank.npz',2,6,12),('coefficient',go.task.BANK/'coarse-bank.npz',2,6,12),
       ('outer',go.task.BANK/'fine-bank.npz',3,6,12),('quadrature',go.task.BANK/'fine-bank.npz',2,8,12),('abscissa',go.task.BANK/'fine-bank.npz',2,6,14)]
pilots={};checks=[]
for label,bank,outer,quad,sigma in cases:
    t=time.monotonic();p=cached.Problem(4,bank,outer,quad);setup=time.monotonic()-t;elapsed=0.
    ids=np.load(OUT/'pilot.npz')['ids'] if label=='baseline' else np.array([1,683,1365,2048])
    saved=np.load(OUT/'pilot.npz')['answers'] if label=='baseline' else None
    for j,k in enumerate(ids):
        z=go.contour(k,sigma)[0];original=saved[j] if saved is not None else cached.OriginalProblem.transform(p,z)
        t=time.monotonic();answer=p.transform(z);elapsed+=time.monotonic()-t
        E=go.transfer.source(p.model,z)
        ref=go.transfer.readout(p.model,original,z,E);actual=go.transfer.readout(p.model,answer,z,E)
        errors=go.transfer.norms(p.model,actual-ref)/go.transfer.norms(p.model,ref)
        assert max(errors)<1e-10,(label,k,errors)
        checks.append(dict(case=label,contour_id=int(k),relative=dict(zip(go.task.FIELDS,errors.astype(float)))))
    pilots[label]=dict(setup_seconds=setup,solve_seconds=elapsed,solves=len(ids),dofs=p.model.size)
    nr=len(p.model.original.native);del p
L=go.laguerre(12);sample=np.ones((go.COUNT//2+1,128),dtype=np.clongdouble);times=[]
for _ in range(3):
    t=time.monotonic();a=go.coefficients(sample)
    for n in go.DEGREES:values=L[:,:n]@a[:n]
    coarse=L@go.coefficients(sample[::2])[:2048];times.append(time.monotonic()-t)
inversion=float(np.median(times)*2*nr/128)
for label,row in pilots.items():
    row['forecast_seconds']=1.2*(row['setup_seconds']+2048*row['solve_seconds']/row['solves']+inversion+20)
prior_cost=sum(json.loads((OUT/name).read_text())['seconds'] for name in ['pilot-budget.json','p4-result.json','spatial-pilot.json','spatial-result.json','operator-check.json','contrasts-pilot.json','profile-result.json'])
prior_cost+=json.loads((OUT/'profile-import-failure.json').read_text())['budget_charge_seconds']
seconds=time.monotonic()-started;remaining=1800-prior_cost-seconds
forecast=sum(row['forecast_seconds'] for label,row in pilots.items() if label!='baseline')
result=dict(classification='Counterexample candidate',transfer_equivalence_passed=True,checks=checks,cases=pilots,
    reconstruction128_seconds=times,reconstruction_forecast_seconds=inversion,forecast_seconds=forecast,
    seconds=seconds,remaining_compute_seconds=remaining,budget_passed=forecast<remaining,
    assumptions='Fixed frequency pilots for changed cases with20 percent margin. Native inversion uses fixed3-sample median. Historical import failure charged its whole60s cap. Remaining-budget alarm still required.')
task.write(OUT/'cached-pilot.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items() if k not in ['checks']},indent=2))
