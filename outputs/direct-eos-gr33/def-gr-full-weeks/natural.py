"""Use the existing banded node ordering for sparse LU; compare saved pilot."""
from pathlib import Path
import json
import time
import signal
import numpy as np
import def_gr_full_weeks as go

OUT=go.OUT
assert not (OUT/'natural-pilot.json').exists()
go.write(OUT/'natural-plan.json',dict(classification='Counterexample candidate',
    reason='Initial577.59s forecast exceeds450s first-case cap. Test native sparse LU NATURAL column ordering on all16 saved pilot points; preserve pivoting and original extended residual.',
    decision='Only use if every four-readout transfer changes by less1e-8 and the revised measured forecast fits450s. No contour or degree reduction; no time evolution in this test.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(go.__file__),OUT/'pilot.npz',OUT/'pilot-budget.json']}))
lu=go.splu;go.splu=lambda A:lu(A,permc_spec='NATURAL')
signal.alarm(90);start=time.monotonic();p=go.Problem();setup=time.monotonic()-start
pilot=np.load(OUT/'pilot.npz');checks=[];times=[]
for k,old in zip(pilot['ids'],pilot['answers']):
    z,_=go.contour(k,12);start=time.monotonic();new=p.transform(z);times.append(time.monotonic()-start)
    E=go.transfer.source(p.model,z)
    a=go.transfer.readout(p.model,old,z,E);b=go.transfer.readout(p.model,new,z,E)
    errors=go.transfer.norms(p.model,b-a)/go.transfer.norms(p.model,a)
    checks.append(dict(contour_id=int(k),relative=dict(zip(go.task.FIELDS,errors.astype(float)))))
old=json.loads((OUT/'pilot-budget.json').read_text());forecast=1.4*(setup+sum(times)/len(times)*2048+old['inversion_forecast_seconds']+20)
row=dict(classification='Counterexample candidate',setup_seconds=setup,solve_seconds=sum(times),first_case_forecast_seconds=forecast,
    transfer_equivalence_passed=all(max(r['relative'].values())<1e-8 for r in checks),checks=checks,linear_residual=p.error)
go.write(OUT/'natural-pilot.json',row);signal.alarm(0);print('NATURAL PILOT',json.dumps(row),flush=True)
