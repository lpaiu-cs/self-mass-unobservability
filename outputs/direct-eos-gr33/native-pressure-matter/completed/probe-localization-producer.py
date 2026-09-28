from pathlib import Path
import json,resource,time
import numpy as np
import return_native_pressure_matter as run

out=run.OUT;start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.original.inf.incident.native.deadline(45)
failed=out/'sweep-1/material/pilot-64.npz';row=run.read(failed.with_suffix('.json'))
assert not row['passed'] and not (out/'probe-localization.json').exists()
files=[Path(__file__),Path(run.__file__),failed,failed.with_suffix('.json'),out/'pilot-receipt.json']
run.write(out/'probe-localization-plan.json',dict(classification='Counterexample candidate',
    question='Localize the already observed half/nominal/double RHS discrepancy on the saved failed material state. Distinguish analytic deep faces from the shared/native atmosphere; no new probe range or evolution.',
    budget_seconds=45,new_evolution_steps=0,existing_probe_factors=[.5,1.,2.],
    bindings={str(p):run.sha(p) for p in files}))
run.initialize();m=run.c.Material(128,64);d=np.load(failed);t=float(d['time']);z=d['delta_scaled'];values=[];probe=[]
for p in [.5,1.,2.]:
    m.min_probe=np.inf;v,l,dt=m.rhs(t,z,p);values.append(v)
    probe.append(dict(factor=p,epsilon=float(m.min_probe),cfl=float(dt)))
v=np.asarray(values);rows=[]
for k,name in enumerate(['B','S','E','H']):
    norm=max(np.sum(abs(v[1,k])),1.)
    for i in [0,2]:
        delta=abs(v[i,k]-v[1,k]);order=np.argsort(delta)[-8:][::-1]
        rows.append(dict(component=name,other_factor=[.5,1.,2.][i],relative=float(sum(delta)/norm),
            deep_relative=float(sum(delta[:m.nb-1])/norm),shared_relative=float(sum(delta[m.nb-1:m.nb+1])/norm),
            outer_relative=float(sum(delta[m.nb+1:])/norm),
            cells=[dict(cell=int(j),nominal=float(v[1,k,j]),other=float(v[i,k,j]),relative=float(delta[j]/norm),radius_E=float(m.rE[j])) for j in order]))
np.savez_compressed(out/'probe-localization-rates.npz',t=t,state=z,factors=[.5,1.,2.],rates=v)
run.write(out/'probe-localization.json',dict(classification='Counterexample candidate',rows=rows,probes=probe,
    deep_cells=m.nb,physical_branch_ratio=float(m.physical_branch_ratio),seconds=time.monotonic()-start,
    evolved_new_state=False,full_goal_complete=False))
print(json.dumps(dict(rows=rows[:4],probes=probe,seconds=time.monotonic()-start)))
