"""Bind the stopped recovery checkpoints and a local array-guard correction."""
from pathlib import Path
import ast,json,shutil
import numpy as np
import resume_coordinate_exact_photons as r
p=r.OUT;assert not (p/'resume-plan.json').exists()
old=p/'array-guard-failed-producer.py';current=Path(r.__file__)
functions=lambda path:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef)}
a,b=functions(old),functions(current)
for name in ['restored_gas','interval_starts','interval_model']:assert a[name]==b[name]
rows=[]
for n in [64,128]:
    x=dict(np.load(p/f'accepted-{n}.npz'));y=dict(np.load(p/f'resume-input-{n}.npz'))
    for key in x:assert np.array_equal(x[key],y[key]),key
    assert int(y['step'])==len(json.loads(str(y['logs'])))
    rows.append(dict(clock=n,checkpoint_steps=int(y['step']),checkpoint_value_identity=True))
first=p/'first-dispatch';first.mkdir()
for name in ['progress-64.json','progress-128.json','expanded-recovery-64.py','expanded-recovery-128.py']:
    shutil.copyfile(p/name,first/name)
r.write(p/'resume-plan.json',dict(classification='Conjectural',
    correction='Compare the array photon normalization with array_equal at original model transitions. The failed scalar truth test changed no equations. Commit accepted endpoints immediately as well as before steps, so a subsequent reconstruction error cannot discard the last accepted photon solve.',
    same_inverse_and_model_functions=True,all652precheck_reused=True,rows=rows,
    reuse='Restore frozen7coarse/11fine photon checkpoints. The eighth accepted coarse photon endpoint was not checkpointed before the failed array guard; repeat only that one conditional photon block. All native-fluid steps remain reused. Reuse no failed state as accepted.',
    gates=r.read(p/'plan.json')['gates'],coarse_cap_seconds=7200,fine_cap_seconds=10800,CPU_threads_per_path=1,virtual_GiB_per_path=6,
    original_failure_preserved=True,producer_sha256=r.sha(current),
    inputs={str(q):r.sha(q) for q in [p/'resume-input-64.npz',p/'resume-input-128.npz',old,p/'check_interval-receipt.json',p/'restart-check.json',p/'coarse-receipt.json']},
    source_sha256=r.sha(__file__),final_charge_conclusion='unadjudicated',full_goal_complete=False))
