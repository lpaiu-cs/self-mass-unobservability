"""Read-only original spatial-L1 source comparison; preserve206failure."""
import json
from pathlib import Path
import numpy as np
import return_resolved_joint_history as r

out=r.OUT;path=out/'original-source-norm.json';assert not path.exists()
a=dict(np.load(out/'gr/source-64.npz'));b=dict(np.load(out/'gr/source-128.npz'))
ids=np.array([np.argmin(abs(b['t']-t)) for t in a['t']]);assert np.max(abs(a['t']-b['t'][ids]))<1e-18
keys=list(r.read(out/'sources.json')['source_time'])
old_relative=r.run.owner.joint.previous.run.c.relative
errors=dict(zip(keys,old_relative(np.stack([a[k] for k in keys],axis=1),np.stack([b[k][ids] for k in keys],axis=1))))
row=dict(classification='Counterexample candidate',passed=max(errors.values())<.02,
    source_time_spatial_L1=errors,original_gate=.02,
    correction='206source used a maximum-cell norm. The original186source norm sums absolute errors over cells and takes the maximum over times. Recompute that unchanged original norm from the same saved source arrays; both norms fail. No integration or input change.',
    source_knots=[len(a['t']),len(b['t'])],horizon_seconds=float(a['t'][-1]),
    old_failure_preserved=True,GR_input_applied=False,final_charge_conclusion='unadjudicated',
    source_sha256=r.sha(__file__),bindings={str(out/'gr'/f'source-{n}.npz'):r.sha(out/'gr'/f'source-{n}.npz') for n in [64,128]})
r.write(path,row);print(json.dumps(row))
