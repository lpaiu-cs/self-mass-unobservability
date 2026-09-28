"""Read saved178 states; no new evolution or extrapolated error certificate."""
from pathlib import Path
import time
import numpy as np
import return_native_stage_collisions as r
start=time.monotonic();rows=[]
for n in [64,128]:
    d=dict(np.load(r.OUT/f'material-{n}.npz'))
    p=dict(np.load(r.OUT/f'sweep-1/photons/pilot-{n}.npz'))
    old=dict(np.load(r.OLD/f'sweep-0/material/steps-{n}-reference-128.npz'))
    j=int(np.argmin(abs(old['t']-d['t'][-1])));assert abs(old['t'][j]-d['t'][-1])<1e-18
    z=d['delta_scaled'];last=old['history_scaled'][j]
    residual=(np.sum(abs(z[:2]-last[:2]),axis=1,dtype=r.LD)/np.maximum(np.sum(abs(z[:2]),axis=1,dtype=r.LD),1.)).astype(float).tolist()
    C,_=r.dense(p);lift=(z-C(d['t'][-1]))-d['lifted_history_scaled'][-1]
    lift_error=float(np.max(abs(lift))/max(np.max(abs(z)),1.))
    assert lift_error<1e-14
    rows.append(dict(steps=n,lagged_B_S_relative=residual,lift_reconstruction_relative=lift_error))
a=np.load(r.OUT/'material-64.npz')['delta_scaled'];b=np.load(r.OUT/'material-128.npz')['delta_scaled']
diff=abs(a[0]-b[0]);ids=np.argsort(diff)[-8:][::-1]
result=dict(classification='Counterexample candidate',rows=rows,new_evolution_steps=0,
    original_lagged_input_gate=.002,lagged_B_S_accepted=bool(max(v for row in rows for v in row['lagged_B_S_relative'])<.002),
    baryon_time_top8=[dict(cell=int(j),absolute_scaled=float(diff[j]),fraction=float(diff[j]/np.sum(diff,dtype=r.LD))) for j in ids],
    scope='Saved prefix state comparison only. Do not count178 neutral-transport acceptance as full block or final charge acceptance.',
    seconds=time.monotonic()-start,source_sha256=r.sha(__file__))
r.write(r.OUT/'saved-state-check.json',result);print(result)
