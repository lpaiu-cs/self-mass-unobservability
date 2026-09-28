"""Reuse saved trajectories to locate the actual spatial velocity discrepancy."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_gr_energy_modes as modes
import def_gr_canonical_regular as regular

OUT=Path(__file__).resolve().parent;assert not (OUT/'inspection.json').exists()
modes.write(OUT/'inspection-plan.json',dict(classification='Counterexample candidate',checkpoint='95d99987',
    claim='Locate spatial velocity differences using saved fine/coarse trajectories and compare cell widths with actual acoustic propagation distance. No new evolution or EOS call.',
    budget=dict(hard_seconds=60,new_evolution_paths=0,new_EOS_calls=0),source=modes.digest(Path(__file__))))
signal.alarm(60);start=time.monotonic();modes.fem.canonical=regular
model=modes.fem.Model(modes.evolution.BANK/'fine-bank.npz');prior=modes.fem.OUT/'regular'
fine=dict(np.load(prior/'fine-512.npz'));coarse=dict(np.load(prior/'spatial-512.npz'))
r=model.original.native;dv=fine['native_velocity']-coarse['native_velocity'];v=fine['native_velocity'];w=fine['weights']
point=model.bg.sample(r);N,a=model.original.radiation.geometry.metric(r)
sound=N/a*np.sqrt(point['gamma']*point['p']/(point['e']+point['p']))
R=model.original.radiation.geometry.R;ids=np.clip(np.searchsorted(model.grid,r)-1,0,len(model.grid)-2);dx=np.diff(model.grid)[ids]
rows=[]
for name,mask in [('global',np.ones(len(r),bool)),('old',fine['masks'][0]),('new',fine['masks'][1])]:
    indices=np.flatnonzero(mask);error=abs(dv[:,mask]);j,k=np.unravel_index(np.argmax(error),error.shape);i=int(indices[k])
    contribution=np.sum(w[mask][None,:]*dv[:,mask]**2,axis=0);rank=np.argsort(contribution)[::-1]
    n90=int(np.searchsorted(np.cumsum(contribution[rank]),.9*contribution.sum())+1)
    rows.append(dict(region=name,peak_tau=j/64,native_cell=i,native_original_cell=len(r)-1-i,
        radius_m=float(r[i]*R),cell_width_m=float(dx[i]*R),acoustic_distance_full_horizon_m=float(sound[i]*R),
        acoustic_cells=float(sound[i]/dx[i]),fine_velocity=float(v[j,i]),coarse_velocity=float(coarse['native_velocity'][j,i]),
        max_point_difference=float(error[j,k]),cells_for90percent_error=n90,window_cells=len(indices),
        acoustic_cells_quantiles=np.quantile(sound[mask]/dx[mask],[0,.1,.5,.9,1]).tolist()))
flux,energy=model.heat.faces(1.);lift_velocity=model.H@(flux*model.original.radiation.geometry.tc)
lift_ratio=float(np.linalg.norm(lift_velocity)/np.linalg.norm(fine['qt']))
result=dict(classification='Counterexample candidate',regions=rows,heat_lift_velocity_ratio=lift_ratio,
    seconds=time.monotonic()-start,original_failure_resolved=False)
modes.write(OUT/'inspection.json',result);signal.alarm(0);print(json.dumps(result),flush=True)
