"""Bounded saved-state check of the heat lift at the remaining spatial failure."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_gr_full_weeks as go

OUT=go.OUT;B=OUT/'beta1024';signal.alarm(60);started=time.monotonic()
assert not (OUT/'lift-check.json').exists()
go.write(OUT/'lift-check-plan.json',dict(classification='Counterexample candidate',
    claim='Test whether the saved new-boundary spatial discrepancy is dominated by the directly subtracted degree-dependent heat momentum lift.',
    decision='Evaluate the unchanged lift at the frozen maximum-discrepancy time on p2/p4. No new trajectory, grid, fitting or acceptance change. A small direct lift difference rules out only this readout-subtraction mechanism, not forcing error.',
    budget=dict(hard_seconds=60,CPU_threads=1,new_EOS_calls=0,new_evolutions=0),
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),OUT/'spatial-location.json',B/'p2-2048.npz',B/'p4-2048.npz',Path(go.space.__file__)]}))
location=json.loads((OUT/'spatial-location.json').read_text());t=location['tau'];values=[];details=[]
for degree in [2,4]:
    m=go.space.Model(degree,go.task.BANK/'fine-bank.npz');r=m.original.native;w=m.original.weights;mask=m.original.masks[1]
    flux,energy=m.heat.faces(t);N,a=m.original.radiation.geometry.metric(r);speed=a/N*r*go.task.fem.base.task.h.gr.C
    values.append(speed*(m.nativeV[0]@(m.H@(flux*m.original.radiation.geometry.tc))))
    faces=m.heat.edges;active=m.heat.face_ids
    details.append(dict(degree=degree,dofs=m.size,source_faces=len(active),first_source_face=int(active.min()),
        last_source_face=int(active.max()),last_source_radius_R=float(faces[active.max()]),
        readout_peak_lift_m_s=float(abs(values[-1][mask]).max())))
    del m
d2=np.load(B/'p2-2048.npz');d4=np.load(B/'p4-2048.npz');j=round(t*64)
dv=d4['native_velocity'][j]-d2['native_velocity'][j];dl=values[1]-values[0]
norm=lambda x:float(np.sqrt(np.sum(w[mask]*x[mask]**2)/w[mask].sum()))
result=dict(classification='Counterexample candidate',tau=t,details=details,
    field_difference_rms_m_s=norm(dv),lift_difference_rms_m_s=norm(dl),
    lift_to_field_difference_ratio=norm(dl)/norm(dv),
    unlifted_response_difference_rms_m_s=norm(dv+dl),seconds=time.monotonic()-started,
    original_failure_resolved=False,new_evolutions=0,
    scope='Only direct output subtraction tested. Spatial dynamics and weak source assembly remain to be resolved.')
go.write(OUT/'lift-check.json',result);signal.alarm(0);print(json.dumps(result,indent=2))
