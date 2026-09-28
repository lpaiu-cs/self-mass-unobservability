"""Saved-only same-solution audit and canonical incident fields for GR reuse."""
import json,resource,time
from pathlib import Path
import numpy as np
import couple_full_incident_fluid as run

start=time.monotonic();cpu=time.process_time();out=run.OUT
run.joint.previous.original.inf.incident.native.deadline(15)
assert run.read(out/'pilot-receipt.json')['seconds']+15<180
p=dict(np.load(out/'sweep-1/photons/pilot-64.npz'));ld=np.longdouble
w=p['joint_stage_weights'].astype(ld)
assert np.array_equal(w,p['accepted_angular_quadrature_weights'])
assert np.array_equal(p['joint_stage_times'],p['accepted_angular_times'])
q=p['conserved_material_history'][-1]/run.AMP
# Eref-Etilde is archived explicitly, so no EOS or normalization re-fit is needed.
et=(q[2]-p['energy_offset_reference'][-1]/run.AMP)
final=np.column_stack([et,q[3],q[0],q[1]])
expected=np.sum(w[:,None,None]*(p['joint_native_rates_scaled']+p['joint_collision_rates_scaled']),axis=0,dtype=ld)
got=final+p['material_floor_discard_scaled']
err=np.sum(abs(got-expected),axis=0,dtype=ld)/np.maximum(np.sum(abs(got)+abs(expected),axis=0,dtype=ld),ld('1e-290'))
assert max(err)<1e-8,err
port=run.joint.previous.run.packets(out/'sweep-1/photons/pilot-64.npz')[2]
driver=run.drive.Driver(8);tt=np.unique(np.r_[p['t'],p['joint_stage_times']]);rows=[driver.at(t) for t in tt]
data={key:np.array([r[key] for r in rows]) for key in rows[0]}
assert np.any(data['delta_log_lapse'])
np.savez_compressed(out/'actual-incident-fields.npz',t=tt,radius_E=p['radius_E'],**data)
run.write(out/'saved-audit.json',dict(classification='Counterexample candidate',passed=True,
    local_material_balance=err.astype(float).tolist(),angular_port_relative=port,
    actual_field_count=len(tt),actual_stage_count=len(w),horizon_seconds=float(p['t'][-1]),
    final_charge_conclusion='unadjudicated',self_GR_return_closed=False,
    seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
    peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
    source_sha256=run.sha(__file__),producer_sha256=run.sha(run.__file__),
    solution_sha256=run.sha(out/'sweep-1/photons/pilot-64.npz'),fields_sha256=run.sha(out/'actual-incident-fields.npz')))
