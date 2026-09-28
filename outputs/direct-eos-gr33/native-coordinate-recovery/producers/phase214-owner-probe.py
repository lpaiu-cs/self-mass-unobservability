"""Replay the original stage call sequence at the first remaining mismatch."""
from types import FunctionType
import json,resource,time
import numpy as np
import resume_coordinate_exact_photons as r
out=r.OUT;result=out/'owner-probe.json';assert not result.exists();start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));r.base.joint.previous.original.inf.incident.native.deadline(300)
FunctionType(r.base.base.prior.initialize.__code__,dict(r.base.base.prior.initialize.__globals__,OUT=out))()
m=r.base.base.owner.Model(64);z=dict(np.load(r.prior.saved(64)));i=16;t=z['joint_stage_times'][i];q=z['joint_stage_conserved_scaled'][i];g=r.restored_gas(m,q);old=z['joint_native_rates_scaled'][i];rows=[]
def measure(label):
    native=m.native(t,g)*m.units
    rows.append(dict(label=label,exact=bool(np.array_equal(native,old)),relative_L1=(np.sum(abs(native-old),axis=0)/np.maximum(np.sum(abs(old),axis=0),r.LD('1e-290'))).astype(float).tolist()))
measure('fresh');m.local(t);m.source(t);measure('after_local_source');m.jacobian(t,g);measure('after_J')
for k in range(16):
    now=z['joint_stage_times'][k];gg=r.restored_gas(m,z['joint_stage_conserved_scaled'][k]);m.local(now);m.source(now);m.native(now,gg);m.pressure(now,gg)
measure('after_prior_calls')
r.write(result,dict(classification='Counterexample candidate',rows=rows,new_physical_steps=0,new_photon_solves=0,seconds=time.monotonic()-start,source_sha256=r.sha(__file__)));print(json.dumps(rows))
