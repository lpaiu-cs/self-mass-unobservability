"""Replay the archived rate with its original scalar time arithmetic."""
from types import FunctionType
import json,resource,time
import numpy as np
import resume_coordinate_exact_photons as r
out=r.OUT;result=out/'time-probe.json';assert not result.exists();start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));r.base.joint.previous.original.inf.incident.native.deadline(300)
FunctionType(r.base.base.prior.initialize.__code__,dict(r.base.base.prior.initialize.__globals__,OUT=out))()
m=r.base.base.owner.Model(64);z=dict(np.load(r.prior.saved(64)));rows=[]
for i in [0,1,14,15,16,17,18,100,221]:
    t=z['joint_stage_times'][i];q=z['joint_stage_conserved_scaled'][i];g=r.restored_gas(m,q);old=z['joint_native_rates_scaled'][i]
    v=[]
    for label,now in [('stored',t),('float64',np.float64(t))]:
        value=m.native(now,g)*m.units
        v.append(dict(label=label,time_identical=bool(now==t),exact=bool(np.array_equal(value,old)),relative_L1=(np.sum(abs(value-old),axis=0)/np.maximum(np.sum(abs(old),axis=0),r.LD('1e-290'))).astype(float).tolist(),max_difference=str(np.max(abs(value-old)))))
    rows.append(dict(stage=i,stored_type=str(t.dtype),background_type=str(m.t.dtype),clock=v))
r.write(result,dict(classification='Counterexample candidate',rows=rows,new_physical_steps=0,new_photon_solves=0,seconds=time.monotonic()-start,source_sha256=r.sha(__file__)));print(json.dumps(rows))
