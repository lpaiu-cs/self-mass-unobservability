"""Isolate the original per-interval native cache lifecycle."""
from types import FunctionType
import json,resource,time
import numpy as np
import resume_coordinate_exact_photons as r
out=r.OUT;result=out/'cache-probe.json';assert not result.exists();start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));r.base.joint.previous.original.inf.incident.native.deadline(300)
FunctionType(r.base.base.prior.initialize.__code__,dict(r.base.base.prior.initialize.__globals__,OUT=out))()
m=r.base.base.owner.Model(64);z=dict(np.load(r.prior.saved(64)));rows=[]
def measure(i,label):
    t=z['joint_stage_times'][i];g=r.restored_gas(m,z['joint_stage_conserved_scaled'][i]);old=z['joint_native_rates_scaled'][i];native=m.native(t,g)*m.units
    rows.append(dict(stage=i,label=label,exact=bool(np.array_equal(native,old)),relative_L1=(np.sum(abs(native-old),axis=0)/np.maximum(np.sum(abs(old),axis=0),r.LD('1e-290'))).astype(float).tolist()))
measure(0,'first_interval');measure(16,'continuous_model');m.material.thermo_jets.clear();measure(16,'clear_thermo_jets_only');measure(17,'same_new_cache');measure(18,'same_new_cache')
r.write(result,dict(classification='Counterexample candidate',rows=rows,material_dicts={k:len(v) for k,v in vars(m.material).items() if isinstance(v,dict)},new_physical_steps=0,new_photon_solves=0,seconds=time.monotonic()-start,source_sha256=r.sha(__file__)));print(json.dumps(r.read(result)))
