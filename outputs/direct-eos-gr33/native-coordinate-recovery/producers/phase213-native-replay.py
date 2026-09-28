"""Locate the exact archival-native difference without another photon solve."""
from pathlib import Path
from types import FunctionType
import json,resource,time
import numpy as np
import recover_remaining_joint_photons as r

out=r.OUT;result=out/'native-replay.json';assert not result.exists();start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));r.base.joint.previous.original.inf.incident.native.deadline(300)
FunctionType(r.base.base.prior.initialize.__code__,dict(r.base.base.prior.initialize.__globals__,OUT=out))()
m=r.base.base.owner.Model(128);z=dict(np.load(r.saved(128)));p=dict(np.load(out/'rejected-original-128.npz'));step=int(p['step'])
rows=[]
for j,(t,g) in enumerate(zip(p['stage_times'],p['gas'])):
    old=z['joint_native_rates_scaled'][2*step+j];first=m.native(t,g)*m.units
    m.jacobian(t,g);after=m.native(t,g)*m.units
    diff=after-old;ids=np.argwhere(diff!=0);q=z['joint_stage_conserved_scaled'][2*step+j];roundtrip=m.conserved(g)-q
    rows.append(dict(stage=j,native_before_J_equals_after=bool(np.array_equal(first,after)),different_entries=len(ids),
        relative_L1=(np.sum(abs(diff),axis=0)/np.maximum(np.sum(abs(old),axis=0),r.LD('1e-290'))).astype(float).tolist(),
        differences=[dict(cell=int(k),component=int(c),saved=str(old[k,c]),replayed=str(after[k,c]),difference=str(diff[k,c]),saved_ulp=str(np.spacing(old[k,c]))) for k,c in ids[:30]],
        conserved_roundtrip=(np.sum(abs(roundtrip),axis=1)/np.maximum(np.sum(abs(q),axis=1),r.LD('1e-290'))).astype(float).tolist(),
        roundtrip_max_abs_by_component=[str(v) for v in np.max(abs(roundtrip),axis=1)]))
r.write(result,dict(classification='Counterexample candidate',rows=rows,seconds=time.monotonic()-start,source_sha256=r.sha(__file__),new_physical_steps=0,new_photon_solves=0));print(json.dumps(rows))
