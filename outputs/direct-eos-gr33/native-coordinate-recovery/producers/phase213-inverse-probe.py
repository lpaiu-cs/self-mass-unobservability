"""Check whether stored conserved coordinates identify native rates uniquely."""
from types import FunctionType
import json,resource,time
import numpy as np
import recover_remaining_joint_photons as r
out=r.OUT;result=out/'inverse-probe.json';assert not result.exists();start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));r.base.joint.previous.original.inf.incident.native.deadline(300)
FunctionType(r.base.base.prior.initialize.__code__,dict(r.base.base.prior.initialize.__globals__,OUT=out))()
m=r.base.base.owner.Model(128);z=dict(np.load(r.saved(128)));p=dict(np.load(out/'rejected-original-128.npz'));step=int(p['step']);g=p['gas'][1];t=p['stage_times'][1]
q=z['joint_stage_conserved_scaled'][2*step+1];old=z['joint_native_rates_scaled'][2*step+1]
qdiff=m.conserved(g)-q;cells=np.flatnonzero(np.any(qdiff!=0,axis=0));rows=[];witness=[]
for k in cells:
    values={0:g[k,0]}
    for direction in [-1,1]:
        value=g[k,0]
        for distance in range(1,9):
            value=np.nextafter(value,r.LD(direction)*r.LD('inf'));values[direction*distance]=value
    matches=[]
    for offset,value in sorted(values.items()):
        candidate=g.copy();candidate[k,0]=value;conserved=m.conserved(candidate)
        exact_cell=bool(np.array_equal(conserved[:,k],q[:,k]))
        if exact_cell:
            native=m.native(t,candidate)*m.units;matches.append((offset,candidate,native))
            rows.append(dict(cell=int(k),offset=offset,normalized_energy=str(value),exact_conserved_cell=True,native_cell_momentum=str(native[k,3]),native_max_difference=str(np.max(abs(native-old)))))
    if len(matches)>1:
        for a in range(len(matches)):
            for b in range(a+1,len(matches)):
                x,y=matches[a],matches[b]
                if not np.array_equal(x[2],y[2]):
                    assert np.array_equal(m.conserved(x[1]),m.conserved(y[1]))
                    witness.append(dict(cell=int(k),offsets=[x[0],y[0]],same_full_conserved=True,different_native=True,native_difference=str(np.max(abs(x[2]-y[2])))))
                    if not (out/'inverse-witness.npz').exists():np.savez_compressed(out/'inverse-witness.npz',time=t,gas_a=x[1],gas_b=y[1],conserved=m.conserved(x[1]),native_a=x[2],native_b=y[2])
                    break
            if witness:break
r.write(result,dict(classification='Counterexample candidate',changed_conserved_cells=cells.tolist(),preimages=rows,witnesses=witness,new_physical_steps=0,new_photon_solves=0,seconds=time.monotonic()-start,source_sha256=r.sha(__file__)));print(json.dumps(r.read(result)))
