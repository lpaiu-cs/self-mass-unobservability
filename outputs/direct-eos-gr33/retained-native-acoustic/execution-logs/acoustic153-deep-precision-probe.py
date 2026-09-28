import numpy as np,time,json
import apply_retained_native_acoustic as a
from retained_deep_collision_precision import difference
start=time.monotonic();a.native.deadline(90);a.initialize_photon();s=a.motion.State();s.select(64);k=2;m=s.m
s.precision(True);z=m.motion[k];c=s.coefficients(k,np.zeros_like(z),0.)
d,center,rounds,root=difference(s,k,z,1.,40)
err=[float(np.max(abs(v-center[j]))/np.max(abs(v))) for j,v in enumerate([c['emit'][:m.nb],c['loss'][:m.nb]])]
np.savez_compressed(a.OUT/'deep-precision-40-probe.npz',delta=d,center=center,rounding=rounds)
print(json.dumps(dict(seconds=time.monotonic()-start,owner=err,root=root)),flush=True)
assert max(err)<1e-9
