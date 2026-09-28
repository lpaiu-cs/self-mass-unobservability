"""Exact nonadditivity witness; no numerical evolution or charge inference."""
import time
started=time.monotonic()
import numpy as np
import complete_native_fluid_time as run
f=run.joint.previous.engine.face.minmod_direction
zero=np.array([0.]);one=np.array([1.])
values=[f(zero,zero,one,zero)[0],f(zero,zero,zero,one)[0],f(zero,zero,one,one)[0]]
assert values==[0.,0.,1.]
run.write(run.OUT/'superposition-boundary.json',dict(classification='Proven',values=values,
    scope='Exact nonadditivity of the original minmod directional map at zero background slopes. This does not quantify the actual full-driver superposition error.',
    seconds=time.monotonic()-started,source_sha256=run.sha(__file__),native_source_sha256=run.sha(run.joint.previous.engine.face.__file__)))
