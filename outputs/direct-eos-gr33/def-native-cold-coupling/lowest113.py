import sys,traceback,json,time
import numpy as np
sys.path.insert(0,'verification')
import def_native_cold_coupling as t
n=t.native_variant('scaled',16);d=np.load(t.OUT/'diagnostic-states.npz');start=time.monotonic();failure=None
try:n.state(float(d['x']),np.log(80.),float(d['y']))
except Exception as e:failure=repr(e);traceback.print_exc()
t.write(t.OUT/'lowest-probe.json',dict(classification='Counterexample candidate',failure=failure,native_calls=n.ion.calls+n.variant_initial_calls,seconds=time.monotonic()-start))
