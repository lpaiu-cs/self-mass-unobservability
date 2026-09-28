import time,json,numpy as np
import apply_retained_native_acoustic as a
from retained_deep_collision_precision import Pressure
start=time.monotonic();a.native.deadline(60);a.initialize_photon();p=Pressure(a.motion.State());a.initialize_material(False);m=a.Material(128,128)
d=np.load(a.MATERIAL/'steps-128-reference-128.npz');rows=[]
for k in [0,8,16]:
 z=d['history_scaled'][int(np.argmin(abs(d['t']-m.t[k])))];field=m.fields(m.t[k])[2]
 _,delta,error=p(m,k,z,field);relative=float(np.sum(abs(error))/max(np.sum(abs(delta)),1.))
 rows.append(dict(**p.rows[-1],arithmetic_over_increment=relative));print(json.dumps(rows[-1]),flush=True)
 np.savez_compressed(a.OUT/f'pressure-pilot-{k}.npz',**p.last)
assert max(r['arithmetic_over_increment'] for r in rows)<.002
upper=2*max(r['seconds'] for r in rows)*34+25
result=dict(classification='Counterexample candidate',passed=True,rows=rows,upper_seconds=upper,eligible=upper<160,seconds=time.monotonic()-start)
a.write(a.OUT/'pressure-pilot.json',result);assert result['eligible']
