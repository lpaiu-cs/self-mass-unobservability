from pathlib import Path
import inspect,json,resource,time
import numpy as np
import finish_stored_time_photons as r
out=r.OUT
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
r.base.joint.previous.original.inf.incident.native.deadline(180)
s=inspect.getsource(r.prior.old.original_equation)
a='        assert row[\'passed\'],row'
assert s.count(a)==1
s=s.replace(a,"        np.savez_compressed(OUT/'actual-defect.npz',defect=defect,rhs=rhs,solution=sol,initial=v)\n"+a)
ns=dict(r.prior.old.original_equation.__globals__);exec(compile(s,__file__,'exec'),ns)
r.prior.old.original_equation=ns['original_equation']
start=time.monotonic()
try:r.check()
except AssertionError as exc:print(str(exc))
z=np.load(out/'actual-defect.npz');m=r.base.base.owner.Model(64)
rows=[]
for vv in z['defect'].reshape(2,-1):
    ph,g=m.unpack(vv);rows.append(dict(photon_norm=str(np.linalg.norm(ph)),gas_norm=[str(np.linalg.norm(g[:,i])) for i in range(4)],gas_max=[str(np.max(abs(g[:,i]))) for i in range(4)],gas_cell=[int(np.argmax(abs(g[:,i]))) for i in range(4)]))
p=np.load(r.OLD/'rejected-original-64.npz');end=np.load('native-momentum-continuation210-work/last-accepted-64.npz')
g=p['gas'][-1];original=end['g'];active=m.material.active(p['stage_times'][-1])
diff=g-original;diff[~active]=0
report=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,defect_blocks=rows,
    restored_end_vs_actual_checkpoint_max=[str(np.max(abs(diff[:,i]))) for i in range(4)],
    restored_end_vs_actual_checkpoint_cells=[np.flatnonzero(diff[:,i]).tolist() for i in range(4)],
    actual_endpoint_conserved_equal=bool(np.array_equal(m.conserved(g)[:,active],m.conserved(original)[:,active])))
r.write(out/'actual-defect-location.json',report);print(json.dumps(report))
