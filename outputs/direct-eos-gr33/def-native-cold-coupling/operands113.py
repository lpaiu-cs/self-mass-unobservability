import sys,ctypes,json,time
import numpy as np
sys.path.insert(0,'verification')
import def_native_cold_coupling as task
n=task.prior.optical.ex.Native(cap=20);d=np.load(task.OUT/'diagnostic-states.npz');rows=[];start=time.monotonic()
for T in [240.,200.]:
    error=None
    try:n.state(float(d['x']),float(np.log(T)),float(d['y']))
    except Exception as exc:error=repr(exc)
    lib=n.ion.gas.gas_lib
    def arr(name,count):return np.ctypeslib.as_array((ctypes.c_double*count).in_dll(lib,'__mod_excitation_block_MOD_'+name)).copy()
    q=arr('qstar',290).reshape(29,10);qt=arr('qstart',290).reshape(29,10)
    rows.append(dict(T=T,error=error,qH2plus=q[-1].tolist(),qtH2plus=qt[-1].tolist(),qH=q[0].tolist(),qh2plus=float(arr('qh2plus',1)[0]),qh2=float(arr('qh2',1)[0]),c2t=float(arr('c2t',1)[0]),x=arr('x',5).tolist()))
    print(json.dumps(rows[-1]),flush=True)
task.write(task.OUT/'partition-operands.json',dict(classification='Counterexample candidate',rows=rows,native_calls=n.ion.calls,seconds=time.monotonic()-start))
