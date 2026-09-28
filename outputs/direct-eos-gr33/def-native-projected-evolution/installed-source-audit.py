"""Bounded read-only check of the actual installed initial source."""
import json
import time
import numpy as np
import def_native_projected_evolution as task

start=time.monotonic();m=task.Coupled();b=m.bulk;f=m.flow;n=b.n
ref=np.load(task.INPUT/'balanced-initial-state.npz');p,u,*_=b.eos.gas(np.zeros(n),np.zeros(n))
rho,v,lt,y=f.primitive(f.initial);f.eos.y=y;pa,ua,*_=f.eos(rho,lt)
density=np.r_[b.d['rho'],rho*f.eos.rho0].astype(task.LD)
pressure=np.r_[p,pa*f.eos.rho0*task.C**2].astype(task.LD)
specific=np.r_[u,ua*task.C**2].astype(task.LD)
energy=density*(task.LD(m.cx)*task.LD(task.C)**2+specific)
mask=density>0
errors=dict(density=float(np.max(abs(density[mask]/ref['density'][mask]-1))),
    pressure=float(np.max(abs(pressure[mask]/ref['gas_pressure'][mask]-1))),
    energy_including_rest=float(np.max(abs(energy[mask]/ref['gas_energy'][mask]-1))))
result=dict(classification='Counterexample candidate',errors=errors,seconds=time.monotonic()-start,
    interpretation='Actual installed material source compared with the accepted local-Gamma initial constraints. This is an input consistency diagnostic, not an added response acceptance gate or uniform native EOS certificate.',
    new_fluid_steps=0,bindings={str(p):task.sha(p) for p in [__file__,task.__file__,task.INPUT/'balanced-initial-state.npz']})
task.write(task.OUT/'installed-source-audit.json',result);print(json.dumps(result))
