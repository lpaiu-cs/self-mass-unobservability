import inspect,numpy as np
import apply_retained_native_acoustic as a
a.initialize_photon();s=a.motion.State();s.select(64);m=s.m;model=m.model;b=model.bulk;e=b.eos
for fn in [model.set_material,model.kinetic,model.energy,b.eos.gas,b.eos.radiation]:
 print(inspect.getsource(fn))
print('EOS',type(e),sorted(vars(e)))
print('BASE',type(e.base),sorted(vars(e.base)))
print('MECH',[(k,np.shape(v)) for k,v in vars(model.mech).items() if k in ['x','xi','div','K']])
print('MODELMASS',model.cx,model.m.a0)
np.savez_compressed(a.OUT/'deep-highprecision-input.npz',mass0=model.mass0,V=b.volume,a=b.d['a'],cx=model.cx,a0=model.m.a0,xi=model.mech.xi)
