import numpy as np,inspect
import apply_retained_native_acoustic as a
a.initialize_material(True);m=a.Material(128,128);k=16;p=m.point(k);field=np.zeros((5,m.n));q=np.zeros_like(p['Q']);q[1,10]=p['Q'][0,10]*a.prior.C**2*1e-8/a.prior.AMP
print('raw',m.raw.__qualname__,'deep_owner_same',m.deep.__self__ is m.model)
print([s for s in m.extended_sources[-1].splitlines() if any(x in s for x in ['sound','Z=','vf=','ps='])])
old=m.model.face_K.copy();deep=m.deep

def traced(*args):
 print('deep_K',float(np.max(abs(m.model.face_K/old-1))), 'bound_K',float(np.max(abs(deep.__self__.face_K/old-1))))
 return deep(*args)
m.deep=traced
z=m.raw(k,q,field,a.prior.AMP,p)[0];m.native_K[k]=old;w=m.raw(k,q,field,a.prior.AMP,p)[0]
print('delta per component',np.max(abs(z-w),axis=1).tolist())