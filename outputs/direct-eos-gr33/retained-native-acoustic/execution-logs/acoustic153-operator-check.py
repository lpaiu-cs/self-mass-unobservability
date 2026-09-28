import numpy as np,json
import apply_retained_native_acoustic as a
from pathlib import Path
a.initialize_material(True);m=a.Material(128,128);k=16;p=m.point(k);z=np.zeros_like(p['Q']);q=z.copy();q[1,10]=p['Q'][0,10]*a.prior.C**2*1e-8/a.prior.AMP;f=np.zeros((5,m.n))
v0=m.raw(k,z,f,0,p)[0];v1=m.raw(k,q,f,a.prior.AMP,p)[0];K=m.native_K[k].copy();m.native_K[k]=m.model.face_K.copy()
w0=m.raw(k,z,f,0,p)[0];w1=m.raw(k,q,f,a.prior.AMP,p)[0];m.native_K[k]=K
rows=[]
for n in [64,128]:
 old=np.load(a.native.OUT/f'response/direct-material/steps-{n}-reference-128.npz');new=np.load(a.DIRECT/f'steps-{n}-reference-128.npz')
 rows.append(dict(n=n,exact=np.array_equal(old['history_scaled'],new['history_scaled']),delta=float(np.max(abs(old['history_scaled']-new['history_scaled']))),deep=float(np.max(abs(new['history_scaled'][:,:,:19])))))
print(json.dumps(dict(Kchange=float(np.max(abs(K/m.model.face_K-1))),probe=float(np.sum(abs((v1-v0)-(w1-w0)))/np.sum(abs(v1-v0))),paths=rows)))