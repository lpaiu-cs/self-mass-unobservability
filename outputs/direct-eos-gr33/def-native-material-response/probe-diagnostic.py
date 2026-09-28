import sys,time,signal,json
from pathlib import Path
import numpy as np
import def_native_material_response as x
signal.signal(signal.SIGALRM,x.flow.old.optical.timeout);signal.alarm(20)
start=time.monotonic();m=x.Material(128);d=np.load(x.OUT/'steps-64-reference-128.npz');t=m.t[8];idx=np.argmin(abs(d['t']-t));z=d['history_scaled'][idx];rates=[];rows=[]
for p in [4.,2.,1.,.5,.25]:
    m.min_probe=np.inf;r,l,dt=m.rhs(t,z,p);rates.append(r)
    rows.append(dict(probe=p,epsilon=m.min_probe,L1=np.sum(abs(r),axis=1).tolist()))
for j in range(4):
    a,b=rates[j:j+2];delta=abs(a-b);ids=np.argmax(delta,axis=1)
    rows[j]['half_difference']=(np.sum(delta,axis=1)/np.sum(abs(b),axis=1)).tolist()
    rows[j]['max_cell']=ids.tolist();rows[j]['differences_at_max']=delta[np.arange(4),ids].tolist()
    rows[j]['rate_at_max']=b[np.arange(4),ids].tolist()
    rows[j]['radius_at_max']=m.R[ids].tolist()
p=m.point(8);zero=np.zeros((4,m.n));a=m.raw(8,zero,np.zeros((5,m.n)),0.)
rows.append(dict(baseline_repeat_flux=np.max(abs(a[0]-p['flux']),axis=1).tolist(),baseline_repeat_gravity=float(np.max(abs(a[1]-p['gravity'])))))
out=dict(classification='Counterexample candidate',time=float(t),rows=rows,seconds=time.monotonic()-start)
x.write(x.OUT/'probe-diagnostic.json',out);print(json.dumps(out))
