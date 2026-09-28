import json,signal,time
import numpy as np
import def_native_feedback_return as run
start=time.monotonic();signal.signal(signal.SIGALRM,run.old.base.flow.old.optical.timeout);signal.alarm(30)
m=run.Material(64,128);d=np.load(run.OUT/'steps-128-reference-64.npz');z=d['delta_scaled'];t=float(d['time'])
rates=[];probes=[.125,.25,.5,1.,2.,4.,8.,16.]
for p in probes:rates.append(m.rhs(t,z,p)[0])
rows=[]
for i in range(len(probes)-1):
    diff=abs(rates[i]-rates[i+1]);den=np.maximum(np.sum(abs(rates[i+1]),axis=1),1.)
    rows.append(dict(probes=probes[i:i+2],relative=(diff.sum(1)/den).tolist(),momentum_cells=np.argsort(diff[1])[-8:].tolist(),momentum_error_fractions=(np.sort(diff[1])[-8:]/den[1]).tolist()))
r=dict(classification='Counterexample candidate',rows=rows,seconds=time.monotonic()-start,physical_branch_ratio=m.physical_branch_ratio)
run.write(run.OUT/'probe.json',r);print(json.dumps(r));signal.alarm(0)
