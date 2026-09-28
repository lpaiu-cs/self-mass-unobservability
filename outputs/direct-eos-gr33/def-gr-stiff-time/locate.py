"""Locate the actual remaining velocity differences; no evolution."""
import sys
from pathlib import Path
import json
import numpy as np
sys.path.insert(0,'verification')
import def_gr_stiff_time as run

out=Path(__file__).resolve().parent
model=run.space.Model(4,run.task.BANK/'fine-bank.npz')
r=model.original.native;ids=np.clip(np.searchsorted(model.cells,r)-1,0,len(model.cells)-2)
K=model.K;M=model.M;omega=np.sqrt(np.maximum(K.diagonal()/M.diagonal(),0))
checks={}
for name,directory in [('stiff',out),('gauss',run.common.OUT/'gauss'),('midpoint',run.task.OUT/'direct-time')]:
    a=np.load(directory/'p4-1024.npz');b=np.load(directory/'p4-2048.npz')
    va=a['native_velocity'];vb=b['native_velocity'];weights=b['weights']
    rows=[]
    for mask in [np.ones(len(r),bool),*b['masks']]:
        delta=va-vb;contribution=weights[None,:]*delta**2*mask[None,:]
        j=int(np.argmax(contribution.sum(1)));order=np.argsort(contribution[j])[::-1]
        total=contribution[j].sum();leading=[]
        for i in order[:12]:
            cell=int(ids[i]);ix=model.indices[cell*4:(cell+1)*4+1,0];ix=ix[ix>=0]
            leading.append(dict(native=int(i),r=float(r[i]),cell=cell,cell_width=float(np.diff(model.cells)[cell]),
                fraction=float(contribution[j,i]/total),velocity=float(vb[j,i]),difference=float(delta[j,i]),
                diagonal_frequency_min=float(min(omega[ix])),diagonal_frequency_max=float(max(omega[ix]))))
        rows.append(dict(tau=j/64,weighted_difference=float(np.sqrt(total/weights[mask].sum())),
            leading12_fraction=float(contribution[j,order[:12]].sum()/total),leading=leading))
    checks[name]=rows
result=dict(classification='Counterexample candidate',meaning='Signed native-vector differences, not RMS-history differences or a true-error bound. Diagonal frequencies are local indicators, not eigenfrequencies.',checks=checks)
run.write(out/'location.json',result)
print(json.dumps(result,indent=2))
