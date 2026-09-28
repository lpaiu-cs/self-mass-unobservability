from pathlib import Path
import json
import numpy as np

p=Path('native-balanced-krylov196-work')
z=dict(np.load(p/'failed-linear-64.npz'))
rows=[]
for key in ['rhs','solution','residual']:
    v=z[key].reshape(2,-1);g=v[:,-4*531:].reshape(2,531,4)
    rows.append(dict(name=key,norm=float(np.linalg.norm(v)),
        photons=float(np.linalg.norm(v[:,:-4*531])),
        gas=np.sqrt(np.sum(g*g,axis=(0,1))).astype(float).tolist(),
        gas_max=np.max(abs(g),axis=(0,1)).astype(float).tolist()))
sol=z['solution'].reshape(2,-1)[:,-4*531:].reshape(2,531,4)
res=z['residual'].reshape(2,-1)[:,-4*531:].reshape(2,531,4)
ids=np.argsort(abs(res[:,:,2]).ravel())[-12:]
out=dict(classification='Counterexample candidate',rows=rows,
    large_B=[dict(stage=int(i//531),cell=int(i%531),x=float(sol[:,:,2].ravel()[i]),
                 residual=float(res[:,:,2].ravel()[i]),ulp=float(abs(np.spacing(sol[:,:,2].ravel()[i])))) for i in ids],
    calls=[{k:c[k] for k in ['info','seconds','iterations','true_relative','weight_min']} for c in json.loads((p/'linear-64.json').read_text())['calls']])
print(json.dumps(out,indent=2))
