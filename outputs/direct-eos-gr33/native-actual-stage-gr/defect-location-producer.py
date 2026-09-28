from pathlib import Path
import json,numpy as np
w=Path('native-integer-continuation229-work');z=np.load(w/'rejected-joint-stage.npz');p=np.load(w/'last-accepted-64.npz')
n=p['g'].shape[0];nx=p['x'].size;dim=nx+4*n
d=z['defect'].reshape(2,dim);sol=z['solution'].reshape(2,dim)
norm=np.linalg.norm(d)/json.loads((w/'rejected-joint-stage.json').read_text())['equations'][-1]['relative']
rows=[]
for stage in range(2):
    gas=d[stage,nx:].reshape(n,4);components=[float(np.linalg.norm(gas[:,c])/norm) for c in range(4)]
    top=[]
    for k in np.argsort(abs(d[stage]))[-10:][::-1]:
        top.append(dict(index=int(k),gas_cell=None if k<nx else int((k-nx)//4),component=None if k<nx else int((k-nx)%4),relative=float(d[stage,k]/norm),solution=float(sol[stage,k]),ULP=float(np.spacing(sol[stage,k]))))
    rows.append(dict(stage=stage,photon=float(np.linalg.norm(d[stage,:nx])/norm),gas_components=components,top=top))
result=dict(classification='Counterexample candidate',shape=list(p['x'].shape),normalizer=float(norm),rows=rows)
Path('native-integer-continuation229-work/defect-location.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
