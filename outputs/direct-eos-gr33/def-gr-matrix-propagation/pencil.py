import json,time,signal
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigs,splu
import def_gr_matrix_propagation as b
o=b.OUT
b.write(o/'pencil-plan.json',dict(classification='Counterexample candidate',claim='Check whether the fast growing projected pair is an actual eigenpair of the unchanged full K,D. No evolution or altered operator.',budget=dict(bands=2,modes_per_band=4,hard_seconds=45),source=b.prior.prior.digest(Path(__file__)),operator_source=b.prior.prior.digest(Path(b.__file__))))
signal.alarm(45);start=time.monotonic()
p=b.prior.Problem(b.prior.prior.OUT/'fine-bank.npz');L,_=b.source_map(p);m=b.Propagator(p,L);f=b.Finite(m)
l=np.load(o/'shift-pilot.npz')['lambda_values'];ids=np.argsort(np.sqrt(l+0j).real)
candidate=l[ids[-1]]
bands=[candidate,complex(candidate.real,0)]
rows=[]
for sigma in bands:
 A=(m.K-sigma*m.D).tocsc();scale=np.asarray(abs(A).sum(1)).ravel();lu=splu(A.multiply((1/scale)[:,None]).tocsc())
 def op(q):return f.restrict(lu.solve((m.D@f.extend(np.asarray(q,complex)))/scale))
 R=LinearOperator((f.size,f.size),matvec=op,dtype=complex)
 v,V=eigs(R,k=4,which='LM',tol=1e-11,maxiter=500)
 for mu,q in zip(v,V.T):
  lam=sigma+1/mu;full=lu.solve((m.D@f.extend(q))/scale)/mu
  residual=m.K@full-lam*(m.D@full)
  scale_res=abs(m.K)@abs(full)+abs(lam)*(abs(m.D)@abs(full))
  row=dict(shift=[sigma.real,sigma.imag],lambda_value=[lam.real,lam.imag],growth=[float(np.sqrt(lam).real),float(np.sqrt(lam).imag)],normwise_residual=float(np.linalg.norm(residual)/np.linalg.norm(scale_res)),component_residual=float(np.max(abs(residual)/(scale_res+1e-100))))
  rows.append(row)
result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,rows=rows,original_failure_resolved=False)
b.write(o/'pencil-result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)
