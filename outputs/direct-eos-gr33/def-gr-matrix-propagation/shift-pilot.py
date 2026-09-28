from pathlib import Path
import json, time, numpy as np
from scipy.sparse.linalg import splu
from scipy.linalg import eig
import def_gr_matrix_propagation as b
o=b.OUT
plan=dict(classification='Counterexample candidate',claim='Project the bounded shifted resolvent rather than K^-1 D: avoid subtracting the very large static inverse when resolving fast modes.',retained='Same finite position constraints, source and128-vector saved kinetic basis. No eigenvalue clipping or symmetrization.',budget=dict(hard_seconds=30,new_evolution_paths=0),source=b.prior.prior.digest(Path(__file__)),operator_source=b.prior.prior.digest(Path(b.__file__)))
b.write(o/'shift-plan.json',plan)
start=time.monotonic()
p=b.prior.Problem(b.prior.prior.OUT/'fine-bank.npz');L,_=b.source_map(p);m=b.Propagator(p,L);f=b.Finite(m)
d=dict(np.load(o/'finite-pilot.npz'));Q,W=d['Q'],d['W'];sigma=256.**2
A=(m.K-sigma*m.D).tocsc();scale=np.asarray(abs(A).sum(1)).ravel();lu=splu(A.multiply((1/scale)[:,None]).tocsc())
RQ=np.column_stack([f.restrict(sigma*lu.solve((m.D@f.extend(q))/scale)) for q in Q.T])
H=W.T@np.column_stack([f.kinetic(q) for q in RQ.T])
values=eig(H,right=False);lam=sigma*(1+values)/values;growth=np.sqrt(lam+0j)
np.savez_compressed(o/'shift-pilot.npz',H=H,eigenvalues=values,lambda_values=lam)
result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,projected_positive_lambda=int(np.sum(lam.real>0)),projected_complex_modes=int(np.sum(abs(lam.imag)>1e-10*abs(lam))),max_growth=float(growth.real.max()),max_frequency=float(abs(growth.imag).max()),nonsymmetry=float(np.linalg.norm(H-H.T)/np.linalg.norm(H)),original_failure_resolved=False)
b.write(o/'shift-pilot.json',result);print(json.dumps(result),flush=True)
