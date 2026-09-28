"""Check extended solution iterates, not merely extended residual evaluation."""
from pathlib import Path
import time
import signal
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import def_gr_multishift as method

task=method.task;OUT=method.OUT.parent
task.write(OUT/'reciprocity-refinement-plan.json',dict(classification='Counterexample candidate',
    claim='Determine whether storing the inverse-solve iterate in double precision, rather than extended precision, limits reciprocity despite tiny backward residuals.',
    methods=['original coefficients, extended iterates','symmetric coefficients, extended iterates'],
    fixed='Same64 columns and input;8 refinement corrections maximum; no time evolution.',
    gates=dict(inverse_symmetry=1e-12),budget=dict(hard_seconds=90,CPU_threads=1,new_EOS_calls=0),
    bindings={str(p):task.digest(p) for p in [Path(__file__),OUT/'reciprocity.json']}))
signal.alarm(90);start=time.monotonic();m=method.space.Model(4,task.BANK/'fine-bank.npz')
_,_,Q,scales,_=method.basis(m,64);P=task.Projection(m);T=P.C@(P.C.T@Q);rhs=T.astype(np.longdouble)
results={};matrices={}
for name in ['original','symmetric']:
    K=m.K.astype(np.longdouble);M=m.M.astype(np.longdouble)
    if name=='symmetric':K=(K+K.T)/2;M=(M+M.T)/2
    A=K+P.sigma*M;scale=np.sqrt(np.asarray(A.diagonal(),float));D=diags(1/scale)
    lu=splu((D@A.astype(float)@D).tocsc())
    X=(lu.solve(T/scale[:,None])/scale[:,None]).astype(np.longdouble);rows=[]
    for j in range(9):
        if j:
            defect=rhs-A@X
            X+=(lu.solve(np.asarray(defect/scale[:,None],float))/scale[:,None]).astype(np.longdouble)
        if j in [0,2,4,8]:
            hi=np.asarray(X,float);lo=np.asarray(X-hi,float)
            R=P.sigma*(T.T@hi+T.T@lo)
            rows.append(dict(corrections=j,inverse_symmetry=float(np.max(abs(R-R.T))/np.max(abs(R))),
                componentwise_residual=float(np.max(abs(rhs-A@X)/(abs(A)@abs(X)+abs(rhs)+1e-100)))))
    results[name]=rows;matrices[name]=R
np.savez_compressed(OUT/'refined-reciprocity-matrices.npz',**matrices)
result=dict(classification='Counterexample candidate',results=results,seconds=time.monotonic()-start,
    original_failure_resolved=False,new_time_paths=0,new_EOS_calls=0)
task.write(OUT/'reciprocity-refinement.json',result);signal.alarm(0);print('REFINED RECIPROCITY',results,flush=True)
