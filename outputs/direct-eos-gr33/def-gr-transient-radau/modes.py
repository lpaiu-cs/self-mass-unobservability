"""Bounded actual-pencil modal diagnosis; no new evolution or EOS calls."""
import json
import signal
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigs,splu
import def_gr_transient_radau as work

out=work.OUT;task=work.prior.old.task;signal.alarm(90);begin=time.monotonic()
ray=dict(np.load(task.coupled.OUT/'fine-rays.npz'));rad=task.Radiation(ray,work.prior.OUT/'fine-bank.npz',False)
bg=task.coupled.Background(rad,2);fn,_=task.coupled.reactive.symbolic();K,D,base=task.coupled.assemble(bg,fn)
forcing=base(1.)+K@rad.heat.lift(1.,bg.nodes)
scale=np.asarray(abs(K).sum(1)).ravel();static=splu(K.multiply((1/scale)[:,None]).tocsc()).solve(forcing/scale)
rows=[];arrays={};N=len(static)
for band in [32,128,512]:
    shift=-float(band)**2;matrix=(K-shift*D).tocsc();scale=np.asarray(abs(matrix).sum(1)).ravel()
    lu=splu(matrix.multiply((1/scale)[:,None]).tocsc())
    forward=lambda x:lu.solve((D@x)/scale)
    transpose=lambda x:D.T@(lu.solve(x,trans='T')/scale)
    op=LinearOperator((N,N),matvec=forward,dtype=float);adjoint=LinearOperator((N,N),matvec=transpose,dtype=float)
    values,right=eigs(op,k=8,which='LM',tol=1e-10,maxiter=300,v0=static/np.linalg.norm(static))
    leftvalues,left=eigs(adjoint,k=8,which='LM',tol=1e-10,maxiter=300,v0=np.ones(N)/np.sqrt(N))
    for j,mu in enumerate(values):
        lam=shift+1/mu;r=right[:,j];l=left[:,np.argmin(abs(leftvalues-mu))]
        match=float(np.min(abs(leftvalues-mu))/abs(mu));den=l@r
        coefficient=(l@static)/den
        residual=np.linalg.norm(K@r-lam*(D@r))/max(np.linalg.norm(K@r)+abs(lam)*np.linalg.norm(D@r),1e-100)
        rows.append(dict(band=band,lambda_real=float(lam.real),lambda_imag=float(lam.imag),omega_abs=float(np.sqrt(abs(lam))),pencil_residual=float(residual),left_match=match,
            coefficient_real=float(coefficient.real),coefficient_imag=float(coefficient.imag),left_right_overlap=float(abs(den))))
    arrays[str(band)+'_values']=shift+1/values;arrays[str(band)+'_right']=right;arrays[str(band)+'_left']=left;arrays[str(band)+'_leftvalues']=leftvalues
    print('BAND',band,'seconds',time.monotonic()-begin,'omega',[round(r['omega_abs'],3) for r in rows[-8:]],'maxres',max(r['pencil_residual'] for r in rows[-8:]),flush=True)
np.savez_compressed(out/'modes.npz',static=static,**arrays)
work.write(out/'modes.json',dict(classification='Counterexample candidate',seconds=time.monotonic()-begin,rows=rows,complete_spectral_certificate=False))
signal.alarm(0)
