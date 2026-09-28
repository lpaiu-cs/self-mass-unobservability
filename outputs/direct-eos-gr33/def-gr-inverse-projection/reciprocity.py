"""Locate the failed inverse reciprocity at fixed64 columns; no time paths."""
from pathlib import Path
from types import SimpleNamespace
import time
import signal
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
from scipy.linalg import cholesky_banded,cho_solve_banded
import def_gr_multishift as method

task=method.task;OUT=method.OUT.parent
task.write(OUT/'reciprocity-plan.json',dict(classification='Counterexample candidate',
    claim='Separate matrix assembly reciprocity from shifted linear solution reciprocity after the field-inverse pilot failed symmetry. Fixed64 trial columns, no actual new evolution or EOS calls.',
    methods=['original LU','explicit symmetric assembly LU','explicit symmetric band Cholesky'],
    gates=dict(inverse_symmetry=1e-12,linear_residual=1e-9),budget=dict(hard_seconds=60,CPU_threads=1),
    bindings={str(p):task.digest(p) for p in [Path(__file__),Path(method.__file__),OUT/'field-inverse/plan.json']}))
signal.alarm(60);start=time.monotonic();m=method.space.Model(4,task.BANK/'fine-bank.npz')
_,_,Q,scales,_=method.basis(m,64);P=task.Projection(m);white=P.C.T@Q;T=P.C@white
results={};matrices={}
for mode in ['original','symmetric-LU','symmetric-Cholesky']:
    K=m.K if mode=='original' else (m.K+m.K.T)/2
    M=m.M if mode=='original' else (m.M+m.M.T)/2
    A=(K+P.sigma*M).tocsc();scale=np.sqrt(A.diagonal());D=diags(1/scale);matrix=(D@A@D).tocsc()
    skew=A-A.T;coo=matrix.tocoo();width=int(max(abs(coo.row-coo.col)))
    meta=dict(absolute_matrix_skew=float(np.max(abs(skew.data),initial=0)),scaled_matrix_skew=float(np.max(abs((matrix-matrix.T).data),initial=0)))
    if mode.endswith('Cholesky'):
        band=np.zeros((width+1,m.size))
        for j in range(width+1):band[j,:m.size-j]=matrix.diagonal(-j)
        try:factor=cholesky_banded(band,lower=True)
        except np.linalg.LinAlgError as exc:
            results[mode]=dict(**meta,factor_failed=str(exc));continue
        solve=lambda rhs:cho_solve_banded((factor,True),rhs/scale[:,None])/scale[:,None]
    else:
        factor=splu(matrix);solve=lambda rhs:factor.solve(rhs/scale[:,None])/scale[:,None]
    S=solve(T);extended=A.astype(np.longdouble)
    for _ in range(2):S+=solve(np.asarray(T.astype(np.longdouble)-extended@S.astype(np.longdouble),float))
    R=P.sigma*T.T@S
    error=float(np.max(abs(T.astype(np.longdouble)-extended@S.astype(np.longdouble))/(abs(extended)@abs(S)+abs(T)+1e-100)))
    reciprocity=float(np.max(abs(R-R.T))/np.max(abs(R)))
    results[mode]=dict(**meta,linear_residual=error,inverse_symmetry=reciprocity)
    matrices[mode]=R
np.savez_compressed(OUT/'reciprocity-matrices.npz',**matrices)
task.write(OUT/'reciprocity.json',dict(classification='Counterexample candidate',results=results,seconds=time.monotonic()-start,
    original_failure_resolved=False,new_time_paths=0,new_EOS_calls=0))
signal.alarm(0);print('RECIPROCITY',results,flush=True)
