"""Numerical full-matrix inertia check for the inversion abscissa."""
from pathlib import Path
import signal
import time
import numpy as np
from scipy.sparse import diags
from scipy.linalg import cholesky_banded
import def_gr_full_weeks as go

OUT=go.OUT;assert not (OUT/'abscissa-check.json').exists();signal.alarm(60);start=time.monotonic()
go.write(OUT/'abscissa-plan.json',dict(classification='Counterexample candidate',
    claim='Check the full finite matrix K+12^2*M has a positive numerical Cholesky factor, rather than inferring the contour side from projected eigenvalues.',
    limitation='Floating factorization check only, not an interval certificate of positive definiteness or a continuum stability theorem.',
    budget_seconds=60,bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(go.__file__),go.task.BANK/'fine-bank.npz']}))
p=go.Problem();A=(p.K+144*p.M).tocsc();s=1/np.sqrt(A.diagonal());D=diags(s);scaled=(D@A@D).tocsc()
coo=scaled.tocoo();w=int(max(abs(coo.row-coo.col)));n=A.shape[0];band=np.zeros((w+1,n))
for j in range(w+1):band[j,:n-j]=scaled.diagonal(-j)
C=cholesky_banded(band,lower=True);L=diags([C[j,:n-j] for j in range(w+1)],-np.arange(w+1),shape=(n,n)).tocsc()
residual=L@L.T-scaled;error=float(np.max(abs(residual.data),initial=0));assert error<1e-12 and C[0].min()>0
row=dict(classification='Counterexample candidate',full_shifted_cholesky_passed=True,shift=144,dofs=n,bandwidth=w,
    minimum_factor_diagonal=float(C[0].min()),scaled_factor_residual=error,seconds=time.monotonic()-start,
    scope='Numerical evidence that sigma12 lies right of the real unstable poles of this finite symmetric pair; not outward-rounded or continuum certification.')
go.write(OUT/'abscissa-check.json',row);signal.alarm(0);print(row)
