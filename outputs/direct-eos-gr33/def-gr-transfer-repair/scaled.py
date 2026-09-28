"""Repair only the reference linear-solve scaling, preserving its operator."""
from pathlib import Path
import json
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import def_gr_transfer_repair as go

OUT=go.OUT
go.write(OUT/'reference-failure.json',dict(classification='Counterexample candidate',
    status='STOP_FULL_REFERENCE_RESIDUAL',stage='first frequency4+0i before any comparison',
    exception='AssertionError at full(): assert residual<1e-9',
    residual_value=None,meaning='The first attempt did not log the number. Do not infer it from a later attempt.',
    propagation_passed=False,original_failure_resolved=False))
go.write(OUT/'scaling-plan.json',dict(classification='Counterexample candidate',
    change='Same K+z^2M and source; use symmetric diagonal scaling of the complex matrix for its LU preconditioner. Retain extended original residual and at most5 corrections. Stop if the reference remains inadmissible; no evolution.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(go.__file__),OUT/'plan.json']}))
checks=[]


def full(K,M,rhs,z):
    matrix=(K+z*z*M).tocsc();scale=np.sqrt(abs(matrix.diagonal()));D=diags(1/scale)
    lu=splu((D@matrix@D).tocsc());x=(lu.solve(rhs/scale)/scale).astype(np.clongdouble)
    extended=matrix.astype(np.clongdouble);b=rhs.astype(np.clongdouble);errors=[]
    for _ in range(6):
        defect=b-extended@x;denom=abs(extended)@abs(x)+abs(b)+1e-100
        i=int(np.argmax(abs(defect)/denom));errors.append(float(abs(defect[i])/denom[i]))
        if errors[-1]<1e-18:break
        x+=(lu.solve(np.asarray(defect/scale,complex))/scale).astype(np.clongdouble)
    checks.append(dict(z=[z.real,z.imag],errors=errors,worst_row=i,rhs_abs=float(abs(b[i])),solution_abs=float(abs(x[i])),denominator=float(denom[i])))
    go.write(OUT/'reference-residuals.json',dict(classification='Counterexample candidate',checks=checks))
    assert errors[-1]<1e-9
    return x,errors[-1]


go.full=full;go.OUT=OUT/'scaled';go.run()
