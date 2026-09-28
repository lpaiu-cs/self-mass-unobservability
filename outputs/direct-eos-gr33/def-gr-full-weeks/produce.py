"""Execute the measured native-order full resolvent, unchanged first cap."""
from pathlib import Path
from types import SimpleNamespace
import json
import signal
import resource
import numpy as np
from scipy.sparse import csc_matrix
import def_gr_full_weeks as go

OUT=go.OUT;pilot=json.loads((OUT/'natural-pilot.json').read_text())
assert pilot['transfer_equivalence_passed'] and pilot['first_case_forecast_seconds']<450
assert not (OUT/'execution-plan.json').exists()
go.write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',
    change='Use measured NATURAL column ordering only. Same original pivoting,linear equations,extended residual,contour,time basis,input and gates. Reuse original pilot resolvents whose equivalence was verified.',
    forecast_seconds=pilot['first_case_forecast_seconds'],hard_seconds=450,
    decision='One full degree4 contour. Stop on any original propagation/contour gate. Do not enlarge counts on failure.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(go.__file__),OUT/'plan.json',OUT/'natural-plan.json',OUT/'natural-pilot.json',OUT/'pilot.npz']}))
for plan in ['plan.json','natural-plan.json','execution-plan.json']:
    for p,h in json.loads((OUT/plan).read_text())['bindings'].items():assert go.task.digest(Path(p))==h,p
lu=go.splu;go.splu=lambda A:lu(A,permc_spec='NATURAL')
# Actual extended transform on a two-field problem with a weak forced field.
p=go.Problem.__new__(go.Problem);K=np.array([[9.,1e-11],[1e-11,25.]])
p.K=csc_matrix(K);p.M=csc_matrix(np.eye(2));p.Kx=p.K.astype(np.clongdouble);p.Mx=p.M.astype(np.clongdouble)
p.load=csc_matrix([[1.],[2e-12]]).astype(np.clongdouble);p.error=0.
heat=SimpleNamespace(rates=np.array([[1e12]]),amplitude=np.array([[1.]]),edges=np.array([0.]),face_ids=np.array([0]))
p.model=SimpleNamespace(heat=heat,original=SimpleNamespace(radiation=SimpleNamespace(geometry=SimpleNamespace(tc=1.))))
errors=[]
for z in [4+0j,4+233j,12+1001j]:
    expected=np.linalg.solve(K+z*z*np.eye(2),np.array([1.,2e-12])*1e12/(z*z*(z+1e12)))
    actual=p.transform(np.clongdouble(z));errors.append(float(np.max(abs(actual-expected)/abs(expected))))
assert max(errors)<1e-12
# Check actual FFT, long-precision Laguerre recurrence and endpoint kernel.
half=np.zeros((go.COUNT//2+1,2),dtype=np.clongdouble)
for k in range(1,len(half)):
    z,f=go.contour(k,12);half[k]=f/(z*(z*z+np.array([9.,640000.])))
half[-1]=half[-1].real;L=go.laguerre(12);coeff=go.coefficients(half);series=L@coeff[:2048]
exact=(1-np.cos(np.array([3.,800.])[None,:]*go.TIMES[:,None]))/np.array([9.,640000.])
inversion=float(np.max(abs(series[1:]-exact[1:])*np.array([9.,640000.])))
assert inversion<1e-6
ell=np.zeros(go.COUNT,dtype=np.longdouble);ell[:2048]=L[-1];kernel=go.fft(ell)[:go.COUNT//2+1]/go.COUNT
endpoint=np.sum((kernel[:,None]*half).real*np.r_[1.,np.full(len(half)-2,2.),1.][:,None],axis=0)
endpoint_error=float(np.max(abs(endpoint-series[-1])));assert endpoint_error<1e-14
go.write(OUT/'control-extended.json',dict(classification='Counterexample candidate',passed=True,weak_field_transform_relative=errors,
    oscillator_scaled_error=inversion,endpoint_kernel_absolute=endpoint_error,
    scope='Actual extended transform and reconstruction controls, not proof of GR response accuracy.'))
signal.alarm(450);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
go.solve(reuse=True);signal.alarm(0)
