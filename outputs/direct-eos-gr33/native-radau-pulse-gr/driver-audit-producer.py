"""Independent live-driver check of the polynomial actually sent to GR."""
from pathlib import Path
from types import FunctionType
import json,resource,time
import numpy as np
from scipy.interpolate import PPoly
import apply_driver_aware_radau_gr as r

out=r.OUT;start=time.monotonic();assert not (out/'driver-polynomial-fixed-audit.json').exists();assert (out/'amplitude-fix.json').exists()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));r.prior.prior.evolution.joint.previous.original.inf.incident.native.deadline(300)
FunctionType(r.prior.prior.initialize.__code__,dict(r.prior.prior.initialize.__globals__,OUT=out))()
model=r.base.run.owner.Model(64);rows=[]
for n in [64,128]:
    d=dict(np.load(out/'gr'/f'source-{n}.npz'));knots,co=r.coefficients(d);times=(knots[:-1]+knots[1:])/2
    phi=np.array([model.driver.wave(float(t),model.driver.xc)[0]/(model.driver.centers*r.AMP) for t in times])
    errors={}
    for k in r.KEYS:
        actual=PPoly(np.asarray(co[k][::-1],float),knots)(times)
        expected=r.prior.poly(d['t'],d['state_coeff_'+k])(times)
        if k not in r.KEYS[-2:]:
            a,b=d['geometry_coeff_'+k];expected=expected+(a[None]+times[:,None]*b)*phi
        errors[k]=float(np.max(np.sum(abs(actual-expected),axis=-1) if actual.ndim>1 else abs(actual-expected))/max(np.max(np.sum(abs(expected),axis=-1) if expected.ndim>1 else abs(expected)),1e-290))
    rows.append(dict(clock=n,interval_midpoints=len(times),errors=errors))
result=dict(classification='Counterexample candidate',passed=max(max(row['errors'].values()) for row in rows)<1e-12,rows=rows,seconds=time.monotonic()-start,
    source_sha256=r.sha(__file__),bindings={str(out/'gr'/f'source-{n}.npz'):r.sha(out/'gr'/f'source-{n}.npz') for n in [64,128]},
    scope='Actual consumer polynomial versus the original live Driver at every support-segment midpoint; no physical-state or final-charge acceptance.')
r.write(out/'driver-polynomial-fixed-audit.json',result);print(json.dumps(result));assert result['passed'],result
