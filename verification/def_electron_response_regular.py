"""Resolve the narrow angular boundary layer without enlarging the grid."""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.laguerre import laggauss
from scipy.integrate import quad
import def_electron_correlated_response as model

h=model.h
original_angular=model.angular
OUT=model.base.OUT/'correlated-regular'


def angular(u,w,v2):
    u,w,v2=np.broadcast_arrays(u,w,v2);result=np.empty(u.shape)
    narrow=(u*w>=100)&(w>50)
    if np.any(~narrow):result[~narrow]=original_angular(u[~narrow],w[~narrow],v2[~narrow])
    if np.any(narrow):
        a,b,v=u[narrow],w[narrow],v2[narrow];sw=a*b
        nodes,weights=laggauss(16)
        # Subtract the exponential boundary layer using y=w*q. The omitted
        # y>w tail is exponentially small here; independent controls remain.
        correction=np.sum(weights*nodes*(1-v[:,None]*nodes/b[:,None])/(1+nodes/sw[:,None])**2,axis=1)/(2*sw*sw)
        result[narrow]=model.base.angular(a,v)-correction
    assert np.min(result)>0
    return result


def main():
    assert not OUT.exists();OUT.mkdir()
    previous=model.OUT;plan=json.loads((previous/'plan.json').read_text())
    paths=[Path(__file__),previous/'angular-controls.json',previous/'resume-plan.json',Path(model.__file__)]
    plan['sources'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths})
    plan['angular_repair']='At u=0.001,w=1e5 the prior whole-interval rule missed a narrow exponential layer (1.7481e-6 relative control error). Integrate that layer in y=w*q with 16 Laguerre nodes; retain the original rule elsewhere. No change to physical equations, thresholds, cohort or total energy-node counts.'
    h.write(OUT/'plan.json',plan)
    h.write(OUT/'pilot.json',json.loads((previous/'pilot.json').read_text()))
    model.OUT=OUT;model.angular=angular
    errors=[]
    for u in [.001,.1,.999,1.,10.,1e6]:
        for w in [1e-8,.01,.1,1.,20.,1e5]:
            for v in [0.,.2,.9]:
                top=np.log1p(1/u)
                f=lambda t:.5*(-np.expm1(-t))*(1-v*u*np.expm1(t))*(-np.expm1(-w*u*np.expm1(t)))
                reference=quad(f,0,top,epsabs=1e-28,epsrel=2e-11,limit=200)[0]
                value=angular(np.array([u]),np.array([w]),np.array([v]))[0]
                errors.append(float(abs(value/reference-1)))
    h.write(OUT/'angular-controls.json',dict(classification='Counterexample candidate',max_relative_error=max(errors),passed=max(errors)<1e-6))
    print('ANGULAR',max(errors),flush=True)
    model.run()


if __name__=='__main__':main()
