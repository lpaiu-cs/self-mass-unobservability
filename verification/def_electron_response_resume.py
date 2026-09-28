"""Resume the unchanged collision model after a NumPy-bool JSON failure."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
import def_electron_correlated_response as model


def main():
    h=model.h;out=model.OUT
    assert not (out/'result.json').exists()
    # The original source and original preregistration remain byte-for-byte.
    h.write(out/'resume-plan.json',dict(classification='Counterexample candidate',
        failure='prepare reached the angular-control serialization and raised TypeError: Object of type bool is not JSON serializable (numpy.bool_). No full-cohort production had run.',
        change='Convert diagnostic scalars to builtin float/bool. Re-run the same angular controls, then the unchanged production function, gates, grid sizes and 60s hard budget.',
        sources={str(p.relative_to(h.ROOT)):h.digest(p) for p in [Path(__file__),Path(model.__file__),out/'plan.json',out/'pilot.json']}))
    errors=[]
    for u in [.001,.1,.999,1.,10.,1e6]:
        for w in [1e-8,.01,.1,1.,20.,1e5]:
            for v in [0.,.2,.9]:
                top=np.log1p(1/u)
                f=lambda t:.5*(-np.expm1(-t))*(1-v*u*np.expm1(t))*(-np.expm1(-w*u*np.expm1(t)))
                exact=quad(f,0,top,epsabs=1e-28,epsrel=2e-11,limit=200)[0]
                actual=model.angular(np.array([u]),np.array([w]),np.array([v]))[0]
                errors.append(float(abs(actual/exact-1)))
    h.write(out/'angular-controls.json',dict(classification='Counterexample candidate',max_relative_error=max(errors),passed=bool(max(errors)<1e-6)))
    print('ANGULAR',max(errors),flush=True)
    model.run()


if __name__=='__main__':main()
