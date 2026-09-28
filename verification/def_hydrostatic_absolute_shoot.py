"""Recover the nonzero solve with absolute-coordinate shooting differences.

The relative finite differences in MINPACK vanish for the near-zero fitted
log-radius/log-mass offsets. Reuse the accepted GR solve and EOS table.
"""
from pathlib import Path
from types import FunctionType
import json
import numpy as np
from scipy.optimize import root as scipy_root
import def_hydrostatic_background as h

OUT=h.OUT/'absolute-shoot'


def root(function,x,method,options):
    steps=np.array([2e-6,2e-6,2e-6,.001,.0001])
    def jacobian(value):
        columns=[]
        for i,step in enumerate(steps):
            direction=np.zeros(5);direction[i]=step
            columns.append((function(value+direction)-function(value-direction))/(2*step))
        return np.asarray(columns).T
    return scipy_root(function,x,method=method,jac=jacobian,options=dict(xtol=options['xtol'],maxfev=45))


def run():
    assert not OUT.exists();OUT.mkdir()
    accepted=json.loads((h.OUT/'result-0.json').read_text())
    assert accepted['GR_reproduction_relative']<1e-8
    failed=json.loads((h.OUT/'progress-0.001.json').read_text())['history']
    differences=np.array([r['residual'] for r in failed])-failed[0]['residual']
    assert np.all(differences[4:6]==0)
    h.write(h.OUT/'shoot-failure.json',dict(classification='Counterexample candidate',passed=False,
        reason='Relative MINPACK differences in near-zero log radius/mass were 1.377e-17 and 9.526e-18; both produced identically zero residual columns, followed by overflow in the unconstrained proposal.',
        recovery='Absolute central differences, same equations and gates; reuse accepted GR background and table. No old failed path overwritten.'))
    plan=json.loads((h.OUT/'plan.json').read_text());plan['phi_infinity']=[.001]
    plan['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),h.OUT/'result-0.json',h.OUT/'extended-table.npz',h.OUT/'shoot-failure.json']})
    plan['numerical_correction']='Explicit absolute central shooting differences [2e-6,2e-6,2e-6,.001,.0001]; no physical or acceptance change.'
    h.write(OUT/'plan.json',plan)
    pilot=json.loads((h.OUT/'shoot-pilot.json').read_text());pilot['parameters']=accepted['parameters']
    h.write(OUT/'shoot-pilot.json',pilot)
    runner=FunctionType(h.run.__code__,dict(h.run.__globals__,OUT=OUT,root=root))
    runner()


if __name__=='__main__':run()
