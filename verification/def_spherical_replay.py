"""Replay saved native inputs and all coupled finite equations without EOS calls."""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import json
import numpy as np
import def_spherical_balanced_wave as balanced
import def_spherical_coupled as original
from gr_implicit_coupled_evolution import CachedOnly

e,m,ld=original.e,original.m,original.ld


def initialize_saved(_):
    data=dict(np.load(original.imported.initial.OUT/'initial.npz'))
    star=original.imported.initial.original.make_star(CachedOnly(),data)
    fn=m.prior.initialize
    return FunctionType(fn.__code__,dict(fn.__globals__,parent=SimpleNamespace(initialize=lambda _:star)))(CachedOnly())


def verify():
    plan=json.loads((balanced.OUT/'plan.json').read_text());p=plan['original_plan']
    for rel,sha in plan['bindings'].items():assert e.digest(e.ROOT/rel)==sha,rel
    result=json.loads((balanced.OUT/'result.json').read_text())
    assert result['passed']
    original.imported.initialize=initialize_saved
    checks=[]
    for path,amplitude in [(original.OUT/'path-0/endpoint.npz',ld(0)),(balanced.OUT/'endpoint.npz',ld('1e-6'))]:
        saved=dict(np.load(path));star=original.initialize(None,amplitude,ld(p['h_seconds']))
        star.__class__=balanced.BalancedStar
        zero=np.zeros_like(star.base);star.amplitude=ld(0);z0=star.evaluate(zero);star.amplitude=amplitude
        delta=saved['delta'];y=star.base+delta;logA=-2*saved['psi']**2
        rows=zip(y[:,0]-3*logA,y[:,1]-logA,y[:,5:])
        star.material_cache.update({e.material_key(row):aux for row,aux in zip(rows,saved['raw'])})
        star.field_guess=saved['psi'],saved['Pi'],saved['Phi']
        value,z=original.residual(star,delta,(zero,z0),(zero,z0),star.h,m.weights(star.h,None))
        norm=float(np.max(abs(value)/m.ATOL));assert norm<=1,norm
        assert np.array_equal(z['raw'],saved['raw'])
        differences={}
        for key in ['rho','T','P','u','rest','m','mf','a','N','Q','aux','dU','psi','Pi','Phi','Ephi']:
            scale=max(np.max(abs(saved[key])),ld('1e-100'))
            difference=float(np.max(abs(z[key]-saved[key]))/scale)
            assert difference<2e-16,(key,difference)
            differences[key]=difference
        assert z['scalar_relative_balance']<p['scalar_relative_budget_tolerance']
        assert z['scalar_residual']<p['scalar_residual_tolerance']
        checks.append(dict(classification='Counterexample candidate',path=path.relative_to(e.ROOT).as_posix(),
            passed=True,native_residual_norm=norm,scalar_residual=z['scalar_residual'],
            scalar_relative_balance=z['scalar_relative_balance'],maximum_array_relative_difference=differences,
            native_arguments_exact_cache_hit=True,new_native_EOS_calls=0))
    report=dict(classification='Counterexample candidate',passed=True,paths=checks,
        scope='Independent reconstruction from the bound initial/native arrays and saved primitive/field state; recompute mass, lapse, frame maps, fields, fluxes and all native finite residuals. Exact native query keys are required. No new EOS evaluations, time-convergence or continuum proof.')
    e.write(balanced.OUT/'replay.json',report)
    print(json.dumps(report),flush=True)


if __name__=='__main__':verify()
