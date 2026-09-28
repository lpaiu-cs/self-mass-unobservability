"""Counterexample candidate: isolate outer exchange in the actual GR evolution.

Mirror radial vector fields at the outer face: the shared baryon, energy and
heat fluxes vanish there, while the even radial momentum flux is retained.
This is a reflecting computational wall, not a stellar atmosphere or vacuum.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import subprocess

import numpy as np

import gr_implicit_second_heat_time as prior

old, e = prior.old, prior.e
OUT = e.g.OUT/'gr-implicit-reflecting-boundary'


class ReflectingStar(old.prior.BaryonStar):
    def faces(self, value, odd=False):
        result = super().faces(value, odd=odd)
        if odd:
            result[-1] = 0
        return result


class ReflectingTangent(old.LocalTangent, ReflectingStar):
    pass


def initialize(pool):
    star = old.initialize(pool)
    star.__class__ = ReflectingStar
    return star


def tangent(star, delta, z):
    model = old.tangent(star, delta, z)
    model.__class__ = ReflectingTangent
    return model


# Reuse the bound native stage, time loop and budget, with the SAME boundary
# in the actual RHS and its approximate Jacobian. No edit of a prior source.
engine = SimpleNamespace(**dict(vars(old), initialize=initialize, tangent=tangent))
solver = SimpleNamespace(**vars(prior.solver))
solver.stage = FunctionType(prior.solver.stage.__code__, dict(vars(prior.solver), old=engine))
namespace = dict(vars(prior), OUT=OUT, old=engine, solver=solver,
                 second_tau=lambda: prior.FIRST_TAU)
run = FunctionType(prior.run.__code__, namespace)


def check():
    e.TAU = prior.FIRST_TAU
    software = FunctionType(old.check.__code__, dict(vars(old), initialize=initialize, tangent=tangent))()
    star = initialize(None)
    delta = np.zeros_like(star.base)
    _, state = star.rhs(delta)
    assert np.array_equal(state['boundary_rates'], np.zeros(2))
    for value in [state['N']*state['D']*state['v'], state['N']*state['S']/state['a'],
                  state['N']*state['Q']/state['a']]:
        assert star.faces(value, odd=True)[-1] == 0
    assert prior.analysis.cones(state)['sampled_cone_inside_light_cone']
    modes = {}
    # One bounded regression: the initial outer 32-cell frozen-tangent heat
    # block. This is not the spectrum or stability proof of the full solution.
    selected = np.ravel(np.column_stack([np.arange(1, 128, 4), np.arange(3, 128, 4)]))
    for name, initial, linearize in [('old_exchange', old.initialize, old.tangent),
                                     ('reflecting_wall', initialize, tangent)]:
        model = initial(None)
        _, z = model.rhs(delta)
        matrix = old.jacobian(linearize(model, delta, z), delta)[-128:, -128:].toarray()
        roots = np.linalg.eigvals(matrix[np.ix_(selected, selected)])
        modes[name] = float(roots.real.max())
    assert modes['old_exchange'] > 100 and modes['reflecting_wall'] < 0, modes
    result = dict(classification='Counterexample candidate', implementation=software,
                  outer_frozen_heat_maximum_real_eigenvalues=modes,
                  full_nonlinear_stability_proven=False, physical_exterior_match=False)
    print('PASS shared reflecting flux, actual tangent and bounded heat-mode regression', result, flush=True)
    return result


def prepare():
    assert not OUT.exists(), 'Keep every earlier boundary and heat-time candidate.'
    checked = check()
    plan = json.loads((prior.OUT/'plan.json').read_text())
    for rel, digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    plan.update(checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        parent_plan_sha256=e.digest(prior.OUT/'plan.json'),
        candidate='Original constant proper heat time and initial star, with a reflecting computational outer wall.',
        proper_tau_seconds=str(prior.FIRST_TAU),
        boundary='Outer face: mirror radial vector fields so shared baryon, coordinate-energy and heat fluxes vanish. Retain even radial momentum flux. The wall can exert stress; this is not an atmosphere, vacuum match or luminosity prediction.',
        implementation_check=checked)
    for path in [Path(__file__), Path(prior.__file__), prior.OUT/'plan.json']:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    print('PREPARED REFLECTING BOUNDARY NATIVE EVOLUTION', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['check', 'prepare', 'run'])
    parser.add_argument('--refinement', type=int, choices=[1, 2, 4], default=1)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.command == 'run':
        run(args.refinement, args.workers)
    else:
        globals()[args.command]()
