"""Use the shared face heat flux in the actual material energy derivative.

Counterexample candidate. Reuses the unchanged EOS, time matrix and SSPRK2
driver. This corrects the heat-divergence mismatch, not every spatial defect.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import subprocess

import numpy as np
import sympy as sp
import gr_conservative_evolution as prior
import gr_coupled_evolution as original

OUT = original.g.OUT/'gr-flux-coupled-evolution'
write, digest = original.write, original.digest


class FluxStar(original.Star):
    def state(self, y):
        z = super().state(y)
        self.current = z
        return z

    def gradient(self, value, odd=False):
        # Only the physical Q object returned by the current state takes this
        # path. Other gradients and both uses of Q_r retain their declared role.
        z = getattr(self, 'current', None)
        if z is not None and value is z['Q']:
            a, N, Q = z['a'], z['N'], z['Q']
            loga_r = a*a*(4*np.pi*original.GRAV*self.r*z['E']-z['m']/self.r**2)
            return a/N*self.divergence(self.faces(N*Q/a, odd=True))-Q*(z['nur']-loga_r+2/self.r)
        return super().gradient(value, odd=odd)


def symbolic():
    r = sp.symbols('r', positive=True)
    a, N, Q = [sp.Function(k)(r) for k in ['a', 'N', 'Q']]
    divergence = sp.diff(r*r*N*Q/a, r)/(r*r)
    reconstructed = a/N*divergence-Q*(sp.diff(N, r)/N-sp.diff(a, r)/a+2/r)
    assert sp.simplify(reconstructed-sp.diff(Q, r)) == 0
    # At v=0, the mass/lapse constraints turn the geometrical heat terms
    # into w*a_t/a. Baryon compression then gives E_t=-c*div(NQ/a).
    nu, ar, w, k, c, q, av, nv = sp.symbols('nu ar w k c q av nv')
    at = -4*sp.pi*k*c*r*nv*av**2*q
    remainder = c*nv*q/av*(nu-ar+2/r)-2*c*nv*q/av*(nu+1/r)-w*at/av
    assert sp.expand(remainder.subs(ar, 4*sp.pi*k*r*av**2*w-nu)) == 0
    assert original.symbolic()['passed']
    return dict(classification='Proven', passed=True,
        scope='Continuous product identity and conditional v=0 cancellation under radial Einstein constraints. No claim of exact finite-velocity discrete conservation.')


def prepare():
    previous = prior.verify_initial()
    assert not OUT.exists()
    OUT.mkdir()
    inputs = [Path(__file__), Path(prior.__file__), Path(original.__file__), prior.OUT/'initial-manifest.json']
    plan = dict(previous, checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        bindings=dict(previous['bindings'], **{p.relative_to(original.ROOT).as_posix(): digest(p) for p in inputs}),
        evolution='Same conservative initial restriction, native EOS, material time matrix, composition advection and SSPRK2. For physical Q_r only, use a/N*FVdiv(NQ/a)-Q*(nu_r-a_r/a+2/r); a_r/a comes from the current radial energy constraint. At v=0 this aligns the heat energy equation with the shared coordinate-energy face flux. Other finite-velocity/composition/spatial chain-rule defects remain measurable.',
        symbolic=symbolic())
    for name in ['initial.npz', 'restriction.json']:
        (OUT/name).write_bytes((prior.OUT/name).read_bytes())
    write(OUT/'plan.json', plan)
    write(OUT/'initial-manifest.json', dict(sha256={p.name: digest(p) for p in
          [OUT/'initial.npz', OUT/'restriction.json', OUT/'plan.json']}))
    verify_initial()


# Reuse the actual driver; only its explicitly bound Star implementation and
# output directory differ. No copied time loop or runtime edit of prior modules.
engine = SimpleNamespace(**vars(original))
engine.Star = FluxStar
namespace = dict(vars(prior), OUT=OUT, evolution=engine)
verify_initial = FunctionType(prior.verify_initial.__code__, namespace)
namespace['verify_initial'] = verify_initial
run = FunctionType(prior.run.__code__, namespace)
compare = FunctionType(prior.compare.__code__, namespace)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'run', 'compare', 'verify_initial'])
    parser.add_argument('--steps', type=int, default=16)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    if args.command == 'run': run(args.steps, args.workers)
    else: globals()[args.command]()
