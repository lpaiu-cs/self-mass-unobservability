"""Counterexample candidate: two-carrier evolution with covariant momentum.

The immutable parent omitted one a_t*S term from the primitive momentum row.
Keep the parent trajectories as numerical controls; restart this correction
from the same initial state, EOS, transport laws, time mesh and tolerances.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import subprocess
import sys

import numpy as np
import sympy as sp

import gr_two_carrier_evolution as parent
import gr_two_carrier_analysis as audit

e, ld, wall = parent.e, parent.ld, parent.wall
OUT = e.g.OUT/'gr-covariant-two-carrier'
cones, energy_budget = parent.cones, parent.energy_budget


class CovariantStar(parent.TwoCarrierStar):
    def rhs(self, delta):
        y = self.base+delta
        z = self.state(y)
        rho,T,P,w,v,Q,W,D,E,S,R,a,N,nur,at = [z[k] for k in
            ['rho','T','P','w','v','Q','W','D','E','S','R','a','N','nur','at']]
        aux = z['aux']
        speed = e.C*N*v/a
        fB, fE, fS = self.faces(N*D*v, odd=True), self.faces(N*S/a, odd=True), self.faces(N*R)
        Bdot = -e.C*self.divergence(fB)-at*D
        Sdot = -e.C*self.divergence(fS)-2*at*S-N*nur*e.C*E+e.C*N*P*self.area_difference_over_volume
        matrix = np.zeros((self.n, 5, 5), dtype=ld)
        right = np.zeros((self.n, 5), dtype=ld)
        matrix[:, 0, 0], matrix[:, 0, 2] = a*D, a*rho*W**3*v
        right[:, 0] = Bdot+a*speed*self.gradient(D)
        matrix[:, 1, :4] = np.column_stack([rho*aux[:, 9]-P, rho*aux[:, 10], 2*Q*W*W, v*self.qscale])
        right[:, 1] = -e.C*N/(a*W*W)*self.gradient(Q, odd=True)-2*Q*(e.C*N*nur/a+v*at/a+e.C*N/(a*self.r))
        er, et = rho*(z['rest']+z['u']+aux[:, 9]), rho*aux[:, 10]
        pr, pt = P*aux[:, 5], P*aux[:, 6]
        matrix[:, 2, :4] = a[:, None]*np.column_stack([W*W*v*(er+pr), W*W*v*(et+pt),
                               W**4*((1+v*v)*w+4*Q*v), (1+v*v)*W*W*self.qscale])
        right[:, 2] = Sdot+a*speed*self.gradient(S, odd=True)
        force = e.C*self.gradient(y[:, 1])/(a*W)+e.C*W*nur/a+W*v*at/(N*a)
        for row, q, K, tau, kr, kt in [
            (3, z['Qcond'], z['Kcond'], e.TAU, -aux[:, 28], 5-aux[:, 29]),
            (4, z['Qrad'], z['Krad'], z['tau_rad'], 1, 5)]:
            h = K*T/e.C**2
            matrix[:, row, :3] = (W/N)[:, None]*np.column_stack(
                [-tau*q*kr/2, h*v-tau*q*kt/2, h*W*W])
            matrix[:, row, row] = W/N*tau*self.qscale
            if row == 3:
                matrix[:, row, 4] = -W/N*tau*self.qscale
            right[:, row] = -q-h*force
        scales = np.max(abs(matrix), axis=2)
        A, b = np.asarray(matrix/scales[:, :, None], float), np.asarray(right/scales, float)
        rates = np.linalg.solve(A, b[:, :, None])[:, :, 0].astype(ld)
        residual = float(np.max(abs(np.einsum('nij,nj->ni', A, rates)-b)/(1+abs(b))))
        assert residual < 1e-10, residual
        out = np.zeros_like(y)
        for j in range(5):
            spatial = (self.gradient(z['Q'] if j == 3 else z['Qrad'], odd=True)/self.qscale
                       if j >= 3 else self.gradient(y[:, j], odd=(j == 2)))
            out[:, j] = rates[:, j]-speed*spatial
        left = np.vstack([y[0, 5:]*0, np.diff(y[:, 5:], axis=0)/np.diff(self.r)[:, None]])
        right_x = np.vstack([np.diff(y[:, 5:], axis=0)/np.diff(self.r)[:, None], y[-1, 5:]*0])
        out[:, 5:] = -speed[:, None]*np.where(speed[:, None] >= 0, left, right_x)
        z['time_matrix_residual'] = residual
        z['boundary_rates'] = np.array([-e.C*4*np.pi*self.rf[-1]**2*fB[-1], -e.C*4*np.pi*self.rf[-1]**2*fE[-1]])
        return out, z


class CovariantTangent(parent.LocalTangent, CovariantStar):
    pass


def initialize(pool):
    star = parent.initialize(pool)
    star.__class__ = CovariantStar
    return star


def tangent(star, delta, z):
    model = parent.tangent(star, delta, z)
    model.__class__ = CovariantTangent
    return model


engine = SimpleNamespace(**dict(vars(parent.engine), initialize=initialize, tangent=tangent))
solver = SimpleNamespace(**dict(vars(parent.solver),
    stage=FunctionType(parent.stage.__code__, dict(vars(parent), engine=engine))))
run = FunctionType(parent.run.__code__, dict(parent.namespace, OUT=OUT, old=engine, solver=solver))


def check():
    # Tensor balance, independently reduced to the primitive momentum row.
    t, r = sp.symbols('t r', positive=True)
    a, N, J, R, P, E = [sp.Function(k)(t, r) for k in ['a','N','J','R','P','E']]
    covariant = (sp.diff(a*a*r*r*J, t)+sp.diff(N*a*r*r*R, r)
        -N*a*r*r*(R*sp.diff(a, r)/a+2*P/r-E*sp.diff(N, r)/N))
    primitive = (a*sp.diff(J, t)+2*sp.diff(a, t)*J
        +sp.diff(r*r*N*R, r)/r**2+E*sp.diff(N, r)-2*N*P/r)
    assert sp.simplify(covariant/(a*r*r)-primitive) == 0
    assert sp.simplify(covariant/(a*r*r)-(primitive-sp.diff(a, t)*J)) == J*sp.diff(a, t)
    # Exercise BOTH production RHS methods on identical native EOS states.
    # The velocity is a manufactured control, not a new physical initial star.
    star = initialize(None)
    star.base[:, 2] = ld('.001')
    star.base[:, 3:5] = 0
    residuals = []
    for method in [parent.TwoCarrierStar.rhs, CovariantStar.rhs]:
        rate, z = method(star, np.zeros_like(star.base))
        a,N,W,v,Q,P,S,at = [z[k] for k in ['a','N','W','v','Q','P','S','at']]
        er = z['rho']*(z['rest']+z['u']+z['aux'][:, 9])
        et = z['rho']*z['aux'][:, 10]
        pr,pt = P*z['aux'][:, 5], P*z['aux'][:, 6]
        time_row = a[:, None]*np.column_stack([W*W*v*(er+pr), W*W*v*(et+pt),
            W**4*((1+v*v)*z['w']+4*Q*v), (1+v*v)*W*W*star.qscale])
        speed = e.C*N*v/a
        material_rate = rate[:, :4]+speed[:, None]*np.column_stack([
            star.gradient(star.base[:, 0]), star.gradient(star.base[:, 1]),
            star.gradient(v, odd=True), star.gradient(Q, odd=True)/star.qscale])
        expected = (-e.C*star.divergence(star.faces(N*z['R']))-2*at*S
            -N*z['nur']*e.C*z['E']+e.C*N*P*star.area_difference_over_volume
            +a*speed*star.gradient(S, odd=True))
        residuals.append(np.sum(time_row*material_rate, axis=1)-expected)
    scale = np.max(abs(at*S))
    corrected = float(np.max(abs(residuals[1]))/scale)
    omitted = float(np.max(abs(residuals[0]-at*S))/scale)
    assert scale > 0 and corrected < 1e-6 and omitted < 1e-6, (corrected, omitted)
    assert np.max(abs(residuals[0]))/scale > .99
    native = FunctionType(parent.check.__code__,
        dict(vars(parent), initialize=initialize, tangent=tangent))()
    result = dict(classification='Counterexample candidate', passed=True,
        tensor_identity=dict(classification='Proven',
            primitive_momentum='a*d_t S = -c*div(N*R)-2*a_t*S-c*N*nu_r*E+2*c*N*P/r',
            scope='Smooth spherical covariant stress balance, independent of an EOS closure.'),
        corrected_production_residual_over_missing_term=corrected,
        legacy_missing_term_reproduction_error=omitted,
        parent_implementation_check=native,
        full_nonlinear_trajectory=False, physical_EOS_certified=False)
    print('PASS covariant momentum production RHS and native tangent', json.dumps(result), flush=True)
    return result


def prepare():
    assert not OUT.exists(), 'Preserve earlier candidates and their frozen comparisons.'
    checked = check()
    plan = json.loads((parent.OUT/'plan.json').read_text())
    for rel, digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    plan.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
        candidate='Two heat carriers with the corrected covariant primitive momentum row.',
        correction='Replace -a_t*S by -2*a_t*S in the primitive momentum RHS only.',
        parent_plan_sha256=e.digest(parent.OUT/'plan.json'),
        implementation_check=checked,
        inherited_numerical_limits='Same primitive spatial approximation and measured, nonzero energy/species budgets. This correction is not an exact finite-volume conservation repair.',
        previous_paths='The immutable parent has a covariant momentum source omission. Preserve its trajectories as controls; do not call them completed physical GR evolutions.')
    for path in [Path(__file__), Path(audit.__file__), parent.OUT/'plan.json']:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    comparison['prepare']()
    print('PREPARED COVARIANT TWO-CARRIER EVOLUTION AND FULL COMPARISON', flush=True)


# Reuse the already frozen five-field full-path verifier and convergence gate.
comparison = dict(vars(audit), evolution=sys.modules[__name__], OUT=OUT, __file__=__file__)
for name in ['bindings', 'completed', 'compare', 'prepare']:
    comparison[name] = FunctionType(getattr(audit, name).__code__, comparison)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['check','prepare','run','completed','compare'])
    parser.add_argument('--refinement', type=int, choices=[1,2,4], default=1)
    parser.add_argument('--workers', type=int, default=15)
    args = parser.parse_args()
    if args.command == 'run':
        run(args.refinement, args.workers)
    elif args.command == 'completed':
        print(json.dumps(comparison['completed'](args.refinement)[1], indent=2))
    elif args.command == 'compare':
        comparison['compare']()
    else:
        globals()[args.command]()
