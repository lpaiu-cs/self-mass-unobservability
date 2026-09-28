"""Counterexample candidate: native-residual SDIRK2 for the existing coupled MOL.

The sparse chord matrix freezes geometry and EOS derivatives only during the
linear solve. Every accepted nonlinear stage uses the unmodified native RHS.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import subprocess
import time

import numpy as np
import sympy as sp
from scipy.linalg import solve_banded
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import splu

import gr_baryon_flux_evolution as prior

e = prior.e
ld = e.ld
OUT = e.g.OUT/'gr-implicit-coupled-evolution'
GAMMA = 1-1/np.sqrt(ld(2))
SCALE = np.array([1e-6, .1, 1e-8, 1e-10], dtype=ld)
ATOL = np.array([1e-13, 1e-9, 1e-14, 1e-17]+[1e-13]*26, dtype=ld)


class CachedOnly:
    def map(self, function, rows, **kwargs):
        assert not list(rows), 'A tangent/software check must never call native EOS.'
        return iter(())


def symbolic():
    q = 1-1/sp.sqrt(2)
    z = sp.symbols('z')
    A = sp.Matrix([[q, 0], [1-q, q]])
    b = sp.Matrix([[1-q, q]])
    assert sp.simplify((b*A*sp.ones(2, 1))[0]-sp.Rational(1, 2)) == 0
    R = sp.cancel(1+z*(b*(sp.eye(2)-z*A).inv()*sp.ones(2, 1))[0])
    assert sp.limit(R, z, -sp.oo) == 0
    # On the imaginary axis, |denominator|^2-|numerator|^2 = q^4*y^4.
    y = sp.symbols('y', real=True)
    assert sp.simplify(abs((1-q*sp.I*y)**2)**2
                       -abs(1+(1-2*q)*sp.I*y)**2-q**4*y**4) == 0
    assert prior.symbolic()['passed']
    return dict(classification='Proven', passed=True,
        scope='SDIRK tableau order two and rational L stability; not nonlinear, spatial, native EOS or continuum error certification.')


def initialize(pool):
    data = np.load(prior.OUT/'initial.npz')
    star = prior.BaryonStar(len(data['base']), CachedOnly() if pool is None else pool)
    assert np.array_equal(star.r, data['radius'])
    assert np.array_equal(star.volume, data['volume'])
    star.base, star.qscale = data['base'].copy(), data['qscale'].copy()
    star.initial_aux = data['aux'].copy()
    star.material_cache = {e.material_key(row): value for row, value in zip(
        zip(star.base[:, 0], star.base[:, 1], star.base[:, 4:]), data['aux'])}
    return star


class LocalTangent(prior.BaryonStar):
    """Approximate Jacobian only; never used to accept a time stage."""
    def state(self, y):
        d = y-self.anchor
        aux = self.reference['aux'].copy()
        aux[:, 1] *= np.exp(aux[:, 5]*d[:, 0]+aux[:, 6]*d[:, 1])
        aux[:, 2] += aux[:, 9]*d[:, 0]+aux[:, 10]*d[:, 1]
        aux[:, 21] *= np.exp(aux[:, 22]*d[:, 0]+aux[:, 23]*d[:, 1])
        self.material_cache = {e.material_key(row): value for row, value in zip(
            zip(y[:, 0], y[:, 1], y[:, 4:]), aux)}
        z = super().state(y)
        for key in ['mf', 'm', 'a', 'N', 'nur']:
            z[key] = self.reference[key]
        z['at'] = -4*np.pi*e.GRAV*e.C*self.r*z['N']*z['a']**2*z['S']
        return z


def tangent(star, delta, z):
    model = LocalTangent.__new__(LocalTangent)
    model.__dict__ = dict(star.__dict__)
    model.pool = CachedOnly()
    model.anchor, model.reference = star.base+delta, z
    return model


def jacobian(model, delta):
    # Five colors include the three-node one-sided outer gradient. The frozen
    # geometry removes global constraint coupling from this preconditioner.
    n = model.n
    rows, columns, values = [], [], []
    step = ld('0.001')
    for field in range(4):
        for color in range(5):
            selected = np.arange(color, n, 5)
            change = np.zeros_like(delta)
            change[selected, field] = step*SCALE[field]
            response = (model.rhs(delta+change)[0][:, :4]
                        -model.rhs(delta-change)[0][:, :4])/(2*step*SCALE)
            for offset in range(-2, 3):
                target = selected+offset
                valid = (target >= 0) & (target < n)
                target, source = target[valid], selected[valid]
                rows.extend((4*target[:, None]+np.arange(4)).ravel())
                columns.extend(np.repeat(4*source+field, 4))
                values.extend(np.asarray(response[target], float).ravel())
    J = coo_matrix((values, (rows, columns)), shape=(4*n, 4*n)).tocsc()
    J.eliminate_zeros()
    return J


def species(star, stage_base, h, y):
    z = star.current
    speed = e.C*z['N']*y[:, 2]/z['a']
    dr = np.diff(star.r)
    left = np.r_[ld(0), np.maximum(speed[1:], 0)/dr]
    right = np.r_[np.maximum(-speed[:-1], 0)/dr, ld(0)]
    band = np.zeros((3, star.n))
    band[1] = 1+h*(left+right)
    band[0, 1:], band[2, :-1] = -h*right[:-1], -h*left[1:]
    differences = np.diff(star.base[:, 4:], axis=0)
    drive = left[:, None]*np.vstack([np.zeros(26), -differences])
    drive += right[:, None]*np.vstack([differences, np.zeros(26)])
    return solve_banded((1, 1), band, np.asarray(stage_base[:, 4:]+h*drive, float)).astype(ld)


def stage(star, stage_base, seed, h, log):
    delta = seed.copy()
    factor = None
    for iteration in range(16):
        rates, z = star.rhs(delta)
        residual = delta-stage_base-h*rates
        tolerance = ATOL+ld('1e-7')*np.maximum(abs(delta), abs(stage_base))
        norm = float(np.max(abs(residual)/tolerance))
        record = dict(iteration=iteration, residual_norm=norm,
            maximum_absolute_residual=np.max(abs(residual), axis=0).astype(float).tolist())
        log(record)
        if norm <= 1:
            return delta, rates, z
        if factor is None or iteration % 4 == 0:
            J = jacobian(tangent(star, delta, z), delta)
            factor = splu(eye(4*star.n, format='csc')-float(h)*J)
        correction = factor.solve(np.asarray(-residual[:, :4]/SCALE, float).ravel())
        correction = correction.reshape(star.n, 4).astype(ld)*SCALE
        # Limit the nonlinear iterate, never clip an accepted physical solution.
        limit = min(1., .1/max(float(np.max(abs(correction[:, :2]))), 1e-300),
                    .01/max(float(np.max(abs(correction[:, 2]))), 1e-300))
        delta[:, :4] += limit*correction
        delta[:, 4:] = species(star, stage_base, h, star.base+delta)
    raise RuntimeError(('Native implicit stage did not converge', norm))


def energy_budget(star, delta, z, mass_flux):
    # Same stable native endpoint increment as the completed explicit paths.
    # Never obtain local thermal errors by subtracting large mass prefixes.
    rho0 = np.exp(star.base[:, 0])
    rest0 = (star.base[:, 4:]/e.g.c.A)@e.g.c.W*e.C**2
    drest = (delta[:, 4:]/e.g.c.A)@e.g.c.W*e.C**2
    aux0 = star.initial_aux
    eps0 = rho0*(rest0+aux0[:, 2])
    deps = rho0*((rest0+aux0[:, 2])*np.expm1(delta[:, 0])
                +np.exp(delta[:, 0])*(drest+z['u']-aux0[:, 2]))
    v, Q = z['v'], z['Q']
    dE = (deps+(z['P']+eps0)*v*v+2*Q*v)/(1-v*v)
    heat = rho0*aux0[:, 10]*star.volume
    defect = dE*star.volume+np.diff(mass_flux)/e.GRAV
    return defect/heat, dict(
        maximum_cell_energy_defect_over_initial_heat_capacity=float(np.max(abs(defect/heat))),
        total_energy_defect_over_initial_heat_capacity=float(defect.sum()/heat.sum()))


def check():
    assert symbolic()['passed']
    star = initialize(None)
    delta = np.zeros_like(star.base)
    _, z = star.rhs(delta)  # Initial exact native outputs are already frozen.
    assert np.all(energy_budget(star, delta, z, np.zeros(star.n+1, dtype=ld))[0] == 0)
    model = tangent(star, delta, z)
    J = jacobian(model, delta)
    rng = np.random.default_rng(33)
    direction = rng.uniform(-1, 1, size=(star.n, 4)).astype(ld)
    shift = np.zeros_like(delta)
    step = ld('0.0001')
    shift[:, :4] = step*direction*SCALE
    finite = (model.rhs(delta+shift)[0][:, :4]-model.rhs(delta-shift)[0][:, :4])/(2*step*SCALE)
    predicted = J.dot(np.asarray(direction, float).ravel()).reshape(star.n, 4)
    error = float(np.max(abs(finite-predicted))/np.max(abs(finite)))
    assert error < 1e-4, error
    # Exercise the tableau on a stiff linear equation; this software control
    # does not test convergence of the actual nonlinear stage iteration.
    lam = ld(-10000)
    errors = []
    for steps in [8, 16, 32]:
        h = ld('.0001')/steps
        value = ld(1)
        for _ in range(steps):
            first = value/(1-GAMMA*h*lam)
            value = (value+(1-GAMMA)*h*lam*first)/(1-GAMMA*h*lam)
        errors.append(float(abs(value-np.exp(lam*ld('.0001')))))
    orders = np.log2(np.array(errors[:-1])/errors[1:])
    assert np.all(orders > 1.9)
    result = dict(classification='Counterexample candidate', passed=True,
        frozen_local_tangent_direction_error=error, linear_test_orders=orders.tolist(),
        actual_nonlinear_trajectory=False, physical_EOS_certified=False)
    print('IMPLICIT IMPLEMENTATION CHECK', json.dumps(result), flush=True)
    return result


def prepare():
    previous = prior.verify_initial()
    assert not OUT.exists(), 'Preserve every old plan and failed attempt.'
    checks = check()
    OUT.mkdir()
    # Fixed mesh in coordinate time. Refinement bisects every interval, keeping
    # the same end time and the entire initial transient in each path.
    edges = [ld(0)]
    while edges[-1] < 1024:
        edges.append(min(ld(1024), edges[-1]+min(ld(32), ld('.25')*2**(len(edges)-1))))
    e.write(OUT/'plan.json', dict(classification='Counterexample candidate',
        checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        bindings=dict(previous['bindings'], **{p.relative_to(e.ROOT).as_posix():e.digest(p)
            for p in [Path(__file__), Path(prior.__file__), prior.OUT/'initial-manifest.json']}),
        runtime_sha256=previous['runtime_sha256'], cells=5735,
        duration_tau=1024, duration_seconds=float(1024*e.TAU),
        time_edges_tau=[float(v) for v in edges], refinements=[1, 2, 4],
        nonlinear_absolute_tolerances=ATOL.astype(float).tolist(), nonlinear_relative_tolerance=1e-7,
        maximum_stage_iterations=16, symbolic=symbolic(), implementation_check=checks,
        method='L-stable order-two SDIRK; native EOS, opacity, all primitive equations, passive composition and reconstructed mass/lapse in each accepted stage residual. Sparse local frozen-geometry/EOS-derivative chord matrix is only a preconditioner. Passive species use the current velocity in an implicit upwind solve and are checked in the full residual.',
        refusal='An unconverged stage stops this immutable path; no residual relaxation or unrecorded dt reduction.',
        time_refinement_gate='All four endpoint maximum differences decrease; lnT and Q observed orders at least 1.5.',
        boundary='Inherited finite-pressure computational outer face; no physical exterior match.',
        limitations='Native nonlinear residual tolerance is not a rigorous error enclosure. Existing spatial operator, heat parameter and advected nonreacting composition are unchanged; no full physical EOS, exterior, stellar lifetime or observational closure.'))
    print('IMPLICIT COUPLED PLAN', len(edges)-1, 'base steps to', float(1024*e.TAU), 's', flush=True)


def run(refinement, workers):
    plan = json.loads((OUT/'plan.json').read_text())
    assert refinement in plan['refinements']
    for rel, digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    for rel, digest in plan['runtime_sha256'].items():
        assert e.digest(Path(rel)) == digest, rel
    folder = OUT/f'path-{refinement}'
    assert not folder.exists()
    folder.mkdir()
    edges = np.array(plan['time_edges_tau'], dtype=ld)*e.TAU
    times = np.r_[ld(0), np.concatenate([np.linspace(a, b, refinement+1, dtype=ld)[1:]
                                        for a, b in zip(edges[:-1], edges[1:])])]
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers, initializer=e.worker_init) as pool:
        star = initialize(pool)
        delta = np.zeros_like(star.base)
        exchange = np.zeros(2, dtype=ld)
        mass_flux = np.zeros(star.n+1, dtype=ld)
        for step, t in enumerate(times):
            rates, z = star.rhs(delta)
            cell_defect, budget = energy_budget(star, delta, z, mass_flux)
            np.savez_compressed(folder/f'step-{step:04d}.npz', delta=delta,
                **{k:z[k] for k in ['m', 'mf', 'a', 'N', 'Q', 'aux']},
                time_seconds=t, boundary_exchange=exchange, integrated_face_mass_flux=mass_flux,
                normalized_cell_energy_defect=cell_defect)
            record = dict(classification='Counterexample candidate', step=step,
                time_seconds=float(t), maximum_changes=np.max(abs(delta[:, :4]), axis=0).astype(float).tolist(),
                elapsed_seconds=time.monotonic()-started, **budget)
            e.write(folder/'progress.json', record)
            print('IMPLICIT COUPLED STEP', refinement, json.dumps(record), flush=True)
            if step == len(times)-1:
                break
            dt = times[step+1]-t
            def log(row):
                with (folder/'iterations.jsonl').open('a') as stream:
                    stream.write(json.dumps(dict(step=step+1, stage=stage_number, **row))+'\n')
                print('NATIVE IMPLICIT STAGE', refinement, step+1, stage_number,
                      row['iteration'], row['residual_norm'], flush=True)
            try:
                stage_number = 1
                first, f1, z1 = stage(star, delta, delta, GAMMA*dt, log)
                stage_number = 2
                delta, _, z2 = stage(star, delta+(1-GAMMA)*dt*f1, first, GAMMA*dt, log)
            except Exception as error:
                e.write(folder/'failure.json', dict(classification='Counterexample candidate',
                    completed=False, step=step+1, stage=stage_number, reason=repr(error)))
                raise
            exchange += dt*((1-GAMMA)*z1['boundary_rates']+GAMMA*z2['boundary_rates'])
            flux = ((1-GAMMA)*star.faces(z1['N']*z1['S']/z1['a'], odd=True)
                    +GAMMA*star.faces(z2['N']*z2['S']/z2['a'], odd=True))
            mass_flux += dt*e.GRAV*e.C*4*np.pi*star.rf**2*flux
        e.write(folder/'result.json', dict(classification='Counterexample candidate', completed=True,
            steps=len(times)-1, duration_seconds=float(times[-1]), actual_native_nonlinear_path=True,
            physical_EOS_certified=False, physical_exterior_match=False, observational_closure=False))
    e.write(folder/'manifest.json', dict(plan_sha256=e.digest(OUT/'plan.json'),
        sha256={p.name:e.digest(p) for p in folder.iterdir() if p.is_file()}))


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
