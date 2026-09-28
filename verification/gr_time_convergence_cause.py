"""Counterexample candidate: diagnose the frozen full-time convergence failure.

This is an affine fluid-sector control, not a replacement GR evolution.  The
actual finite-volume fluxes, gravity reconstruction and native rho/T EOS
derivatives are retained. Heat flux and composition perturbations are frozen.
The control must be compared against the saved full-GR central response.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import os
import subprocess

import numpy as np
from scipy.linalg import eig, lu_factor, lu_solve

import gr_conservative_composition_tangent as gr

ROOT, ld = gr.e.ROOT, gr.ld
SOURCE = ROOT/'outputs/direct-eos-gr33/gr-conservative-composition-recovery-20260914'
FIELDS = ['lnrho_B', 'lnT', 'v/c']
SCALE = gr.SCALE[:3]
EPS = ld('.001')
STAR, ZERO, Z0, NC = None, None, None, None


class FrozenThermodynamics(gr.prior.ConservativeStar):
    def evaluate(self, delta):
        aux = self.initial_aux.copy()
        aux[:, 1] *= np.exp(aux[:, 5]*delta[:, 0]+aux[:, 6]*delta[:, 1])
        aux[:, 2] += aux[:, 9]*delta[:, 0]+aux[:, 10]*delta[:, 1]
        for j in [24, 27]:
            aux[:, j] *= np.exp(aux[:, j+1]*delta[:, 0]+aux[:, j+2]*delta[:, 1])
        aux[:, 21] = aux[:, 24]*aux[:, 27]/(aux[:, 24]+aux[:, 27])
        y = self.base+delta
        self.material_cache = {gr.e.material_key(row): value for row, value in
                               zip(zip(y[:, 0], y[:, 1], y[:, 5:]), aux)}
        return super().evaluate(delta)


def initialize(cells):
    global STAR, ZERO, Z0, NC
    gr.OUT = SOURCE
    gr.bindings()
    STAR = gr.initialize(None)
    STAR.__class__ = FrozenThermodynamics
    ZERO = np.zeros_like(STAR.base)
    Z0 = STAR.evaluate(ZERO)
    NC = cells


def rows(value):
    answer = value[:NC, :3].copy()
    # Exact row operation removes the large common rest-energy contribution
    # before conversion to binary64; it does not change the linear system.
    rest_ratio = Z0['rho']*(Z0['rest']+Z0['u'])/STAR.heat0
    answer[:, 1] -= rest_ratio[:NC]*answer[:, 0]
    return (answer/SCALE).ravel()


def evaluate(x):
    delta = ZERO.copy()
    delta[:NC, :3] = np.asarray(x, dtype=ld).reshape(NC, 3)*SCALE
    z = STAR.evaluate(delta)
    spatial, _ = gr.residual(STAR, delta, (delta, z), (delta, z), ld(1), gr.weights(ld(1), None))
    mass, _ = gr.residual(STAR, delta, (ZERO, Z0), (ZERO, Z0), ld('1e-40'), gr.weights(ld(1), None))
    return rows(spatial), rows(mass)


def column(j):
    x = np.zeros(3*NC, dtype=ld)
    x[j] = EPS
    plus, mp = evaluate(x)
    minus, mm = evaluate(-x)
    return j, np.asarray((plus-minus)/(2*EPS), float), np.asarray((mp-mm)/(2*EPS), float)


def native_probe():
    gr.e.worker_init()
    samples = []
    selected = [0, 1, 15, 3119, 3142, 3143, 5704, 5734]
    for refinement, step in [(4, 0), (1, 39), (2, 78), (4, 156)]:
        with np.load(SOURCE/f'path-{refinement}'/f'step-{step:04d}.npz') as cp:
            y = STAR.base+cp['delta']
            for i in selected:
                actual = gr.parent.parent.material((y[i, 0], y[i, 1], y[i, 5:]))
                saved = cp['aux'][i]
                samples.append(dict(refinement=refinement, step=step, cell=i,
                    exact_equal=bool(np.array_equal(actual, saved)),
                    P_relative=float(abs(actual[1]-saved[1])/abs(saved[1])),
                    u_over_initial_cvT=float(abs(actual[2]-saved[2])/STAR.initial_aux[i, 10]),
                    maximum_relative=float(np.max(abs(actual-saved)/np.maximum(abs(saved), ld('1e-300'))))))
    return samples


def integrate(M, K, b, plan, refinement):
    times = gr.parent.wall.prior.time_nodes(plan, refinement)
    previous = older = np.zeros(len(b))
    cache, history, worst = {}, [], 0.
    for step in range(1, len(times)):
        h = times[step]-times[step-1]
        c0, c1, c2 = map(float, gr.weights(h, None if step == 1 else times[step-1]-times[step-2]))
        hf = float(h)
        key = c0, hf
        if key not in cache:
            matrix = c0*M+hf*K
            scale = np.max(abs(matrix), axis=1)
            assert np.all(scale > 0)
            equilibrated = matrix/scale[:, None]
            cache[key] = lu_factor(equilibrated), scale, equilibrated
        factor, scale, matrix = cache[key]
        rhs = (-c1*(M@previous)-c2*(M@older)-hf*b)/scale
        current = lu_solve(factor, rhs)
        relative = float(np.max(abs(matrix@current-rhs))/(1+np.max(abs(rhs))))
        worst = max(worst, relative)
        assert np.isfinite(current).all() and relative < 1e-10
        older, previous = previous, current
        if step % refinement == 0:
            history.append(current.reshape(NC, 3)[:16].copy()*np.asarray(SCALE, float))
    return previous.reshape(NC, 3)*np.asarray(SCALE, float), dict(
        refinement=refinement, steps=len(times)-1, linear_solve_residual=worst,
        distinct_factors=len(cache)), np.asarray(history)


def run(cells, workers):
    output = SOURCE.parent/f'gr-time-cause-central-{cells}'
    assert not output.exists(), 'Preserve every diagnostic run.'
    initialize(cells)
    original = json.loads((SOURCE/'plan.json').read_text())
    comparison = json.loads((SOURCE/'time-refinement.json').read_text())
    assert comparison['passed'] is False
    plan = dict(classification='Counterexample candidate', checkpoint=subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], text=True).strip(), cells=cells, core_cells=16,
        workers=workers, fd_scaled_step=str(EPS), fields=FIELDS, refinements=[1, 2, 4, 8, 16, 32],
        assumptions='Affine rho/T/velocity sector about the actual initial data; frozen heat flux and composition perturbations; exact original gravity reconstruction; symmetric derivative of the zero-velocity upwind branch; zero perturbations outside the selected window. This control is not a new full-GR or physical-EOS certificate.',
        diagnostic_scope='Post-failure selected central region. Test reproduction of actual 1/2/4 errors, finer-step response and native EOS/cache agreement. Do not overwrite or relax the original convergence gate.',
        sha256={str(p.relative_to(ROOT)): gr.e.digest(p) for p in
                [Path(__file__), SOURCE/'plan.json', SOURCE/'time-refinement.json', SOURCE/'time-refinement-manifest.json']})
    output.mkdir()
    gr.e.write(output/'plan.json', plan)
    probe = native_probe()
    gr.e.write(output/'native-probe.json', dict(classification='Counterexample candidate', samples=probe))
    print('NATIVE PROBE', json.dumps(dict(exact=sum(x['exact_equal'] for x in probe), samples=len(probe))), flush=True)
    beta = (Z0['P']/Z0['rho']-Z0['aux'][:, 9])/Z0['aux'][:, 10]
    cs2 = Z0['P']/Z0['w']*(Z0['aux'][:, 5]+Z0['aux'][:, 6]*beta)
    speed = gr.e.C*Z0['N']/Z0['a']*np.sqrt(cs2)
    travel = np.cumsum(np.diff(STAR.rf)/speed)
    geometry = dict(sound_time_window_to_core=float(travel[cells-1]-travel[15]),
                    duration_seconds=float(original['duration_seconds']))
    print('CONTROL GEOMETRY', json.dumps(geometry), flush=True)
    b, _ = evaluate(np.zeros(3*cells))
    b = np.asarray(b, float)
    K, M = np.empty((3*cells, 3*cells)), np.empty((3*cells, 3*cells))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for j, k, m in pool.map(column, range(3*cells), chunksize=1):
            K[:, j], M[:, j] = k, m
            if (j+1) % 96 == 0:
                print('CONTROL COLUMNS', j+1, 3*cells, flush=True)
    np.savez_compressed(output/'operator.npz', M=M, K=K, b=b)
    direction = np.sin(np.arange(3*cells)+ld('.3'))*ld('1e-4')
    p, mp = evaluate(direction)
    q, mq = evaluate(-direction)
    checks = []
    for name, matrix, actual in [('spatial', K, (p-q)/2), ('mass', M, (mp-mq)/2)]:
        actual = actual.reshape(cells, 3)
        predicted = (matrix@np.asarray(direction, float)).reshape(cells, 3)
        errors = np.max(abs(predicted-actual), axis=0)/np.maximum(np.max(abs(actual), axis=0), ld('1e-30'))
        checks.append(dict(kind=name, relative=errors.astype(float).tolist()))
        assert np.max(errors) < .005, (name, errors)
    runs, final, histories = [], {}, {}
    for refinement in plan['refinements']:
        final[refinement], record, histories[refinement] = integrate(M, K, b, original, refinement)
        runs.append(record)
        print('CONTROL PATH', json.dumps(record), flush=True)
    actual = {}
    for refinement in [1, 2, 4]:
        with np.load(SOURCE/f'path-{refinement}'/f'step-{39*refinement:04d}.npz') as cp:
            actual[refinement] = np.asarray(cp['delta'][:16, :3], float)
    pairs = []
    for a, c in zip(plan['refinements'][:-1], plan['refinements'][1:]):
        delta = final[c][:16]-final[a][:16]
        record = dict(refinements=[a, c], errors=np.max(abs(delta), axis=0).tolist())
        if c <= 4:
            truth = actual[c]-actual[a]
            record.update(full_GR_errors=np.max(abs(truth), axis=0).tolist(),
                          difference_reproduction_relative=(np.max(abs(delta-truth), axis=0)/np.max(abs(truth), axis=0)).tolist())
        pairs.append(record)
    values, vectors = eig(-K, M)
    coupling = np.linalg.solve(vectors, np.linalg.solve(M, -b))
    finite = np.isfinite(values) & (abs(values) > 1e-12)
    amplitude = np.zeros(len(values))
    amplitude[finite] = abs(vectors[2, finite]*coupling[finite]/values[finite])*float(SCALE[2])
    selected = np.argsort(amplitude)[-12:][::-1]
    modes = [dict(real=float(values[j].real), imag=float(values[j].imag),
                  center_velocity_amplitude=float(amplitude[j])) for j in selected]
    result = dict(classification='Counterexample candidate', geometry=geometry,
        derivative_checks=checks, runs=runs, pairs=pairs, dominant_center_velocity_modes=modes,
        physical_EOS_certified=False, full_GR_time_convergence_fixed=False,
        scope=plan['assumptions'])
    np.savez_compressed(output/'trajectories.npz', **{f'path_{r}': h for r, h in histories.items()})
    gr.e.write(output/'result.json', result)
    gr.e.write(output/'manifest.json', dict(sha256={p.name: gr.e.digest(p) for p in output.iterdir() if p.is_file()}))
    print('CAUSE CONTROL RESULT', json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cells', type=int, default=128)
    parser.add_argument('--workers', type=int, default=15)
    args = parser.parse_args()
    assert 16 < args.cells <= 512 and 1 <= args.workers <= 16
    assert set(os.sched_getaffinity(0)) <= set(range(16))
    run(args.cells, args.workers)
