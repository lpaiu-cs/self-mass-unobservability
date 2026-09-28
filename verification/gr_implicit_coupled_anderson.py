"""Native-residual acceleration and continuation of the frozen SDIRK paths.

Counterexample candidate. The physical RHS, time nodes and acceptance tolerance
are unchanged. Completed parent steps are copied and replayed, never refitted.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import subprocess
import time

import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import splu

import gr_implicit_coupled_evolution as old

e, ld = old.e, old.ld
OUT = e.g.OUT/'gr-implicit-coupled-anderson'


def direction(x, f, factor, history):
    raw = factor.solve(-np.asarray(f, float))
    if not history:
        return raw, raw
    pairs = history[-6:]+[(x, f)]
    dx = np.column_stack([b[0]-a[0] for a, b in zip(pairs[:-1], pairs[1:])])
    df = np.column_stack([b[1]-a[1] for a, b in zip(pairs[:-1], pairs[1:])])
    # Apply the CURRENT inverse to every secant, including across refreshes.
    w = factor.solve(np.asarray(df, float))
    coefficients = np.linalg.lstsq(w, -raw, rcond=1e-8)[0]
    mixed = raw-(dx-w)@coefficients
    return np.asarray(mixed, float), raw


def stage(star, stage_base, seed, h, log):
    delta, history, factor = seed.copy(), [], None
    rates, z = star.rhs(delta)
    evaluations, trials = 1, []
    for iteration in range(24):
        residual = delta-stage_base-h*rates
        tolerance = old.ATOL+ld('1e-7')*np.maximum(abs(delta), abs(stage_base))
        norm = float(np.max(abs(residual)/tolerance))
        log(dict(iteration=iteration, residual_norm=norm, native_evaluations=evaluations,
            maximum_absolute_residual=np.max(abs(residual), axis=0).astype(float).tolist(),
            preceding_trials=trials))
        if norm <= 1:
            return delta, rates, z
        if iteration == 23:
            break
        if factor is None or iteration % 4 == 0:
            J = old.jacobian(old.tangent(star, delta, z), delta)
            factor = splu(eye(4*star.n, format='csc')-float(h)*J)
        x = (delta[:, :4]/old.SCALE).ravel().copy()
        f = (residual[:, :4]/old.SCALE).ravel().copy()
        mixed, raw = direction(x, f, factor, history)
        trials, accepted = [], None
        proposals = [mixed, raw] if history else [raw]
        for proposal in proposals:
            correction = proposal.reshape(star.n, 4).astype(ld)*old.SCALE
            fraction = min(1., .1/max(float(np.max(abs(correction[:, :2]))), 1e-300),
                           .01/max(float(np.max(abs(correction[:, 2]))), 1e-300))
            for backtrack in range(8):
                candidate = delta.copy()
                candidate[:, :4] += fraction*correction
                star.current = z
                candidate[:, 4:] = old.species(star, stage_base, h, star.base+candidate)
                evaluations += 1
                try:
                    candidate_rates, candidate_z = star.rhs(candidate)
                    candidate_residual = candidate-stage_base-h*candidate_rates
                    candidate_tolerance = old.ATOL+ld('1e-7')*np.maximum(abs(candidate), abs(stage_base))
                    score = float(np.max(abs(candidate_residual)/candidate_tolerance))
                    trials.append(dict(fraction=fraction, norm=score))
                    if score <= 1 or score < norm*(1-1e-4*fraction):
                        accepted = candidate, candidate_rates, candidate_z
                        break
                except (AssertionError, np.linalg.LinAlgError, FloatingPointError) as error:
                    trials.append(dict(fraction=fraction, invalid_trial=repr(error)))
                fraction /= 2
            if accepted is not None:
                break
        if accepted is None:
            raise RuntimeError(('Native residual line search failed', norm, trials))
        history = (history+[(x, f)])[-6:]
        delta, rates, z = accepted
    raise RuntimeError(('Native accelerated implicit stage did not converge', norm, trials))


def check():
    # Secants of a coupled linear residual must recover its exact Newton step
    # when the history spans the space, even with a different current inverse.
    from scipy.linalg import lu_factor, lu_solve
    class Factor:
        def __init__(self, matrix):
            self.factor = lu_factor(matrix)
        def solve(self, rhs):
            return lu_solve(self.factor, rhs)
    matrix = np.array([[3., 2.], [-1., 4.]])
    inverse = Factor(np.diag([2., 7.]))
    target = np.array([2., -3.])
    points = [np.array([0., 0.]), np.array([1., 0.]), np.array([1., 1.])]
    pairs = [(x, matrix@x-target) for x in points]
    step, _ = direction(*pairs[-1], inverse, pairs[:-1])
    assert np.max(abs(matrix@(points[-1]+step)-target)) < 1e-12
    # Exercise the actual stage/line-search code on a nonlinear stiff residual.
    # This is a software control, never a native EOS or stellar validation.
    from types import FunctionType, SimpleNamespace
    from scipy.sparse import diags
    class Toy:
        n = 1
        base = np.zeros((1, 30), dtype=ld)
        def rhs(self, value):
            rates = np.zeros_like(value)
            rates[:, :4] = -20*value[:, :4]-1e6*value[:, :4]**3
            self.current = {}
            return rates, self.current
    surrogate = SimpleNamespace(ATOL=np.full(30, 1e-14, dtype=ld), SCALE=np.ones(4, dtype=ld),
        tangent=lambda *args: None, jacobian=lambda *args: diags(np.full(4, -.1)),
        species=lambda *args: np.zeros((1, 26), dtype=ld))
    tested = FunctionType(stage.__code__, dict(globals(), old=surrogate))
    base = np.zeros((1, 30), dtype=ld)
    base[0, :4] = [.2, .1, .001, .0001]
    logs = []
    value, rates, _ = tested(Toy(), base, np.zeros_like(base), ld(1), logs.append)
    assert np.max(abs(value-base-rates)/(surrogate.ATOL+ld('1e-7')*np.maximum(abs(value), abs(base)))) <= 1
    assert any(t.get('fraction', 1) < 1 for row in logs for t in row['preceding_trials'])
    print('PASS current-inverse multisecant root and nonlinear stage/backtracking controls', len(logs), flush=True)


def prepare(refinement):
    parent = old.OUT/f'path-{refinement}'
    assert not (parent/'result.json').exists()
    failed = (parent/'failure.json').exists()
    receipt = json.loads((parent/('failure.json' if failed else 'accelerator-stop.json')).read_text())
    assert (receipt['completed'] is False) if failed else (receipt['stopped'] is True)
    for process in Path('/proc').glob('[0-9]*'):
        try:
            command = (process/'cmdline').read_bytes().replace(bytes([0]), b' ').decode()
        except OSError:
            continue
        assert f'python3 verification/gr_implicit_coupled_evolution.py run --refinement {refinement} ' not in command, 'Parent producer remains live.'
    plan = json.loads((old.OUT/'plan.json').read_text())
    for rel, digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    folder = OUT/f'path-{refinement}'
    assert not folder.exists()
    if not OUT.exists():
        OUT.mkdir()
        plan.update(maximum_stage_iterations=24,
            parent_plan_sha256=e.digest(old.OUT/'plan.json'),
            nonlinear_solver='Six secants of actual native residuals, preconditioned by the current sparse inverse; backtrack against the unchanged residual norm. Preserve failed or explicitly stopped 16-iteration parent attempts with their actual status. No time-node or physical-RHS change.',
            checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip())
        plan['bindings'][Path(__file__).relative_to(e.ROOT).as_posix()] = e.digest(Path(__file__))
        e.write(OUT/'plan.json', plan)
    else:
        assert json.loads((OUT/'plan.json').read_text())['bindings'][Path(__file__).relative_to(e.ROOT).as_posix()] == e.digest(Path(__file__))
    folder.mkdir()
    progress = json.loads((parent/'progress.json').read_text())
    start = progress['step']
    assert start == (receipt['step']-1 if failed else receipt['last_completed_step'])
    files = [p for p in parent.iterdir() if p.is_file()]
    e.write(folder/'resume.json', dict(parent=parent.relative_to(e.ROOT).as_posix(), start=start,
        parent_files_sha256={p.name:e.digest(p) for p in files},
        rule='Copy completed steps only; replay the last checkpoint with freshly evaluated native EOS before continuing. The uncompleted parent stage and its failed/stopped status remain in the original parent record.'))
    for step in range(start+1):
        name = f'step-{step:04d}.npz'
        (folder/name).write_bytes((parent/name).read_bytes())
    accepted_lines = [line for line in (parent/'iterations.jsonl').read_text().splitlines()
                      if json.loads(line)['step'] <= start]
    (folder/'iterations.jsonl').write_text('\n'.join(accepted_lines)+('\n' if accepted_lines else ''))
    print('PREPARED ACCELERATED CONTINUATION', refinement, start, flush=True)


def run(refinement, workers):
    plan = json.loads((OUT/'plan.json').read_text())
    folder = OUT/f'path-{refinement}'
    assert not (folder/'progress.json').exists() and not (folder/'failure.json').exists()
    for rel, digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    for rel, digest in plan['runtime_sha256'].items():
        assert e.digest(Path(rel)) == digest, rel
    resume = json.loads((folder/'resume.json').read_text())
    parent = e.ROOT/resume['parent']
    for rel, digest in resume['parent_files_sha256'].items():
        assert e.digest(parent/rel) == digest, rel
    start = resume['start']
    checkpoint = np.load(folder/f'step-{start:04d}.npz')
    edges = np.array(plan['time_edges_tau'], dtype=ld)*e.TAU
    times = np.r_[ld(0), np.concatenate([np.linspace(a, b, refinement+1, dtype=ld)[1:]
                                        for a, b in zip(edges[:-1], edges[1:])])]
    assert times[start] == checkpoint['time_seconds']
    began = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers, initializer=e.worker_init) as pool:
        star = old.initialize(pool)
        star.material_cache = {}  # Fresh native replay, not a cached assertion.
        delta = checkpoint['delta'].copy()
        exchange = checkpoint['boundary_exchange'].copy()
        mass_flux = checkpoint['integrated_face_mass_flux'].copy()
        for step in range(start, len(times)):
            _, z = star.rhs(delta)
            cell_defect, budget = old.energy_budget(star, delta, z, mass_flux)
            if step == start:
                for key in ['m', 'mf', 'a', 'N', 'Q', 'aux']:
                    assert np.array_equal(z[key], checkpoint[key]), ('Fresh native checkpoint replay', key)
                assert np.array_equal(cell_defect, checkpoint['normalized_cell_energy_defect'])
                e.write(folder/'checkpoint-replay.json', dict(passed=True, fresh_native=True, step=start))
            else:
                np.savez_compressed(folder/f'step-{step:04d}.npz', delta=delta,
                    **{k:z[k] for k in ['m', 'mf', 'a', 'N', 'Q', 'aux']}, time_seconds=times[step],
                    boundary_exchange=exchange, integrated_face_mass_flux=mass_flux,
                    normalized_cell_energy_defect=cell_defect)
            row = dict(classification='Counterexample candidate', step=step, time_seconds=float(times[step]),
                maximum_changes=np.max(abs(delta[:, :4]), axis=0).astype(float).tolist(),
                elapsed_seconds=time.monotonic()-began, **budget)
            e.write(folder/'progress.json', row)
            print('ACCELERATED COUPLED STEP', refinement, json.dumps(row), flush=True)
            if step == len(times)-1:
                break
            dt = times[step+1]-times[step]
            def log(record):
                with (folder/'iterations.jsonl').open('a') as stream:
                    stream.write(json.dumps(dict(step=step+1, stage=stage_number, **record))+'\n')
                print('ACCELERATED NATIVE STAGE', refinement, step+1, stage_number,
                      record['iteration'], record['residual_norm'], record['native_evaluations'], flush=True)
            try:
                stage_number = 1
                first, f1, z1 = stage(star, delta, delta, old.GAMMA*dt, log)
                stage_number = 2
                delta, _, z2 = stage(star, delta+(1-old.GAMMA)*dt*f1, first, old.GAMMA*dt, log)
            except Exception as error:
                e.write(folder/'failure.json', dict(classification='Counterexample candidate', completed=False,
                    step=step+1, stage=stage_number, reason=repr(error)))
                raise
            exchange += dt*((1-old.GAMMA)*z1['boundary_rates']+old.GAMMA*z2['boundary_rates'])
            flux = ((1-old.GAMMA)*star.faces(z1['N']*z1['S']/z1['a'], odd=True)
                    +old.GAMMA*star.faces(z2['N']*z2['S']/z2['a'], odd=True))
            mass_flux += dt*e.GRAV*e.C*4*np.pi*star.rf**2*flux
        e.write(folder/'result.json', dict(classification='Counterexample candidate', completed=True,
            steps=len(times)-1, duration_seconds=float(times[-1]), actual_native_nonlinear_path=True,
            physical_EOS_certified=False, physical_exterior_match=False, observational_closure=False))
    e.write(folder/'manifest.json', dict(plan_sha256=e.digest(OUT/'plan.json'),
        sha256={p.name:e.digest(p) for p in folder.iterdir() if p.is_file()}))
    assert not (parent/'continuation.json').exists()
    e.write(parent/'continuation.json', dict(path=folder.relative_to(e.ROOT).as_posix(),
        manifest_sha256=e.digest(folder/'manifest.json')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['check', 'prepare', 'run'])
    parser.add_argument('--refinement', type=int, choices=[1, 2, 4], default=1)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.command == 'check':
        check()
    elif args.command == 'prepare':
        prepare(args.refinement)
    else:
        run(args.refinement, args.workers)
