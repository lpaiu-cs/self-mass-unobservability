"""Counterexample candidate: continue the passing pilot to the full GR interval.

Replay each complete pilot prefix on its actual time grid, including both BDF
histories and accumulated fluxes. Compute only the suffix with unchanged native
equations, exact material cache, residual/budget/cone gates and endpoint gate.
"""
import argparse
import json
import os
from pathlib import Path
from types import FunctionType, SimpleNamespace

import numpy as np
import gr_cached_refinement_extension as pilot
from gr_execution_acceleration import state as load_state

cached, e = pilot.cached, pilot.e
ld = cached.ld
OUT = cached.BASE/'production-cached-prefix'
STATE_KEYS = ['m', 'mf', 'a', 'N', 'Q', 'aux', 'dU']


def time_nodes(plan, refinement):
    return np.array(plan['exact_time_nodes_seconds'][str(refinement)], dtype=ld)


# Only the time-grid provider changes in the runner/replay namespace. The
# physical residual and native stage still use the frozen v2 implementation.
parent = SimpleNamespace(**vars(cached.parent))
parent.wall = SimpleNamespace(**vars(parent.wall))
parent.wall.prior = SimpleNamespace(**dict(vars(parent.wall.prior), time_nodes=time_nodes))
context = dict(vars(cached), OUT=OUT, parent=parent)
bindings = FunctionType(cached.prior.bindings.__code__, context)
context['bindings'] = bindings
completed = FunctionType(cached.prior.completed.__code__, context)
context['completed'] = completed


def authenticate(folder):
    manifest = json.loads((folder/'manifest.json').read_text())
    assert not (folder/'failure.json').exists()
    for name, digest in manifest['sha256'].items():
        assert e.digest(folder/name) == digest, (folder, name)
    return manifest


def same_arrays(saved, values):
    assert set(saved.files) == set(values)
    for key, value in values.items():
        assert np.array_equal(saved[key], value), key


def restored(star, folder, step, expected_time):
    star.material_cache.clear()
    delta, state = load_state(star, folder, step)
    with np.load(folder/f'step-{step:04d}.npz') as cp:
        assert cp['time_seconds'] == expected_time
        for key in STATE_KEYS:
            assert np.array_equal(state[key], cp[key]), (step, key)
    return delta, state


def prepare():
    assert not OUT.exists(), 'Preserve every attempted full-time experiment.'
    old = pilot.bindings()
    verdict = json.loads((pilot.OUT/'time-refinement.json').read_text())
    assert verdict['passed'] and [p['steps'] for p in verdict['paths']] == [32, 64, 128]
    for rel, digest in json.loads((pilot.OUT/'time-refinement-manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    full = json.loads((cached.original.SOURCE/'plan.json').read_text())
    stop = ld(full['coordinate_edges_seconds'][-1])
    join = ld(old['coordinate_edges_seconds'][-1])
    cap = np.diff(np.array(full['coordinate_edges_seconds'], dtype=ld)).max()
    h, suffix = join/32, [join]
    while suffix[-1] < stop:
        h = min(2*h, cap, stop-suffix[-1])
        suffix.append(suffix[-1]+h)
    paths = [(1, 32, cached.PARENT/'path-4', cached.PARENT/'plan.json', 4),
             (2, 64, cached.OUT/'path-4', cached.OUT/'plan.json', 4),
             (4, 128, pilot.OUT/'path-8', pilot.OUT/'plan.json', 8)]
    grids, imports, files = {}, {}, []
    for r, count, folder, source_plan, source_r in paths:
        manifest = authenticate(folder)
        assert manifest['plan_sha256'] == e.digest(source_plan)
        prefix = cached.parent.wall.prior.time_nodes(json.loads(source_plan.read_text()), source_r)
        assert len(prefix) == count+1 and prefix[-1] == join
        for step, t in enumerate(prefix):
            with np.load(folder/f'step-{step:04d}.npz') as cp:
                assert cp['time_seconds'] == t
        tail = np.concatenate([np.linspace(a, b, r+1, dtype=ld)[1:]
                               for a, b in zip(suffix[:-1], suffix[1:])])
        times = np.r_[prefix, tail]
        assert times[-1] == stop and np.all(np.diff(times) > 0)
        assert np.max(np.diff(times)[1:]/np.diff(times)[:-1]) <= ld('2.000000000001')
        grids[str(r)] = list(map(str, times))
        imports[str(r)] = dict(steps=count, directory=folder.relative_to(e.ROOT).as_posix(),
                              plan=source_plan.relative_to(e.ROOT).as_posix(), refinement=source_r)
        files.extend([folder/'manifest.json', source_plan])
    plan = dict(old, phase='full original duration after passing native pilot',
        candidate='Frozen compatible v2 native GR with exact material cache and authenticated same-equation pilot prefixes.',
        duration_seconds=float(stop), refinements=[1, 2, 4],
        coordinate_edges_seconds=grids['1'], exact_time_nodes_seconds=grids,
        imported_paths=imports, imported_prefix_steps_by_refinement={r:v['steps'] for r,v in imports.items()},
        new_time_refinement_steps=[len(grids[str(r)])-1 for r in (1, 2, 4)],
        additional_steps_to_compute=sum(len(grids[r])-1-imports[r]['steps'] for r in imports),
        imported_prefix='Replay each complete 32/64/128 pilot from the same original initial state. Preserve both BDF state histories, prior step sizes and integrated energy/baryon flux histories. No backward-Euler restart at the join.',
        reuse_scope='Only the same-equation passing pilot histories are accepted. Original different-spatial-operator 39/78/156 states are not imported.',
        continuation_grid='Keep each executed pilot grid exactly. After the common pilot endpoint, double coarse step sizes up to the original full-grid maximum; split each new coarse interval by 1/2/4. Freeze all nodes before computing the suffix.',
        maximum_coarse_step_seconds=str(cap), suffix_base_steps=len(suffix)-1,
        native_pilot_passed=True, native_pilot_manifest=(pilot.OUT/'time-refinement-manifest.json').relative_to(e.ROOT).as_posix(),
        stopping_rule='Run all three full paths and strict native replay, then the unchanged five-field endpoint convergence gate. Save every accepted step; preserve and stop on a native, budget, cone or final convergence failure.',
        common_time_scope='Report all coarsest common time indices, whole-domain and core/outer differences. A few historical prefix timestamps differ by one longdouble ULP; replay uses each actual timestamp. Endpoint equality is exact. The frozen endpoint gate is unchanged.')
    for key in ['imported_prefix_steps', 'imported_path_validation', 'predecessor_failed_verdict_sha256']:
        plan.pop(key, None)
    plan['bindings'] = dict(old['bindings'])
    files += [Path(__file__), Path(load_state.__code__.co_filename),
              pilot.OUT/'time-refinement-manifest.json', cached.original.SOURCE/'plan.json']
    for path in files:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    bindings()
    print('FULL DURATION PREPARED', plan['new_time_refinement_steps'],
          'new steps', plan['additional_steps_to_compute'], flush=True)


def check():
    plan = bindings()
    assert cached.prior.symbolic()['passed']
    rows = []
    for r in (1, 2, 4):
        info = plan['imported_paths'][str(r)]
        folder, n = e.ROOT/info['directory'], info['steps']
        authenticate(folder)
        star, times = cached.initialize(None), time_nodes(plan, r)
        older, previous, current = [restored(star, folder, k, times[k]) for k in (n-2, n-1, n)]
        h = times[n]-times[n-1]
        coefficients = cached.weights(h, times[n-1]-times[n-2])
        native, z = cached.residual(star, current[0], previous, older, h, coefficients)
        score = float(np.max(abs(native)/cached.ATOL))
        assert score <= 1
        with np.load(folder/f'step-{n:04d}.npz') as cp, \
             np.load(folder/f'step-{n-1:04d}.npz') as prev, \
             np.load(folder/f'step-{n-2:04d}.npz') as before:
            c0, c1, c2 = coefficients
            area = 4*np.pi*star.rf**2
            for key, flux in [('integrated_energy_flux', 1), ('integrated_baryon_flux', 0)]:
                rebuilt = (-c1*prev[key]-c2*before[key]+h*e.C*area*z['fluxes'][flux])/c0
                assert np.array_equal(rebuilt, cp[key]), (r, key)
            values = {key:cp[key].copy() for key in cp.files}
            same_arrays(cp, values)
            values['integrated_energy_flux'][1] = np.nextafter(values['integrated_energy_flux'][1], ld(np.inf))
            try:
                same_arrays(cp, values)
            except AssertionError:
                pass
            else:
                raise AssertionError('Altered flux history was accepted')
        rows.append(dict(refinement=r, prefix_steps=n, native_score=score,
                         restored_flux_history_exact=True, altered_history_rejected=True))
    result = dict(classification='Counterexample candidate', passed=True, paths=rows,
                  symbolic_passed=True, scope='Exact end-of-prefix BDF state and flux restoration; not full-duration convergence.')
    e.write(OUT/'restart-check.json', result)
    print('FULL DURATION RESTORE CHECK', json.dumps(result), flush=True)


def run(refinement, workers):
    plan = bindings()
    info = plan['imported_paths'][str(refinement)]
    folder, count = e.ROOT/info['directory'], info['steps']
    authenticate(folder)
    times = time_nodes(plan, refinement)
    records = {}
    for line in (folder/'iterations.jsonl').read_text().splitlines():
        row = json.loads(line)
        records.setdefault(row['step'], []).append(row)

    def replay_stage(star, previous, older, h, coefficients, log):
        step = star.operator_step+1
        if step > count:
            return cached.stage(star, previous, older, h, coefficients, log)
        star.operator_step = step
        star.operator_time += h
        assert star.operator_time == times[step]
        delta, state = restored(star, folder, step, times[step])
        native, state = cached.residual(star, delta, previous, older, h, coefficients)
        assert np.max(abs(native)/cached.ATOL) <= 1, step
        for row in records[step]:
            log({key:value for key,value in row.items() if key != 'step'})
        return delta, state

    def replay_save(path, delta, state, t, energy, baryon, energy_defect, baryon_defect):
        step = int(path.stem.split('-')[-1]) if path.stem.startswith('step-') else None
        if step is None or step > count:
            return cached.save_state(path, delta, state, t, energy, baryon, energy_defect, baryon_defect)
        source = folder/path.name
        with np.load(source) as cp:
            same_arrays(cp, dict(delta=delta, **{k:state[k] for k in STATE_KEYS}, time_seconds=t,
                integrated_energy_flux=energy, integrated_baryon_flux=baryon,
                normalized_energy_defect=energy_defect, normalized_baryon_defect=baryon_defect))
        assert not path.exists()
        os.link(source, path)

    execute = FunctionType(cached.prior.run.__code__,
                           dict(context, stage=replay_stage, save_state=replay_save))
    execute(refinement, workers)


def common_times():
    plan = bindings()
    grids = {r:time_nodes(plan, r) for r in (1, 2, 4)}
    rows = []
    for step, t in enumerate(grids[1]):
        states = []
        for r in (1, 2, 4):
            with np.load(OUT/f'path-{r}/step-{step*r:04d}.npz') as cp:
                assert cp['time_seconds'] == grids[r][step*r]
                states.append(cp['delta'][:, :5])
        errors = [abs(b-a) for a,b in zip(states[:-1], states[1:])]
        rows.append(dict(step=step, time_seconds=str(t),
            timestamp_difference_seconds=[str(grids[r][step*r]-t) for r in (2, 4)],
            maximum_differences=[d.max(axis=0).astype(float).tolist() for d in errors],
            maximum_cells=[d.argmax(axis=0).tolist() for d in errors],
            core16_maximum_differences=[d[:16].max(axis=0).astype(float).tolist() for d in errors],
            outer32_maximum_differences=[d[-32:].max(axis=0).astype(float).tolist() for d in errors]))
    e.write(OUT/'common-time-comparison.json', dict(classification='Counterexample candidate',
        endpoint_gate_unchanged=True, rows=rows, rigorous_error_bound=False))


def compare():
    try:
        FunctionType(cached.prior.compare.__code__, context)()
    finally:
        if (OUT/'time-refinement-manifest.json').exists():
            common_times()
            files = [OUT/'time-refinement-manifest.json', OUT/'common-time-comparison.json', OUT/'restart-check.json']
            e.write(OUT/'execution-manifest.json', dict(sha256={
                p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'check', 'run', 'chain', 'compare'])
    parser.add_argument('--refinement', type=int, choices=[1, 2, 4], default=1)
    parser.add_argument('--workers', type=int, default=15)
    args = parser.parse_args()
    assert set(os.sched_getaffinity(0)) <= set(range(16)) and 1 <= args.workers <= 16
    if args.command == 'prepare':
        prepare()
        check()
    elif args.command == 'chain':
        assert json.loads((OUT/'restart-check.json').read_text())['passed']
        for r in (1, 2, 4):
            run(r, args.workers)
        compare()
    elif args.command == 'run':
        run(args.refinement, args.workers)
    else:
        globals()[args.command]()
