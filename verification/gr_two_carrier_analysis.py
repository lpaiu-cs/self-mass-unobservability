"""Counterexample candidate: finish the full five-field time comparison.

Only complete native paths are compared. The radiative field has the same
finite convergence gate as the temperature and total heat; no endpoint rescue.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import gr_two_carrier_evolution as evolution

e, OUT = evolution.e, evolution.OUT


def gate(errors):
    errors = np.asarray(errors, float)
    assert errors.shape == (2, 5) and np.all(np.isfinite(errors)) and np.all(errors > 0)
    orders = np.log2(errors[0]/errors[1])
    return bool(np.all(errors[1] < errors[0]) and np.min(orders[[1, 3, 4]]) >= 1.5), orders


def prepare():
    target = OUT/'comparison-plan.json'
    assert not target.exists()
    assert not any((OUT/f'path-{r}/result.json').exists() for r in [1, 2, 4])
    assert gate([[4]*5, [1]*5])[0]
    assert not gate([[4]*5, [1, 1, 1, 1, 2]])[0]
    e.write(target, dict(classification='Counterexample candidate',
        sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in [Path(__file__), OUT/'plan.json']},
        refinements=[1, 2, 4], fields=['lnrho_B', 'lnT', 'v/c', 'Qtotal/w0', 'Qrad/w0'],
        gate='Every endpoint difference decreases; lnT,Qtotal,Qrad observed orders at least 1.5. All saved states and all logged internal stages must remain inside the declared cone. Require the full original duration; an early prefix is insufficient.',
        physical_or_rigorous_time_error_certificate=False))


def bindings():
    for name, key in [('comparison-plan.json', 'sha256'), ('plan.json', 'bindings')]:
        for rel, digest in json.loads((OUT/name).read_text())[key].items():
            assert e.digest(e.ROOT/rel) == digest, rel


def completed(refinement):
    bindings()
    folder, plan = OUT/f'path-{refinement}', json.loads((OUT/'plan.json').read_text())
    result = json.loads((folder/'result.json').read_text())
    assert result['completed'] and not (folder/'failure.json').exists()
    manifest = json.loads((folder/'manifest.json').read_text())
    assert manifest['plan_sha256'] == e.digest(OUT/'plan.json')
    for rel, digest in manifest['sha256'].items():
        assert e.digest(folder/rel) == digest, rel
    times = evolution.wall.prior.time_nodes(plan, refinement)
    assert result['steps'] == len(times)-1 and result['duration_seconds'] == float(times[-1])
    stages = {}
    for line in (folder/'iterations.jsonl').read_text().splitlines():
        row = json.loads(line)
        stages.setdefault((row['step'], row['stage']), []).append(row)
    expected = {(i, j) for i in range(1, len(times)) for j in [1, 2]}
    assert set(stages) == expected
    for rows in stages.values():
        assert [r['iteration'] for r in rows] == list(range(len(rows)))
        assert rows[-1]['residual_norm'] <= 1 and len(rows) <= 24
        assert len(rows[-1]['maximum_absolute_residual']) == 31
    internal = [json.loads(line) for line in (folder/'stage-cones.jsonl').read_text().splitlines()]
    assert len(internal) == len(expected) and {(r['step'], r['stage']) for r in internal} == expected
    assert all(r['sampled_cone_inside_light_cone'] for r in internal)
    speed = max(r['maximum_local_rest_characteristic_speed_over_c'] for r in internal)
    star = evolution.initialize(None)
    maximum_defect = 0.
    for step, t in enumerate(times):
        with np.load(folder/f'step-{step:04d}.npz') as saved:
            assert saved['time_seconds'] == t
            delta = saved['delta'].copy()
            if step == 0:
                assert np.all(delta == 0)
                initial_baryons = np.sum(saved['a']*np.exp(star.base[:, 0])*star.volume)
            y = star.base+delta
            star.material_cache = {e.material_key(row): aux for row, aux in zip(
                zip(y[:, 0], y[:, 1], y[:, 5:]), saved['aux'])}
            z = star.state(y)
            for key in ['m', 'mf', 'a', 'N', 'Q', 'aux']:
                assert np.array_equal(z[key], saved[key]), (step, key)
            defect, budget = evolution.energy_budget(star, delta, z, saved['integrated_face_mass_flux'])
            assert np.array_equal(defect, saved['normalized_cell_energy_defect'])
            assert np.array_equal(saved['boundary_exchange'], np.zeros(2))
            assert saved['integrated_face_mass_flux'][0] == saved['integrated_face_mass_flux'][-1] == 0
            cone = evolution.cones(z)
            assert cone['sampled_cone_inside_light_cone'], (step, cone)
            speed = max(speed, cone['maximum_local_rest_characteristic_speed_over_c'])
            maximum_defect = max(maximum_defect, float(abs(defect).max()))
    row = dict(classification='Counterexample candidate', refinement=refinement, steps=len(times)-1,
        duration_seconds=float(times[-1]), maximum_path_speed_over_c=speed,
        maximum_saved_local_energy_defect=maximum_defect,
        maximum_accepted_native_residual=max(r[-1]['residual_norm'] for r in stages.values()),
        maximum_changes=np.max(abs(delta[:, :5]), axis=0).astype(float).tolist(),
        maximum_composition_change=float(abs(delta[:, 5:]).max()),
        relative_baryon_budget_defect=float(np.sum(z['a']*z['D']*star.volume)/initial_baryons-1), **budget)
    return delta, row


def compare():
    target = OUT/'time-refinement.json'
    assert not target.exists(), 'Preserve any failed finite convergence verdict.'
    paths = [completed(r) for r in [1, 2, 4]]
    errors = np.array([np.max(abs(b[0][:, :5]-a[0][:, :5]), axis=0) for a, b in zip(paths[:-1], paths[1:])])
    passed, orders = gate(errors)
    result = dict(classification='Counterexample candidate', passed=passed,
        endpoint_maximum_differences=errors.astype(float).tolist(), observed_orders=orders.tolist(),
        paths=[p[1] for p in paths], rigorous_time_error_bound=False,
        physical_EOS_certified=False, physical_exterior_match=False, observational_closure=False)
    e.write(target, result)
    files = [Path(__file__), target, OUT/'comparison-plan.json', OUT/'plan.json']
    files += [OUT/f'path-{r}/manifest.json' for r in [1, 2, 4]]
    e.write(OUT/'time-refinement-manifest.json', dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files}))
    print('TWO CARRIER FULL TIME COMPARISON', json.dumps(result), flush=True)
    assert passed, 'The frozen five-field finite convergence gate failed.'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'completed', 'compare'])
    parser.add_argument('--refinement', type=int, choices=[1, 2, 4], default=1)
    args = parser.parse_args()
    if args.command == 'completed':
        print(json.dumps(completed(args.refinement)[1], indent=2))
    else:
        globals()[args.command]()
