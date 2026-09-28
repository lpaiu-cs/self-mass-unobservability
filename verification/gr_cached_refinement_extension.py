"""Counterexample candidate: one further native 32/64/128 time comparison.

Keep both earlier failed verdicts. Reference the authenticated complete 32/64
histories, compute only 128 steps, and apply the unchanged five-field gate.
"""
import argparse
import json
import os
from pathlib import Path
from types import FunctionType

import numpy as np
import gr_cached_material_evolution as cached
import gr_cached_time_comparison as replay

OUT = cached.BASE/'pilot-cached-32-64-128'
e = cached.e
context = dict(cached.__dict__, OUT=OUT)
bindings = FunctionType(cached.prior.bindings.__code__, context)
context['bindings'] = bindings
run = FunctionType(cached.prior.run.__code__, context)
completed = FunctionType(cached.prior.completed.__code__, context)


def prepare():
    assert not OUT.exists(), 'Preserve every attempted path.'
    old = cached.bindings()
    replay.check()
    verdict = json.loads((cached.OUT/'time-refinement.json').read_text())
    assert not verdict['passed']
    for name in ('replay-manifest.json', 'time-refinement-manifest.json'):
        for relative, digest in json.loads((cached.OUT/name).read_text())['sha256'].items():
            assert e.digest(e.ROOT/relative) == digest, relative
    # Keep the exact base grid; do not redivide rounded intermediate edges.
    plan = dict(old, phase='one additional 128-step native temporal refinement',
        refinements=[2, 4, 8], new_time_refinement_steps=[32, 64, 128],
        additional_steps_to_compute=128, imported_paths=[
            dict(refinement=2, directory=(cached.PARENT/'path-4').relative_to(e.ROOT).as_posix(),
                 time_plan=(cached.PARENT/'plan.json').relative_to(e.ROOT).as_posix(), source_refinement=4),
            dict(refinement=4, directory=(cached.OUT/'path-4').relative_to(e.ROOT).as_posix(),
                 time_plan=(cached.OUT/'plan.json').relative_to(e.ROOT).as_posix(), source_refinement=4)],
        candidate='Same native v2 space and time equations with the validated exact material cache. Original 32/64 histories stay in place; only the missing 128 trajectory is computed.',
        rationale='The 16/32/64 endpoint differences all decrease and heat-flux orders reach 1.91196, but temperature order 1.10825 fails at outer cell 5733. Test one further halving of time steps without changing the EOS, boundary, initial state, spatial equations, precision or tolerances.',
        stopping_rule='Compare 32/64/128 with the unchanged gate. Preserve failure and stop if it fails. Do not automatically launch a longer physical interval or a further refinement.',
        predecessor_failed_verdict_sha256=e.digest(cached.OUT/'time-refinement.json'))
    plan['bindings'] = dict(old['bindings'])
    for path in (Path(__file__), cached.OUT/'replay-manifest.json',
                 cached.OUT/'time-refinement-manifest.json'):
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    times = cached.parent.wall.prior.time_nodes(plan, 8)
    assert len(times) == 129 and times[-1] == cached.ld(old['coordinate_edges_seconds'][-1])
    assert np.array_equal(times[::2], cached.parent.wall.prior.time_nodes(old, 4))
    assert cached.prior.symbolic()['passed']
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    bindings()
    print('PREPARED ONLY 128 NEW STEPS', flush=True)


def compare():
    assert not (OUT/'time-refinement.json').exists(), 'Preserve a frozen verdict.'
    bindings()
    replay.check()
    paths = [replay.completed(2), replay.completed(4), completed(8)]
    errors = np.array([np.max(abs(b[0][:, :5]-a[0][:, :5]), axis=0)
                       for a, b in zip(paths[:-1], paths[1:])])
    passed, orders = cached.parent.audit.gate(errors)
    result = dict(classification='Counterexample candidate', passed=passed,
        endpoint_maximum_differences=errors.astype(float).tolist(), observed_orders=orders.tolist(),
        paths=[p[1] for p in paths], rigorous_time_error_bound=False,
        physical_EOS_certified=False, physical_exterior_match=False, observational_closure=False)
    e.write(OUT/'time-refinement.json', result)
    files = [Path(__file__), OUT/'plan.json', OUT/'time-refinement.json', OUT/'path-8/manifest.json',
             cached.PARENT/'path-4/manifest.json', cached.OUT/'path-4/manifest.json',
             cached.OUT/'replay-manifest.json']
    e.write(OUT/'time-refinement-manifest.json', dict(sha256={
        p.relative_to(e.ROOT).as_posix(): e.digest(p) for p in files}))
    print('NATIVE 32/64/128 COMPARISON', json.dumps(result), flush=True)
    assert passed, 'The unchanged time convergence gate failed.'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'chain', 'compare'])
    parser.add_argument('--workers', type=int, default=15)
    args = parser.parse_args()
    assert set(os.sched_getaffinity(0)) <= set(range(16)) and 1 <= args.workers <= 16
    if args.command == 'prepare':
        prepare()
    elif args.command == 'chain':
        run(8, args.workers)
        compare()
    else:
        compare()
