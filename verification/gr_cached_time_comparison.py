"""Counterexample candidate: replay reused paths on their original time grids.

The cached runner and all executed plans remain immutable. Re-subdividing a
rounded longdouble grid is not bitwise identical to the original subdivision.
Authenticate copied histories, then reuse their original strict native replay.
"""
import argparse
import json
from pathlib import Path
from types import FunctionType

import numpy as np
import gr_cached_material_evolution as cached

OUT = cached.OUT
e = cached.e


def check():
    plan = cached.bindings()
    source_plan = cached.original.bindings()
    assert plan['predeclared_new_gate'] == source_plan['time_refinement_gate']
    nodes = cached.parent.wall.prior.time_nodes
    assert np.array_equal(np.array(plan['coordinate_edges_seconds'], dtype=cached.ld),
                          nodes(source_plan, 2))
    records = []
    for refinement in (1, 2, 4):
        folder = OUT/f'path-{refinement}'
        manifest = json.loads((folder/'manifest.json').read_text())
        assert manifest['plan_sha256'] == e.digest(OUT/'plan.json')
        assert not (folder/'failure.json').exists()
        for name, digest in manifest['sha256'].items():
            assert e.digest(folder/name) == digest, (refinement, name)
        rebuilt = nodes(plan, refinement)
        if refinement in (1, 2):
            source = cached.PARENT/f'path-{2*refinement}'
            origin = json.loads((source/'manifest.json').read_text())
            assert manifest['imported_from'] == source.relative_to(e.ROOT).as_posix()
            assert manifest['source_manifest_sha256'] == e.digest(source/'manifest.json')
            assert manifest['sha256'] == origin['sha256']
            for name, digest in origin['sha256'].items():
                assert e.digest(source/name) == digest, (source, name)
            expected = nodes(source_plan, 2*refinement)
        else:
            expected = rebuilt
        saved = []
        for step in range(len(expected)):
            with np.load(folder/f'step-{step:04d}.npz') as cp:
                saved.append(cp['time_seconds'][()])
        saved = np.array(saved)
        assert np.array_equal(saved, expected), 'Saved timestamps must match the executed plan exactly.'
        assert saved[-1] == rebuilt[-1]
        difference = abs(saved-rebuilt)
        bad = np.flatnonzero(difference)
        records.append(dict(refinement=refinement, steps=len(saved)-1,
            exact_executed_grid=True, redivided_grid_mismatch_steps=bad.tolist(),
            maximum_redivided_time_difference_seconds=str(difference.max()),
            maximum_redivided_ulp_difference=str(np.max(difference[1:]/np.spacing(rebuilt[1:]))),
            maximum_relative_step_size_difference=str(np.max(abs(np.diff(saved)-np.diff(rebuilt))/np.diff(saved)))))
    # Regression: the former comparison rejects the real 32-step copied grid.
    assert records[0]['redivided_grid_mismatch_steps'] == []
    assert records[1]['redivided_grid_mismatch_steps'] == [23, 25, 31]
    assert records[1]['maximum_redivided_ulp_difference'] == '1.0'
    assert records[2]['redivided_grid_mismatch_steps'] == []
    return dict(classification='Counterexample candidate', passed=True, paths=records,
        original_timestamp_equality_preserved=True, original_native_gates_preserved=True,
        physical_model_changed=False,
        cause='The imported 32-step history uses the original direct subdivision. Subdividing its already rounded 16-step intermediate grid differs by one longdouble ULP at three nodes.')


def completed(refinement):
    if refinement == 4:
        result = cached.completed(refinement)
    else:
        assert refinement in (1, 2)
        # check() authenticates every copied byte against this original history.
        delta, report = cached.original.completed(2*refinement)
        result = delta, dict(report, refinement=refinement, source_refinement=2*refinement)
    print('STRICT NATIVE REPLAY COMPLETE', refinement, flush=True)
    return result


def compare():
    assert not (OUT/'time-refinement.json').exists(), 'Preserve an existing scientific verdict.'
    audit = check()
    assert cached.prior.symbolic()['passed']
    e.write(OUT/'replay-check.json', audit)
    print('PROVENANCE TIME GRID CHECK', json.dumps(audit), flush=True)
    # Reuse the unchanged endpoint differences, order calculation and gate.
    compare_native = FunctionType(cached.prior.compare.__code__,
                                  dict(cached.__dict__, completed=completed))
    try:
        compare_native()
    finally:
        if (OUT/'time-refinement-manifest.json').exists():
            files = [Path(__file__), OUT/'replay-check.json', OUT/'time-refinement-manifest.json',
                     cached.PARENT/'plan.json', cached.PARENT/'time-refinement-manifest.json']
            files += [cached.PARENT/f'path-{r}/manifest.json' for r in (2, 4)]
            e.write(OUT/'replay-manifest.json', dict(sha256={
                p.relative_to(e.ROOT).as_posix(): e.digest(p) for p in files}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['check', 'compare'])
    args = parser.parse_args()
    if args.command == 'check':
        print(json.dumps(check(), indent=2))
    else:
        compare()
