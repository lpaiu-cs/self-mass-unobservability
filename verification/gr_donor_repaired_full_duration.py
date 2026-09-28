"""Counterexample candidate: resume the original GR grid after native repair.

Reuse the accepted original 0..41 states and the native repaired step 42.
Keep all physical equations and gates; use a consistent donor pivot only when
the original line search fails, within its original 24 outer iterations.
"""
import argparse
import json
import os
from pathlib import Path
from types import FunctionType, SimpleNamespace

import numpy as np
import gr_compatible_full_duration as full
import gr_consistent_donor_offset as repair

m, e, ld = full.cached, full.e, full.ld
OUT = m.BASE/'production-donor-repaired'
stage = FunctionType(repair.core.stage.__code__, dict(vars(repair.core), finish=repair.finish))
cached = SimpleNamespace(**dict(vars(m), stage=stage))
context = dict(full.context, OUT=OUT, stage=stage)
bindings = FunctionType(m.prior.bindings.__code__, context)
context['bindings'] = bindings
completed = FunctionType(m.prior.completed.__code__, context)
context['completed'] = completed
local = dict(vars(full), OUT=OUT, cached=cached, bindings=bindings, context=context)
check = FunctionType(full.check.__code__, local)
run = FunctionType(full.run.__code__, local)
common_times = FunctionType(full.common_times.__code__, local)
compare = FunctionType(full.compare.__code__, dict(local, common_times=common_times))


def prepare():
    assert not OUT.exists(), 'Preserve every attempted full-duration path.'
    old = full.bindings()
    result = json.loads((repair.OUT/'result.json').read_text())
    assert result['passed'] and result['same_native_equations'] and result['same_time_step']
    assert result['same_tolerances'] and result['final_score'] <= 1
    for rel, digest in json.loads((repair.OUT/'manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    times = full.time_nodes(old, 1)
    OUT.mkdir()
    prefix = OUT/'accepted-prefix-1'
    prefix.mkdir()
    for step in range(43):
        source = (full.OUT/'path-1' if step <= 41 else repair.OUT)/f'step-{step:04d}.npz'
        with np.load(source) as cp:
            assert cp['time_seconds'] == times[step], step
        os.link(source, prefix/source.name)
    rows = [json.loads(s) for s in (full.OUT/'path-1/iterations.jsonl').read_text().splitlines()]
    assert {r['step'] for r in rows} == set(range(1, 43))
    rows += [dict(step=42, **r) for r in result['iterations']]
    for step in range(1, 43):
        records = [r for r in rows if r['step'] == step]
        assert [r['iteration'] for r in records] == list(range(len(records)))
        assert len(records) <= 24 and records[-1]['residual_norm'] <= 1
    (prefix/'iterations.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    e.write(prefix/'manifest.json', dict(classification='Counterexample candidate', accepted_prefix_steps=42,
        source_plan_sha256=e.digest(full.OUT/'plan.json'), repair_manifest_sha256=e.digest(repair.OUT/'manifest.json'),
        sha256={p.name:e.digest(p) for p in prefix.iterdir() if p.is_file()}))
    plan = dict(old, phase='Original full duration after native step-42 donor repair',
        candidate='Same native GR equations and exact 70/140/280 grid with a consistent species/donor affine iteration.',
        imported_prefix='Replay accepted 42/64/128 prefixes with both BDF states and accumulated flux histories. Coarse step 42 is the validated repair of the preserved native failure.',
        additional_steps_to_compute=256, imported_prefix_steps_by_refinement={'1':42, '2':64, '4':128},
        nonlinear_operator='Run the original native stage. If its line search fails, pivot the iteration donor, differentiating species on that same branch while retaining the original native composition anchor and consistent affine residual. At most eight donor pivots per remaining outer iteration, total outer iterations at most 24. Actual native residuals, tolerances, conservation, cones, and the full-time endpoint gate are unchanged.',
        native_equations_unchanged=True, original_time_grid_unchanged=True,
        failed_full_plan=(full.OUT/'plan.json').relative_to(e.ROOT).as_posix(),
        repaired_step_manifest=(repair.OUT/'manifest.json').relative_to(e.ROOT).as_posix(),
        subdivision_states_used=False)
    plan['imported_paths'] = dict(old['imported_paths'])
    plan['imported_paths']['1'] = dict(steps=42, directory=prefix.relative_to(e.ROOT).as_posix(), refinement=1,
        source_plans=[(full.OUT/'plan.json').relative_to(e.ROOT).as_posix(), (repair.OUT/'plan.json').relative_to(e.ROOT).as_posix()])
    plan['bindings'] = dict(old['bindings'])
    files = [Path(__file__), Path(repair.__file__), Path(repair.prior.__file__), Path(repair.core.__file__),
             repair.OUT/'manifest.json', prefix/'manifest.json', full.OUT/'plan.json', full.OUT/'path-1/failure.json']
    for path in files:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    assert plan['exact_time_nodes_seconds'] == old['exact_time_nodes_seconds']
    e.write(OUT/'plan.json', plan)
    check()
    print('NATIVE DONOR REPAIR FULL DURATION PREPARED', plan['new_time_refinement_steps'], 'new steps', 256, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'chain', 'compare'])
    parser.add_argument('--workers', type=int, default=15)
    args = parser.parse_args()
    assert set(os.sched_getaffinity(0)) <= set(range(16)) and 1 <= args.workers <= 16
    if args.command == 'prepare':
        prepare()
    elif args.command == 'chain':
        assert json.loads((OUT/'restart-check.json').read_text())['passed']
        for r in (1, 2, 4):
            run(r, args.workers)
        compare()
    else:
        compare()
