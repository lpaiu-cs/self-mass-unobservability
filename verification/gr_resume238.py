"""Counterexample candidate: finish the original GR grid after WSL termination.

Replay completed 70/140 paths and the accepted 238-step fine prefix. The native
solver, exact clocks, iteration budget and final convergence gate are unchanged.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
from types import FunctionType

import numpy as np
import gr_reboot_recovery as stopped

m, e, full = stopped.m, stopped.e, stopped.full
OUT = m.BASE/'production-resume238'
context = dict(stopped.context, OUT=OUT)
bindings = FunctionType(m.prior.bindings.__code__, context)
context['bindings'] = bindings
completed = FunctionType(m.prior.completed.__code__, context)
context['completed'] = completed
local = dict(stopped.local, OUT=OUT, bindings=bindings, context=context)
check = FunctionType(full.check.__code__, local)
run = FunctionType(full.run.__code__, local)
common_times = FunctionType(full.common_times.__code__, local)
compare = FunctionType(full.compare.__code__, dict(local, common_times=common_times))


def prepare():
    assert not OUT.exists(), 'Preserve every attempted run.'
    old = stopped.bindings()
    evidence = stopped.OUT/'interruption-20260917.json'
    interrupted = json.loads(evidence.read_text())
    assert interrupted['accepted_steps'] == {'1':70, '2':140, '4':238}
    assert not interrupted['native_failure_file_present']
    for rel, digest in interrupted['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    for refinement, count in [(1,70), (2,140)]:
        folder = stopped.OUT/f'path-{refinement}'
        manifest = full.authenticate(folder)
        result = json.loads((folder/'result.json').read_text())
        assert manifest['plan_sha256'] == e.digest(stopped.OUT/'plan.json')
        assert result['completed'] and result['steps'] == count
    source = stopped.OUT/'path-4'
    assert not any((source/name).exists() for name in ['failure.json', 'result.json', 'step-0239.npz'])
    OUT.mkdir()
    prefix = OUT/'accepted-prefix-4'
    prefix.mkdir()
    times = full.time_nodes(old, 4)
    for step in range(239):
        path = source/f'step-{step:04d}.npz'
        with np.load(path) as saved:
            assert saved['time_seconds'] == times[step]
        os.link(path, prefix/path.name)
    rows = [json.loads(line) for line in (source/'iterations.jsonl').read_text().splitlines()]
    rows = [row for row in rows if row['step'] <= 238]
    assert {row['step'] for row in rows} == set(range(1, 239))
    for step in range(1, 239):
        records = [row for row in rows if row['step'] == step]
        assert [row['iteration'] for row in records] == list(range(len(records)))
        assert len(records) <= 24 and records[-1]['residual_norm'] <= 1
        assert len(records[-1]['maximum_absolute_residual']) == 31
    (prefix/'iterations.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in rows))
    e.write(prefix/'manifest.json', dict(classification='Counterexample candidate', accepted_prefix_steps=238,
        interruption_sha256=e.digest(evidence), sha256={p.name:e.digest(p) for p in prefix.iterdir()}))
    plan = dict(old, phase='Finish original GR interval after WSL lifetime interruption',
        additional_steps_to_compute=42, imported_prefix_steps_by_refinement={'1':70, '2':140, '4':238},
        imported_prefix='Replay completed 70/140 paths without new integration and accepted 238 fine steps with exact original clocks and both BDF state/flux histories. Compute only original fine steps 239 through 280, then the unchanged final comparison.',
        numerical_solver_unchanged_from_budget_repair=True,
        interruption=evidence.relative_to(e.ROOT).as_posix(),
        resource_budget=dict(new_steps_by_refinement={'1':0, '2':0, '4':42}, workers=15,
            cpu_affinity=list(range(16)), blas_threads=1, gpu=False,
            measured_fine_minutes=interrupted['measured_fine_minutes'],
            assumed_fine_minutes_per_step=[10, 15], assumed_new_integration_hours=[7, 10.5],
            estimate_limit='Recent 10-step mean is about 10.02 minutes. Later nonlinear costs, replay/final validation and another interruption are not bounded by this scenario estimate.',
            scope='Only the original 70/140/280 paths and unchanged final gate. No automatic finer grid, longer duration or additional experiment. Preserve and stop on frozen native/budget/cone failures; reassess overruns before new computation.'))
    plan['imported_paths'] = dict(old['imported_paths'])
    for refinement, count, folder in [(1,70,stopped.OUT/'path-1'), (2,140,stopped.OUT/'path-2'), (4,238,prefix)]:
        plan['imported_paths'][str(refinement)] = dict(steps=count,
            directory=folder.relative_to(e.ROOT).as_posix(), refinement=refinement,
            source_plan=(stopped.OUT/'plan.json').relative_to(e.ROOT).as_posix())
    plan['bindings'] = dict(old['bindings'])
    for path in [Path(__file__), Path(__file__).with_name('run_gr_resume238.ps1'), evidence,
                 prefix/'manifest.json', stopped.OUT/'path-1/manifest.json',
                 stopped.OUT/'path-2/manifest.json', stopped.OUT/'plan.json']:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    assert plan['exact_time_nodes_seconds'] == old['exact_time_nodes_seconds']
    assert sum(len(full.time_nodes(plan,r))-1-plan['imported_paths'][str(r)]['steps'] for r in (1,2,4)) == 42
    assert local['cached'].stage is stopped.local['cached'].stage
    e.write(OUT/'plan.json', plan)
    check()
    print('RECOVERY PREPARED; ORIGINAL GRID, 42 NEW STEPS', flush=True)


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
        assert not (OUT/'path-1').exists() and (OUT/'windows-launcher.json').exists()
        bindings()
        pid = os.getpid()
        with (OUT/'launch.json').open('x') as stream:
            json.dump(dict(pid=pid, started_utc=datetime.now(timezone.utc).isoformat(), cwd=str(e.ROOT),
                command=['verification/gr_resume238.py', 'chain', '--workers', str(args.workers)],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                start_ticks=(Path('/proc')/str(pid)/'stat').read_text().rsplit(') ',1)[1].split()[19],
                plan_sha256=e.digest(OUT/'plan.json'), cpu_affinity=sorted(os.sched_getaffinity(0))), stream, indent=2)
        for refinement in (1, 2, 4):
            run(refinement, args.workers)
        compare()
    else:
        compare()
