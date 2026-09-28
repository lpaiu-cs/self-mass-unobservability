"""Counterexample candidate: resume accepted GR histories after Windows reboot.

No new solver, grid or acceptance gate. Reuse the complete 70-step path and
accepted 134/128 prefixes; recompute only the unfinished 135th step onward.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
from types import FunctionType

import numpy as np
import gr_budget_repaired_full_duration as stopped

m, e, full = stopped.m, stopped.e, stopped.full
OUT = m.BASE/'production-reboot-recovered'
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
    evidence = stopped.OUT/'interruption-20260916.json'
    interrupted = json.loads(evidence.read_text())
    assert interrupted['accepted_steps'] == {'1':70, '2':134, '4':128}
    assert not interrupted['native_failure_file_present']
    for rel, digest in interrupted['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    complete = full.authenticate(stopped.OUT/'path-1')
    assert complete['plan_sha256'] == e.digest(stopped.OUT/'plan.json')
    assert json.loads((stopped.OUT/'path-1/result.json').read_text())['steps'] == 70
    source = stopped.OUT/'path-2'
    assert not any((source/name).exists() for name in ['failure.json', 'result.json', 'step-0135.npz'])
    OUT.mkdir()
    prefix = OUT/'accepted-prefix-2'
    prefix.mkdir()
    times = full.time_nodes(old, 2)
    for step in range(135):
        path = source/f'step-{step:04d}.npz'
        with np.load(path) as saved:
            assert saved['time_seconds'] == times[step]
        os.link(path, prefix/path.name)
    rows = [json.loads(line) for line in (source/'iterations.jsonl').read_text().splitlines()]
    rows = [row for row in rows if row['step'] <= 134]
    assert {row['step'] for row in rows} == set(range(1, 135))
    for step in range(1, 135):
        records = [row for row in rows if row['step'] == step]
        assert [row['iteration'] for row in records] == list(range(len(records)))
        assert len(records) <= 24 and records[-1]['residual_norm'] <= 1
        assert len(records[-1]['maximum_absolute_residual']) == 31
    (prefix/'iterations.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in rows))
    e.write(prefix/'manifest.json', dict(classification='Counterexample candidate', accepted_prefix_steps=134,
        interruption_sha256=e.digest(evidence), sha256={p.name:e.digest(p) for p in prefix.iterdir()}))
    plan = dict(old, phase='Original full duration after Windows upgrade reboot',
        additional_steps_to_compute=158, imported_prefix_steps_by_refinement={'1':70, '2':134, '4':128},
        imported_prefix='Replay the completed 70-step path without numerical integration, then the accepted 134/128 prefixes with exact original clocks and both BDF state/flux histories. Recompute unfinished mid-grid step 135 from accepted step 134.',
        numerical_solver_unchanged_from_budget_repair=True,
        interruption=evidence.relative_to(e.ROOT).as_posix(),
        resource_budget=dict(new_steps_by_refinement={'1':0, '2':6, '4':152}, workers=15,
            cpu_affinity=list(range(16)), blas_threads=1, gpu=False,
            measured_mid_last10_minutes=17.220522256556652,
            assumed_mid_minutes_per_step=[17.2, 60], assumed_fine_minutes_per_step=[15, 20],
            assumed_new_integration_hours=[39.72, 56.67],
            estimate_limit='Fine suffix speed is unmeasured. The interrupted mid step had already taken over 52 minutes; allow up to 60 minutes per remaining mid step in the upper scenario. Replay/final validation and another interruption are excluded.',
            scope='Only the previously authorized 70/140/280 paths and unchanged final gate. No automatic finer grid, longer duration or additional experiment. Preserve and stop on frozen native/budget/cone failures; reassess overruns before any new computation.'))
    plan['imported_paths'] = dict(old['imported_paths'])
    for r, count, folder in [(1,70,stopped.OUT/'path-1'), (2,134,prefix)]:
        plan['imported_paths'][str(r)] = dict(steps=count, directory=folder.relative_to(e.ROOT).as_posix(),
            refinement=r, source_plan=(stopped.OUT/'plan.json').relative_to(e.ROOT).as_posix())
    plan['bindings'] = dict(old['bindings'])
    for p in [Path(__file__), Path(__file__).with_name('run_gr_reboot_recovery.ps1'), evidence,
              prefix/'manifest.json', stopped.OUT/'path-1/manifest.json', stopped.OUT/'plan.json']:
        plan['bindings'][p.relative_to(e.ROOT).as_posix()] = e.digest(p)
    assert plan['exact_time_nodes_seconds'] == old['exact_time_nodes_seconds']
    assert sum(len(full.time_nodes(plan,r))-1-plan['imported_paths'][str(r)]['steps'] for r in (1,2,4)) == 158
    assert local['cached'].stage is stopped.repair.stage
    e.write(OUT/'plan.json', plan)
    check()
    print('REBOOT RECOVERY PREPARED; ORIGINAL GRID, 158 NEW STEPS', flush=True)


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
                command=['verification/gr_reboot_recovery.py', 'chain', '--workers', str(args.workers)],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                start_ticks=(Path('/proc')/str(pid)/'stat').read_text().rsplit(') ',1)[1].split()[19],
                plan_sha256=e.digest(OUT/'plan.json'), cpu_affinity=sorted(os.sched_getaffinity(0))), stream, indent=2)
        for refinement in (1, 2, 4):
            run(refinement, args.workers)
        compare()
    else:
        compare()
