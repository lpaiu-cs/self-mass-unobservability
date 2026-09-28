"""Counterexample candidate: replay the accepted prefix after WSL termination.

The validated donor solver, original clocks and scientific gates are unchanged.
Only accepted states are imported; the unfinished step is recomputed.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
from types import FunctionType

import numpy as np
import gr_donor_repaired_full_duration as stopped

m, e, ld, full = stopped.m, stopped.e, stopped.ld, stopped.full
OUT = m.BASE/'production-harness-recovered'
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
    assert not OUT.exists()
    old = stopped.bindings()
    evidence = stopped.OUT/'interruption-20260915.json'
    interrupted = json.loads(evidence.read_text())
    assert interrupted['accepted_step'] == 53 and not interrupted['native_failure_file_present']
    for rel, digest in interrupted['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    source = stopped.OUT/'path-1'
    assert not (source/'failure.json').exists() and not (source/'result.json').exists()
    assert not (source/'step-0054.npz').exists()
    OUT.mkdir()
    prefix = OUT/'accepted-prefix-1'
    prefix.mkdir()
    times = full.time_nodes(old, 1)
    for step in range(54):
        path = source/f'step-{step:04d}.npz'
        with np.load(path) as saved:
            assert saved['time_seconds'] == times[step]
        os.link(path, prefix/path.name)
    rows = [json.loads(line) for line in (source/'iterations.jsonl').read_text().splitlines()]
    rows = [row for row in rows if row['step'] <= 53]
    assert {row['step'] for row in rows} == set(range(1, 54))
    for step in range(1, 54):
        records = [row for row in rows if row['step'] == step]
        assert [row['iteration'] for row in records] == list(range(len(records)))
        assert len(records) <= 24 and records[-1]['residual_norm'] <= 1
    (prefix/'iterations.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in rows))
    e.write(prefix/'manifest.json', dict(classification='Counterexample candidate', accepted_prefix_steps=53,
        interruption_sha256=e.digest(evidence), sha256={p.name:e.digest(p) for p in prefix.iterdir()}))
    plan = dict(old, phase='Original full duration after a WSL lifetime interruption',
        additional_steps_to_compute=245, imported_prefix_steps_by_refinement={'1':53, '2':64, '4':128},
        imported_prefix='Replay accepted 53/64/128 states with exact original clocks and both BDF state/flux histories. Recompute unfinished coarse step 54 from accepted step 53.',
        numerical_solver_unchanged_from_donor_repair=True, interruption=evidence.relative_to(e.ROOT).as_posix())
    plan['imported_paths'] = dict(old['imported_paths'])
    plan['imported_paths']['1'] = dict(steps=53, directory=prefix.relative_to(e.ROOT).as_posix(), refinement=1,
        source_plan=(stopped.OUT/'plan.json').relative_to(e.ROOT).as_posix())
    plan['bindings'] = dict(old['bindings'])
    for p in [Path(__file__), Path(__file__).with_name('run_gr_harness_recovery.ps1'),
              evidence, prefix/'manifest.json', stopped.OUT/'plan.json']:
        plan['bindings'][p.relative_to(e.ROOT).as_posix()] = e.digest(p)
    assert plan['exact_time_nodes_seconds'] == old['exact_time_nodes_seconds']
    assert local['cached'].stage is stopped.stage
    e.write(OUT/'plan.json', plan)
    check()
    print('HARNESS RECOVERY PREPARED; ORIGINAL GRID, 245 NEW STEPS', flush=True)


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
                command=['verification/gr_harness_recovery.py', 'chain', '--workers', str(args.workers)],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                start_ticks=(Path('/proc')/str(pid)/'stat').read_text().rsplit(') ',1)[1].split()[19],
                plan_sha256=e.digest(OUT/'plan.json'), cpu_affinity=sorted(os.sched_getaffinity(0))), stream, indent=2)
        for r in (1, 2, 4):
            run(r, args.workers)
        compare()
    else:
        compare()
