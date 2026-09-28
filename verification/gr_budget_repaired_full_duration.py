"""Counterexample candidate: restart the full GR interval within 24 iterations.

Reuse the exact accepted histories, never the out-of-budget probe state.
Only the last-correction dispatch changes; all native scientific gates remain.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
from types import FunctionType, SimpleNamespace

import gr_iteration_budget_recovery as repair

stopped, m, e, full = repair.stopped, repair.m, repair.e, repair.full
OUT = repair.OUT
cached = SimpleNamespace(**dict(vars(stopped.local['cached']), stage=repair.stage))
context = dict(stopped.context, OUT=OUT, stage=repair.stage)
bindings = FunctionType(m.prior.bindings.__code__, context)
context['bindings'] = bindings
completed = FunctionType(m.prior.completed.__code__, context)
context['completed'] = completed
local = dict(stopped.local, OUT=OUT, cached=cached, bindings=bindings, context=context)
check = FunctionType(full.check.__code__, local)
run = FunctionType(full.run.__code__, local)
common_times = FunctionType(full.common_times.__code__, local)
compare = FunctionType(full.compare.__code__, dict(local, common_times=common_times))


def prepare():
    assert not OUT.exists()
    repair.budget_check()
    old = stopped.bindings()
    probe = json.loads((repair.PROBE/'result.json').read_text())
    assert probe['correction_passed'] and not probe['original_failed_step_accepted']
    for manifest in [repair.PROBE/'manifest.json', stopped.OUT/'failure-manifest.json']:
        for rel, digest in json.loads(manifest.read_text())['sha256'].items():
            assert e.digest(e.ROOT/rel) == digest, rel
    assert old['imported_prefix_steps_by_refinement'] == {'1':53, '2':64, '4':128}
    plan = dict(old, phase='Original full duration with a reserved final donor correction',
        candidate='Same native GR equations and exact 70/140/280 grid; consistent donor correction is available before the original iteration budget is exhausted.',
        nonlinear_operator='Use the validated donor fallback on a native line-search failure. Additionally, after logging unconverged row 22, reserve the last allowed correction for that same donor solver and log its result as row 23. Initial row 0 plus at most 23 corrections, original 24-row limit. No additional correction after row 23. All native residual, budget, cone and endpoint gates are unchanged.',
        numerical_solver_unchanged_from_donor_repair=False,
        native_equations_unchanged=True, original_time_grid_unchanged=True,
        failed_runtime_plan=(stopped.OUT/'plan.json').relative_to(e.ROOT).as_posix(),
        failed_native_manifest=(stopped.OUT/'failure-manifest.json').relative_to(e.ROOT).as_posix(),
        counterfactual_probe=(repair.PROBE/'manifest.json').relative_to(e.ROOT).as_posix(),
        counterfactual_state_imported=False,
        imported_prefix='Replay the accepted 53/64/128 histories with their exact original clocks and BDF state/flux histories. Recompute coarse step 54 from accepted step 53 with the original 24-row limit. The counterfactual probe state is never imported.')
    plan['bindings'] = dict(old['bindings'])
    for p in [Path(__file__), Path(repair.__file__), Path(__file__).with_name('run_gr_budget_recovery.ps1'),
              repair.PROBE/'manifest.json', repair.PROBE/'result.json',
              stopped.OUT/'failure-manifest.json', stopped.OUT/'plan.json']:
        plan['bindings'][p.relative_to(e.ROOT).as_posix()] = e.digest(p)
    assert plan['exact_time_nodes_seconds'] == old['exact_time_nodes_seconds']
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    check()
    print('BUDGET RECOVERY PREPARED; RECOMPUTE 54 FROM ACCEPTED 53', flush=True)


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
                command=['verification/gr_budget_repaired_full_duration.py', 'chain', '--workers', str(args.workers)],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                start_ticks=(Path('/proc')/str(pid)/'stat').read_text().rsplit(') ',1)[1].split()[19],
                plan_sha256=e.digest(OUT/'plan.json'), cpu_affinity=sorted(os.sched_getaffinity(0))), stream, indent=2)
        for refinement in (1, 2, 4):
            run(refinement, args.workers)
        compare()
    else:
        compare()
