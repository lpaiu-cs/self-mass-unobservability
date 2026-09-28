"""Counterexample candidate: reserve the final native iteration for donor repair.

The failed run stays failed. A counterfactual correction of its saved iterate
is a mechanism check only; production must start again from accepted step 53
and meet the original 24-iteration limit and all native scientific gates.
"""
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from types import FunctionType, SimpleNamespace

import numpy as np
import gr_harness_recovery as stopped

m, e, ld, full = stopped.m, stopped.e, stopped.ld, stopped.full
repair = stopped.stopped.repair
PROBE = m.BASE/'step54-budget-check'
OUT = m.BASE/'production-budget-repaired'


class ReserveDonor(RuntimeError):
    pass


def stage(star, previous, older, h, coefficients, log):
    def record(row):
        log(row)
        # Rows 0..22 have used 22 corrections. Reserve correction 23 and
        # its final row for the consistent donor solver, never a 25th row.
        if row['iteration'] == 22 and row['residual_norm'] > 1:
            raise ReserveDonor()
    try:
        return stopped.stopped.stage(star, previous, older, h, coefficients, record)
    except ReserveDonor:
        return repair.finish(star, star.last_delta, star.last_state,
                             previous, older, h, coefficients, log, 23)


def budget_check():
    calls, rows = [], []
    def original(star, previous, older, h, coefficients, log):
        for i in range(24):
            log(dict(iteration=i, residual_norm=2 if i < 23 else 1))
        raise AssertionError('The final correction was not reserved')
    def finish(star, delta, state, previous, older, h, coefficients, log, start):
        calls.append(start)
        log(dict(iteration=start, residual_norm=0.5))
        return delta, state
    test = FunctionType(stage.__code__, dict(globals(),
        stopped=SimpleNamespace(stopped=SimpleNamespace(stage=original)),
        repair=SimpleNamespace(finish=finish)))
    test(SimpleNamespace(last_delta=None, last_state=None), None, None, None, None, rows.append)
    assert calls == [23] and [r['iteration'] for r in rows] == list(range(24))


def probe():
    budget_check()
    assert not PROBE.exists()
    plan = stopped.bindings()
    failure_manifest = stopped.OUT/'failure-manifest.json'
    for rel, digest in json.loads(failure_manifest.read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    PROBE.mkdir()
    e.write(PROBE/'plan.json', dict(classification='Counterexample candidate',
        bindings={str(Path(__file__).relative_to(e.ROOT)):e.digest(Path(__file__)),
                  str(failure_manifest.relative_to(e.ROOT)):e.digest(failure_manifest)},
        rule='Test the already validated donor correction on the exhausted native iterate. This cannot accept the failed step or add iterations to its old budget. A new execution must replay accepted step 53 and reserve its last permitted correction.'))
    with ProcessPoolExecutor(max_workers=15, initializer=e.worker_init) as pool:
        star = m.initialize(pool)
        times = full.time_nodes(plan, 1)
        older, previous = [full.restored(star, stopped.OUT/'path-1', k, times[k]) for k in (52, 53)]
        h = times[54]-times[53]
        coefficients = m.weights(h, times[53]-times[52])
        star.composition_context = previous, older, h, coefficients
        star.composition_data = m.composition_response(star, *previous)
        with np.load(stopped.OUT/'path-1/failed-iterate.npz') as saved:
            assert saved['time_seconds'] == times[54]
            delta = saved['delta'].copy()
            y = star.base+delta
            star.material_cache.clear()
            star.material_cache.update({e.material_key(row):aux.copy() for row,aux in
                zip(zip(y[:,0],y[:,1],y[:,5:]),saved['aux'])})
        value, state = m.residual(star, delta, previous, older, h, coefficients)
        before = float(np.max(abs(value)/m.ATOL))
        assert before == 5.241030758002932, before
        cell, field = np.unravel_index(np.argmax(abs(value)/m.ATOL), value.shape)
        print('STEP54 FAILURE REPRODUCED', before, int(cell), int(field), flush=True)
        try:
            delta, value, state, trials = repair.correction(star, delta, state, previous, older, h, coefficients)
        except Exception as error:
            e.write(PROBE/'failure.json', dict(classification='Counterexample candidate', reason=repr(error)))
            raise
        after = float(np.max(abs(value)/m.ATOL))
        result = dict(classification='Counterexample candidate', before=before, after=after,
            cell=int(cell), field=int(field), donor_pivots=trials, correction_passed=after<=1,
            original_failed_step_accepted=False, actual_step_within_budget_verified=False,
            scope='Counterfactual mechanism check after an exhausted budget; do not import this state as accepted.')
        np.savez_compressed(PROBE/'counterfactual-state.npz', delta=delta, aux=state['aux'], residual=value)
        e.write(PROBE/'result.json', result)
        print('STEP54 DONOR CHECK', json.dumps(result), flush=True)
    e.write(PROBE/'manifest.json', dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p)
        for p in [Path(__file__), *PROBE.iterdir()] if p.is_file()}))


if __name__ == '__main__':
    probe()
