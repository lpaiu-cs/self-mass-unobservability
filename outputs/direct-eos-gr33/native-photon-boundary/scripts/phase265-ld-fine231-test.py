"""Scratch replay of the first long-double attempt's unrecorded fallback: fine step 231, full production stage.

Counterexample candidate (diagnostic only; no production output). The first long-double attempt
(native-photon-boundary-ld265-work) accepted fine steps 216..230; at step 231 both double solvers
failed (last double iterate: vector 1.48e-13, moments up to 4.6e-9) and the long-double refine ran.
Its log write then raised on a numpy bool, so the outcome of that fallback was never recorded. The
numpy bool can only come from a correction row whose vector residual was already below 1e-14 while a
physical moment was not, so at least one long-double correction had been applied.

This replay rebuilds that stage from the attempt's saved last accepted pair (step 230) in a scratch
directory and runs the production stage of the fixed module (short proposals, original double solver,
long-double refine). It records whether the solver calls reproduce the recorded ones bitwise (and which
recorded solve they match), every long-double event, and whether the stage passes its unchanged gates.
The second attempt computes exactly this stage at its fine step 231.
"""
import inspect, json, os, shutil, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
import apply_photon_geometric_boundary_ld as ldmod
first = ldmod.first; joint = first.joint
SRC = ldmod.LD1
scratch = Path('scratch-ld-fine231-265')
result = dict(classification='Counterexample candidate', started=time.strftime('%H:%M:%S'), stage_passed=None)
def dump(): Path('.phase265-ld-fine231-test.json').write_text(json.dumps(result, indent=1, default=float) + '\n')
if scratch.exists(): shutil.rmtree(scratch)
files = list((SRC/'sweep-0').rglob('*.npz')) + [p for part in ['gr', 'metric'] for p in (SRC/part).iterdir() if p.is_file()]
files += [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json', 'symbolic.json', 'endpoint-check.json', 'metric-result.json']]
for p in files:
    dst = scratch/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
for part in ['sweep-1/photons', 'sweep-1/material']: (scratch/part).mkdir(parents=True, exist_ok=True)
first.OUT = ldmod.OUT = ldmod.hgate.OUT = scratch; first.TAG[0] = 'fine231test'
recorded = json.loads((SRC/'right-calls-fine.json').read_text())
calls = []
try:
    s0 = time.monotonic(); Model = ldmod.initialize(); m = Model(128); result['init_seconds'] = time.monotonic() - s0
    last = dict(np.load(SRC/'last-pair-128.npz')); flags = np.load(SRC/'sweep-1/photons/interval-15-128.npz')['split_macro_steps']
    base_h = m.t[-1]/128; clock = []
    for k in range(128):
        parts = 1 + int(flags[k]); step = base_h/parts; clock.extend((k*base_h + sub*step, step) for sub in range(parts))
    t, hh = clock[230]
    result.update(clock_steps=len(clock), t=float(t), h=float(hh), last_end=float(last['time'] + last['step']))
    assert len(clock) == 231 and t == m.anchor['actual_step_edges'][230] and abs(last['time'] + last['step'] - t) < 1e-18
    x = last['photons'][-1].copy(); g = last['gas'][-1].copy(); m.guide_g = g.copy(); g[~m.material.active(t)] = 0.
    lus = [joint.splu(joint.sparse.eye(m.n*m.q, format='csc') - a*hh*m.A) for a in [5/12, 1/4]]
    stage = Model.run.__globals__['stages']
    full = inspect.getclosurevars(stage.__globals__['solve']).nonlocals['full']
    calls = inspect.getclosurevars(full.__globals__['gmres']).nonlocals['log']
    dump(); s0 = time.monotonic()
    try:
        stage(m, t, hh, x, g, lus)
        result.update(stage_passed=True, stage_audit=m.newton_iterations[-1] if m.newton_iterations else None)
    except AssertionError as exc:
        result.update(stage_passed=False, stage_error=repr(exc)[:1500])
    result['stage_seconds'] = time.monotonic() - s0
except BaseException as exc:
    result.update(error=repr(exc)[:1500])
finally:
    rc = recorded['calls']
    same = lambda a, b: a['history'] == b['history'] and a['true_relative'] == b['true_relative']
    match = [j for j in range(len(rc)) if calls and same(calls[0], rc[j])]
    result.update(scratch_calls=len(calls), first_call_matches_recorded=match, recorded_solves=recorded['solves'][-3:])
    if match:
        j = match[0]; n = min(len(calls), len(rc) - j)
        result.update(bitwise_calls_from_match=sum(1 for i in range(n) if same(calls[i], rc[j + i])), compared_calls=n)
    result.update(long_double_events=ldmod.EVENTS, finished=time.strftime('%H:%M:%S')); dump()
    shutil.rmtree(scratch, ignore_errors=True)
