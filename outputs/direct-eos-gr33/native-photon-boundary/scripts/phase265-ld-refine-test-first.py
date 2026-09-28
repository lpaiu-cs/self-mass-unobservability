"""Scratch test of the production long-double refine (apply_photon_geometric_boundary_ld) on the rebuilt coarse step-119 system.

Exercises the exact fallback code, including its JSON log write, in a scratch directory. No production output.
"""
import inspect, json, os, shutil, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
import apply_photon_geometric_boundary_ld as ldmod
first = ldmod.first; joint = first.joint; LD = np.longdouble
SRC = ldmod.HGATE
scratch = Path('scratch-ld-refine265')
result = dict(started=time.strftime('%H:%M:%S'), passed=False, json_written=False)
def dump(): Path('.phase265-ld-refine-test.json').write_text(json.dumps(result, indent=1) + '\n')
if scratch.exists(): shutil.rmtree(scratch)
files = list((SRC/'sweep-0').rglob('*.npz')) + [p for part in ['gr', 'metric'] for p in (SRC/part).iterdir() if p.is_file()]
files += [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json', 'symbolic.json', 'endpoint-check.json', 'metric-result.json']]
for p in files:
    dst = scratch/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
for part in ['sweep-1/photons', 'sweep-1/material']: (scratch/part).mkdir(parents=True, exist_ok=True)
first.OUT = ldmod.OUT = ldmod.hgate.OUT = scratch; first.TAG[0] = 'test'
try:
    Model = ldmod.initialize(); m = Model(64)
    last = dict(np.load(SRC/'last-pair-64.npz')); flags = np.load(SRC/'sweep-1/photons/interval-15-64.npz')['split_macro_steps']
    base_h = m.t[-1]/64; clock = []
    for k in range(64):
        parts = 1 + int(flags[k]); step = base_h/parts; clock.extend((k*base_h + sub*step, step) for sub in range(parts))
    t, hh = clock[118]; assert t == m.anchor['actual_step_edges'][118] and abs(last['time'] + last['step'] - t) < 1e-18
    x = last['photons'][-1].copy(); g = last['gas'][-1].copy(); m.guide_g = g.copy(); g[~m.material.active(t)] = 0.
    lus = [joint.splu(joint.sparse.eye(m.n*m.q, format='csc') - a*hh*m.A) for a in [5/12, 1/4]]
    captured = {}
    class Captured(Exception): pass
    def intercept(m_, op, P, rhs, guess): captured.update(op=op, P=P, rhs=rhs, guess=guess); raise Captured()
    stage = Model.run.__globals__['stages']
    try: first.bind(stage, solve=intercept)(m, t, hh, x, g, lus)
    except Captured: pass
    refine = inspect.getclosurevars(stage.__globals__['solve']).nonlocals['refine']
    s0 = time.monotonic(); sol = refine(m, captured['op'], captured['P'], captured['rhs'], captured['guess'])
    log = json.loads((scratch/'long-double-test.json').read_text())
    event = log['events'][-1]
    result.update(passed=bool(event['accepted']), json_written=True, corrections=event['corrections'], seconds=time.monotonic() - s0, final=event['rows'][-1],
                  solution_dtype=str(sol.dtype), linear_entries=len(inspect.getclosurevars(stage.__globals__['solve']).nonlocals['linear']))
except BaseException as exc:
    result.update(error=repr(exc)[:800])
finally:
    result['finished'] = time.strftime('%H:%M:%S'); dump(); shutil.rmtree(scratch, ignore_errors=True)
