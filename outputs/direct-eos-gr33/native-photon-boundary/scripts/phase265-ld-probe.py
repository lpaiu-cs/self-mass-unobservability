"""Scratch probe: rebuild the failed coarse step-119 system of the hgate run; test long-double FGMRES on it.

Read-only for the failed run. A scratch directory with hard-linked inputs hosts the model.
Checks: (1) the first recorded short GMRES call is reproduced bitwise; (2) long-double FGMRES
corrections reach the unchanged linear gates; (3) the stage with that solver passes its
unchanged nonlinear gates. Nothing here is a production result.
"""
import inspect, json, os, shutil, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
import apply_photon_geometric_boundary_hgate as h
import solve_long_double_fgmres as ld
first = h.first; joint = first.joint; material = h.material; LD = np.longdouble
SRC = h.OUT
scratch = Path('scratch-ld-probe265')
log = dict(started=time.strftime('%H:%M:%S'))
def dump():
    Path('.phase265-ld-probe.json').write_text(json.dumps(log, indent=1, default=float) + '\n')
if scratch.exists(): shutil.rmtree(scratch)
files = list((SRC/'sweep-0').rglob('*.npz')) + [p for p in (SRC/'gr').iterdir() if p.is_file()] + [p for p in (SRC/'metric').iterdir() if p.is_file()]
files += [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json', 'symbolic.json', 'endpoint-check.json', 'metric-result.json']]
for p in files:
    dst = scratch/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
for part in ['sweep-1/photons', 'sweep-1/material']: (scratch/part).mkdir(parents=True, exist_ok=True)
first.OUT = h.OUT = scratch; first.TAG[0] = 'probe'
t0 = time.monotonic(); Model = first.initialize(); m = Model(64); log['init_seconds'] = time.monotonic() - t0
last = dict(np.load(SRC/'last-pair-64.npz'))
flags = np.load(SRC/'sweep-1/photons/interval-15-64.npz')['split_macro_steps']
base_h = m.t[-1]/64; clock = []
for k in range(64):
    parts = 1 + int(flags[k]); step = base_h/parts
    clock.extend((k*base_h + sub*step, step) for sub in range(parts))
t, hh = clock[118]
log.update(clock_steps=len(clock), t=float(t), h=float(hh), anchor_edge=float(m.anchor['actual_step_edges'][118]), last_end=float(last['time'] + last['step']))
assert len(clock) == 119 and t == m.anchor['actual_step_edges'][118] and abs(last['time'] + last['step'] - t) < 1e-18
x = last['photons'][-1].copy(); g = last['gas'][-1].copy(); m.guide_g = g.copy(); g[~m.material.active(t)] = 0.
lus = [joint.splu(joint.sparse.eye(m.n*m.q, format='csc') - a*hh*m.A) for a in [5/12, 1/4]]
captured = {}
class Captured(Exception): pass
def intercept(m_, op, P, rhs, guess): captured.update(op=op, P=P, rhs=rhs, guess=guess); raise Captured()
stage = Model.run.__globals__['stages']
try: first.bind(stage, solve=intercept)(m, t, hh, x, g, lus)
except Captured: pass
op, P, rhs, guess = captured['op'], captured['P'], captured['rhs'], captured['guess']
wrapper = stage.__globals__['solve']; cv = inspect.getclosurevars(wrapper).nonlocals
full = cv['full']; inner = full.__globals__['gmres']; calls = inspect.getclosurevars(inner).nonlocals['log']
recorded = json.loads((SRC/'right-calls-coarse.json').read_text())['calls'][312]
it = []
inner(op, np.asarray(rhs, float), x0=np.asarray(guess, float), M=P, rtol=1e-14, atol=0., restart=80, maxiter=1, callback=it.append, callback_type='pr_norm')
mine = calls[-1]
log.update(reproduced_history=mine['history'] == recorded['history'], reproduced_true=mine['true_relative'] == recorded['true_relative'],
           recorded_true=recorded['true_relative'], probe_true=mine['true_relative'], dimension=len(rhs))
dump()
def gates(sol):
    residual = rhs - op.matvec(sol); relative = float(ld.norm(residual)/ld.norm(rhs))
    moments = joint.physical_norm(m, residual)/joint.scales(m, rhs, sol); gas = material.gas_relative(m, residual, sol)
    return residual, relative, [float(v) for v in moments], gas
s0 = time.monotonic(); op.matvec(np.asarray(guess, LD)); log['LD_matvec_seconds'] = time.monotonic() - s0
s0 = time.monotonic(); P.matvec(np.asarray(guess, float)); log['double_precondition_seconds'] = time.monotonic() - s0
sol = np.asarray(guess, LD); rows = []
for k in range(8):
    residual, relative, moments, gas = gates(sol)
    ok = relative < 1e-14 and max(moments) < 1e-13 and h.gas_gate(gas)
    rows.append(dict(correction=k, relative=relative, moments=moments, gas=gas, passed=ok)); log['ld_refinement'] = rows; dump()
    if ok: break
    inner_log = []; s0 = time.monotonic()
    sol = sol + ld.fgmres(op.matvec, P.matvec, residual, restart=80, cycles=2, rtol=LD('1e-10'), log=inner_log)
    rows[-1].update(fgmres_iterations=len(inner_log), fgmres_seconds=time.monotonic() - s0, fgmres_last=inner_log[-1] if inner_log else None,
                    fgmres_trace=[r['relative'] for r in inner_log][::max(1, len(inner_log)//12)])
    dump()
log['linear_gates_passed'] = rows[-1]['passed']; dump()
# Nonlinear stage acceptance with long-double solves for every Newton proposal.
def ld_solve(m_, op_, P_, rhs_, guess_):
    s = np.asarray(guess_, LD)
    for k in range(8):
        r = rhs_ - op_.matvec(s); rel = float(ld.norm(r)/ld.norm(rhs_))
        mo = joint.physical_norm(m_, r)/joint.scales(m_, rhs_, s); ga = material.gas_relative(m_, r, s)
        log.setdefault('stage_linear', []).append(dict(correction=k, relative=rel, moments=[float(v) for v in mo], gas=ga)); dump()
        if rel < 1e-14 and max(mo) < 1e-13 and h.gas_gate(ga): m_.max_residual = max(m_.max_residual, rel); return s
        s = s + ld.fgmres(op_.matvec, P_.matvec, r, restart=80, cycles=2, rtol=LD('1e-10'))
    raise AssertionError(('long-double linear refinement', rel))
if log['linear_gates_passed']:
    s0 = time.monotonic()
    try:
        first.bind(stage, solve=ld_solve)(m, t, hh, x, g, lus)
        log.update(stage_passed=True, stage_audit=m.newton_iterations[-1])
    except AssertionError as exc:
        log.update(stage_passed=False, stage_error=repr(exc)[:600])
    log['stage_seconds'] = time.monotonic() - s0
log['finished'] = time.strftime('%H:%M:%S'); dump()
shutil.rmtree(scratch)
