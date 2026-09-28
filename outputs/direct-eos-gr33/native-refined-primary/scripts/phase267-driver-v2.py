"""Phase267 step G (v2): the final driven primary from t=0, in restartable segments.

Counterexample candidate. continue_true_momentum.initialize(False) installs the phase-238 equations and arithmetic
(true-native S at 80 digits in the Newton RHS and defect) without the phase-238 resume seed. Model.run keeps every
original gate. Linear solves use the model's own joint solve; only if it fails its original gate is the phase-238
production solver applied (final unseeded factory on the consistent stable operator, as in phases 232-238), started
from the joint solve's result. Each stage starts Newton from the current gas (phase-232 guide update). Every
accepted stage is captured as in phases 235-238 (photon moments, radial ports, collision rates) for the readout.
v1 (phase267-driver-evolve-v1.py) replayed the production evolve, whose solver skips the first GMRES and failed the
original linear gate at t=0 (3.16e-14); that failure is preserved.
Usage: python3 .phase267-driver.py <work> <label> <limit|-> <restart|-> [<audit json|->] [<deadline seconds>]
"""
import json, os, resource, sys, time
from pathlib import Path
from types import FunctionType
work, label = Path(sys.argv[1]), sys.argv[2]
limit = None if sys.argv[3] == '-' else int(sys.argv[3]); restart = None if sys.argv[4] == '-' else sys.argv[4]
audit = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] != '-' else None; deadline = int(sys.argv[6]) if len(sys.argv) > 6 else 14400
opened = {}
if audit:
    def hook(event, args):
        if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
            p = os.path.abspath(os.fsdecode(args[0])); m = args[1] if len(args) > 1 else None
            e = opened.setdefault(p, dict(reads=0, writes=0)); e['writes' if isinstance(m, str) and any(c in m for c in 'wax+') else 'reads'] += 1
    sys.addaudithook(hook)
if os.environ.get('PHASE267_EXISTS'):  # existence checks switch arithmetic paths but raise no audit event
    import atexit, pathlib
    checked = {}
    def recorder(original, first_arg_is_self):
        def wrapped(*a, **k):
            r = original(*a, **k); checked.setdefault(os.path.abspath(os.fsdecode(a[0] if not first_arg_is_self else str(a[0]))), bool(r)); return r
        return wrapped
    for name in ['exists', 'is_file', 'is_dir']: setattr(pathlib.Path, name, recorder(getattr(pathlib.Path, name), True))
    for name in ['exists', 'isfile', 'isdir']: setattr(os.path, name, recorder(getattr(os.path, name), False))
    atexit.register(lambda: Path(os.environ['PHASE267_EXISTS']).write_text(json.dumps(checked, indent=1)))
sys.path.insert(0, 'verification')
SRC = Path('native-true-momentum238-work')
if not work.exists():
    for p in list((SRC/'sweep-0').rglob('*.npz')) + [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
        dst = work/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
    for part in ['sweep-1/photons', 'sweep-1/material', 'captures']: (work/part).mkdir(parents=True)
resource.setrlimit(resource.RLIMIT_AS, (12*1024**3, 12*1024**3))
import numpy as np
import continue_true_momentum as ctm
import continue_precise_momentum as cpm
ctm.joint.previous.original.inf.incident.native.deadline(deadline)
owner = ctm.owner; fallbacks = []; cache = {}; arithmetic = []; logs = []
import stabilize_full_interval_residual as stable; stable.OUT = work  # as the production evolve does before initialize
FunctionType(ctm.initialize.__code__, dict(ctm.initialize.__globals__, OUT=work))(False)
factory = FunctionType(ctm.factory.__code__, dict(ctm.factory.__globals__, OUT=work))
run = owner.Model.run; stage = run.__globals__['stages']; plain = stage.__globals__['solve']
def solve(m, op, P, rhs, guess):
    try: return plain(m, op, P, rhs, guess)
    except AssertionError as exc:
        tb, sol = exc.__traceback__, None
        while tb:
            if tb.tb_frame.f_code is plain.__code__: sol = tb.tb_frame.f_locals.get('sol')
            tb = tb.tb_next
        if 'solver' not in cache: cache['solver'] = factory(logs, m.n_steps)
        fallbacks.append(dict(stage_count=len(m.stage_t), joint_solve_error=repr(exc.args)[:300]))
        (work/f'{label}-fallbacks.json').write_text(json.dumps(fallbacks, indent=1) + '\n')
        correct = cpm.stable_operator(m, op, arithmetic, True)
        return cache['solver'](m, correct, P, rhs, np.asarray(guess if sol is None else sol, float))
actual = FunctionType(stage.__code__, dict(stage.__globals__, solve=solve))
def guided(m, t, h, x, g, lus):
    m.guide_g = g.copy(); return actual(m, t, h, x, g, lus)
owner.Model.run = FunctionType(run.__code__, dict(run.__globals__, stages=guided), argdefs=run.__defaults__)
boundary = owner.Model.boundary_ports; AMP = ctm.joint.AMP
def captured(m, t, x):  # phases 235-238 (continue_integer_native.initialize, seed=True)
    value = boundary(m, t, x); indices = [i for i in range(len(m.stage_t)-2, len(m.stage_t)) if abs(m.stage_t[i]-t) < 1e-18]
    assert len(indices) == 1, ('Capture must be an actual accepted stage', t); i = indices[0]
    moments = np.array([np.sum(x*m.Eweight, axis=(1, 2)), np.sum(x*m.Eweight*m.model.bulk.mu2[None, :, None], axis=(1, 2)), np.sum(x*m.Nweight, axis=(1, 2))])*AMP
    np.savez_compressed(work/'captures'/f'captured-{m.n_steps}-{i:03d}.npz', time=t, weight=m.stage_h[i], photon_moments=moments,
                        radial_ports=value*AMP, collision_rates=m.stage_collision[i])
    return value
owner.Model.boundary_ports = captured
start = time.monotonic(); result = dict(label=label, limit=limit, restart=restart)
try:
    m = owner.Model(64); m.n_steps = 64; result['construct_seconds'] = time.monotonic() - start
    row = m.run(64, label, limit, restart)
    result.update(passed=bool(row['passed']), fallbacks=len(fallbacks), row={k: v for k, v in row.items() if not isinstance(v, (list, dict))})
except BaseException as exc:
    result.update(passed=False, fallbacks=len(fallbacks), error=repr(exc)[:3000]); raise
finally:
    result['seconds'] = time.monotonic() - start; result['peak_RSS_bytes'] = 1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (work/f'{label}-driver.json').write_text(json.dumps(result, indent=1) + '\n')
    if audit:
        base = os.path.abspath('.') + '/'
        files = [dict(path=p.replace(base, ''), **e) for p, e in sorted(opened.items())
                 if not ('/verification/' in p or p.endswith('.pyc') or '/usr/lib/' in p or 'site-packages' in p or p.startswith('/proc') or p.startswith('/dev'))]
        Path(audit).write_text(json.dumps(dict(files=files), indent=1) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'row'}), flush=True)
sys.exit(0 if result.get('passed') else 3)
