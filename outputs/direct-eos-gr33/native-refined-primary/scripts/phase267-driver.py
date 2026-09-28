"""Phase267 step G (v2): the final driven primary from t=0, in restartable segments.

Counterexample candidate. continue_true_momentum.initialize(False) installs the phase-238 equations and arithmetic
(true-native S at 80 digits in the Newton RHS and defect) without the phase-238 resume seed. Model.run keeps every
original gate. Linear solves use the model's own joint solve; only if it fails its original gate is the phase-238
production solver applied (final unseeded factory on the consistent stable operator, as in phases 232-238), started
from the joint solve's result. Each stage starts Newton from the current gas (phase-232 guide update). Every
accepted stage is captured as in phases 235-238 (photon moments, radial ports, collision rates) for the readout.
v3: when the joint solve fails its linear gate, the phase-265 long-double FGMRES refinement (12 corrections,
restart 80, 2 cycles, inner 1e-10; accepted only at vector 1e-14, physical moments 1e-13 and every material component
1e-13) runs first; the phase-238 production solver is the last resort. v2 (phase267-driver-v2.py) tried the production
solver directly; at macro 61 it stalled near 5e-11 for more than 20 minutes (attempt preserved). v3b starts the
production solver from the best long-double solution (v3 restarted it from the joint solve's rough result).
v3c evaluates the long-double refinement on the production solver's consistent operator: with the double operator the
vector residual floor at macro 61 was 9.6e-12 (attempt 3 preserved).
v4 (user approval 2026-09-27, last segment only, PHASE267_VECTOR_EXCEPTION=1e-10): after the joint solve, the long-double
refinement on the model operator and a 5-minute production polish all miss the original gate, the best long-double
solution with physical moments and every material component below 1e-13 and vector residual <= 1e-10 is accepted and
recorded; nonlinear and all other gates are unchanged, and m.max_residual keeps only original-gate solves.
v5 (user approval 2, same segment, PHASE267_NONLINEAR_EXCEPTION=1e-10): after at least three Newton iterations miss the
1e-12 vector defect, an iterate with vector defect <= 1e-10, physical moments < 1e-13 and every material component of the
defect < 1e-13 is accepted and recorded (the installed stage source is patched after a bytecode identity check).
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
import solve_long_double_fgmres as ldm
import finish_returned_material_accuracy as matacc
LD = np.longdouble; physical_norm, scales = plain.__globals__['physical_norm'], plain.__globals__['scales']
LINEAR = plain.__globals__['radau'].prior.owner.reuse.LINEAR
def linear_failure(exc): return bool(exc.args) and isinstance(exc.args[0], tuple) and exc.args[0][0] == 'Four-moment linear residual'
def record(): (work/f'{label}-fallbacks.json').write_text(json.dumps(fallbacks, indent=1))
best = []  # smallest vector residual reached by the long-double refinement
candidate = []  # best long-double solution meeting the physical gates within the user-approved vector bound
EXCEPTION = float(os.environ.get('PHASE267_VECTOR_EXCEPTION') or 0)  # 2026-09-27 user approval, last segment only
import signal
class ProductionTimeout(Exception): pass
def capped(seconds, fn, *args):  # time cap for the production polish; the segment deadline alarm is restored
    remaining = signal.alarm(0); previous = signal.getsignal(signal.SIGALRM); began = time.monotonic()
    def expire(signum, frame): raise ProductionTimeout(seconds)
    signal.signal(signal.SIGALRM, expire); signal.alarm(seconds)
    try: return fn(*args)
    finally:
        signal.alarm(0); signal.signal(signal.SIGALRM, previous)
        if remaining: signal.alarm(max(1, int(remaining - (time.monotonic() - began))))
exceptions = []
def refine(m, op, P, rhs, start):  # phase 265 (apply_photon_geometric_boundary_ld.refine) with the primary's 1e-13 material gate
    sol = np.asarray(start, LD).copy(); rows = []; began = time.monotonic(); best.clear(); candidate.clear()
    A = op  # v4: the model operator (linear in long double); the consistent operator's linearity defect was 1.2e-11 at macro 61
    if os.environ.get('PHASE267_OPERATOR_PROBE') and len(m.stage_t) >= int(os.environ['PHASE267_OPERATOR_PROBE']):
        x = np.asarray(start, LD); split = np.random.default_rng(7).random(len(x)).astype(LD); nb = ldm.norm(np.asarray(rhs, LD))
        probe = {}
        for name, O in [('plain', op), ('consistent', A)]:
            y = np.asarray(O.matvec(x)); d = y - np.asarray(O.matvec(x*split)) - np.asarray(O.matvec(x - x*split))
            parts = [m.unpack(v) for v in d.reshape(2, -1)]
            probe[name] = dict(dtype=str(y.dtype), linearity_over_rhs=float(ldm.norm(np.asarray(d, LD))/nb),
                photon_part=float(np.sqrt(sum(np.sum(np.asarray(q[0], LD)**2) for q in parts))/nb), gas_part=float(np.sqrt(sum(np.sum(np.asarray(q[1], LD)**2) for q in parts))/nb),
                float_input_over_rhs=float(ldm.norm(y - np.asarray(O.matvec(np.asarray(x, float)), LD))/nb), size_over_rhs=float(ldm.norm(np.asarray(y, LD))/nb))
        fallbacks[-1]['operator_probe'] = probe; record()
    for k in range(13):
        residual = rhs - np.asarray(A.matvec(sol), LD); relative = float(ldm.norm(residual)/max(ldm.norm(rhs), LD('1e-290')))
        moments = physical_norm(m, residual)/scales(m, rhs, sol); gas = matacc.gas_relative(m, residual, sol)
        passed = bool(relative < 1e-14 and max(moments) < 1e-13 and max(gas) < 1e-13)
        rows.append(dict(correction=k, relative=relative, moments=[float(v) for v in moments], material_relative=[float(v) for v in gas], passed=passed))
        if not best or relative < best[0]: best[:] = [relative, sol.copy()]
        if EXCEPTION and relative <= EXCEPTION and max(moments) < 1e-13 and max(gas) < 1e-13 and (not candidate or relative < candidate[0]):
            candidate[:] = [relative, sol.copy(), float(max(moments)), float(max(gas))]
        if passed or k == 12: break
        inner = []; sol = sol + ldm.fgmres(A.matvec, P.matvec, residual, restart=80, cycles=2, rtol=LD('1e-10'), log=inner)
        rows[-1].update(fgmres_iterations=len(inner), fgmres_seconds=inner[-1]['seconds'] if inner else 0.)
    fallbacks[-1]['long_double'] = dict(accepted=passed, seconds=time.monotonic() - began, rows=rows); record()
    if not passed: raise AssertionError(('Four-moment linear residual', relative, [float(v) for v in moments]))
    LINEAR.append(dict(initial_info=1, corrections=len(rows)-1, extended_residual=relative, long_double_fgmres=True))
    m.max_residual = max(m.max_residual, relative); return sol
def solve(m, op, P, rhs, guess):
    try: return plain(m, op, P, rhs, guess)
    except AssertionError as exc:
        if not linear_failure(exc): raise
        tb, sol = exc.__traceback__, None
        while tb:
            if tb.tb_frame.f_code is plain.__code__: sol = tb.tb_frame.f_locals.get('sol')
            tb = tb.tb_next
        start = guess if sol is None else sol
        fallbacks.append(dict(stage_count=len(m.stage_t), joint_solve_error=repr(exc.args)[:300])); record()
        try: return refine(m, op, P, rhs, start)
        except AssertionError as exc2:
            if not linear_failure(exc2): raise
        if 'solver' not in cache: cache['solver'] = factory(logs, m.n_steps)
        fallbacks[-1]['production_solver'] = True; record()
        correct = cpm.stable_operator(m, op, arithmetic, True)
        fallbacks[-1]['production_start_relative'] = float(best[0]) if best else None; record()
        if not EXCEPTION: return cache['solver'](m, correct, P, rhs, np.asarray(best[1] if best else start, float))
        kept = list(candidate)
        try: return capped(int(os.environ.get('PHASE267_POLISH_CAP') or 300), cache['solver'], m, correct, P, rhs, np.asarray(best[1] if best else start, float))
        except (AssertionError, ProductionTimeout) as exc3:
            if isinstance(exc3, AssertionError) and not linear_failure(exc3): raise
            fallbacks[-1]['production_result'] = repr(exc3)[:200]; record()
        if not kept: raise AssertionError(('Four-moment linear residual', 'no long-double solution within the approved vector bound'))
        relative, solution, moment, gas = kept
        exceptions.append(dict(stage_count=len(m.stage_t), vector_relative=relative, physical_moment_max=moment, material_component_max=gas))
        fallbacks[-1]['vector_gate_exception'] = exceptions[-1]; record()
        LINEAR.append(dict(initial_info=1, corrections=-1, extended_residual=relative, long_double_fgmres=True, vector_gate_exception=True))
        return solution
NL_EXCEPTION = float(os.environ.get('PHASE267_NONLINEAR_EXCEPTION') or 0)  # 2026-09-27 user approval 2, last segment only
nonlinear_exceptions = []
if NL_EXCEPTION:
    source = (work/'expanded-high-native-stage.py').read_text()  # the installed stage (bytecode match; true-momentum precise defect)
    check = dict(stage.__globals__); exec(compile(source, stage.__code__.co_filename, 'exec'), check)
    assert check['stages'].__code__.co_code == stage.__code__.co_code, 'expanded stage source differs from the installed stage'
    gate = "        if relative<1e-12 and max(moments)<1e-13:break\n"
    assert source.count(gate) == 1
    source = source.replace(gate, gate + "        if NL_EXCEPTION and newton>=2 and relative<=NL_EXCEPTION and max(moments)<1e-13 and max(gas_relative(m,defect,sol))<1e-13:nl_accept(m,t,h,newton,relative,moments,defect,sol);break\n")
    def nl_accept(m, t, h, newton, relative, moments, defect, sol):
        nonlinear_exceptions.append(dict(time=float(t), step=float(h), newton=int(newton), vector_defect=float(relative), physical_moment_max=float(max(moments)),
                                         material_component_max=float(max(matacc.gas_relative(m, defect, sol)))))
        (work/f'{label}-nonlinear-exceptions.json').write_text(json.dumps(nonlinear_exceptions, indent=1))
    namespace = dict(stage.__globals__, NL_EXCEPTION=NL_EXCEPTION, gas_relative=matacc.gas_relative, nl_accept=nl_accept)
    exec(compile(source, stage.__code__.co_filename, 'exec'), namespace); stage = namespace['stages']
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
    result.update(nonlinear_exception=NL_EXCEPTION or None, nonlinear_exceptions=nonlinear_exceptions)
    result.update(vector_gate_exception=EXCEPTION or None, vector_gate_exceptions=exceptions, vector_gate_exception_max=max((e['vector_relative'] for e in exceptions), default=None))
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
