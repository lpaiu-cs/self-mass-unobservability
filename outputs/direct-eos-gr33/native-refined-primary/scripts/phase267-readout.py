"""Phase267 step H: endpoint compact charge of a 64-clock primary history (phase 244 source, phase 248/266 endpoint field).

Counterexample candidate. Uses the unchanged producers: read_full_captured_history.joined (capture assembly),
read_complete_radau_history.endpoints/source (dense GR source), extend_retarded_history Response.setup/at (retarded
field at T). Only the clock list is cut to [64]; the 64/128 source-time comparison is therefore not evaluated.
The 248 output clock is built from the 64 run alone; the endpoint value needs only t=T. Each stage is its own
process, as in production (the owners' initialize edits shared sources once per process).
Usage: python3 .phase267-readout.py endpoints <out> <complete-64 npz> <recovered-64 npz | captures:<folder>>
       python3 .phase267-readout.py source <out> [<reference 244 source npz>]
       python3 .phase267-readout.py field <out> [<evaluation time s; default the history's last time>]
"""
import inspect, json, os, sys, time
from pathlib import Path
sys.path.insert(0, 'verification')
import numpy as np
from scipy.interpolate import PPoly
import read_full_captured_history as cap
import extend_retarded_history as ext
rch = cap.prior
stage, out = sys.argv[1], Path(sys.argv[2]); report = out/(f'readout-{stage}.json' if not (stage == 'field' and len(sys.argv) > 3) else f'readout-field-{sys.argv[3]}.json'); result = dict(classification='Counterexample candidate', stage=stage)
def dump(): report.write_text(json.dumps(result, indent=1) + '\n')
# Stage states accepted by the phase-238 production polish can hold conserved coordinates finer than any long-double gas
# state (no exact preimage within the owner's 8-ulp search). Only then, the nearest preimage from the same search is used
# and its relative coordinate mismatch is recorded and required to stay below 1e-16; exact states are unchanged.
exact_inverse = rch.recovery.prior.restored_gas; inexact_inverses = []
def restored_gas(m, q):
    try: return exact_inverse(m, q)
    except AssertionError as exc:
        if not (exc.args and isinstance(exc.args[0], tuple) and exc.args[0][0] == 'Coordinate inverse has no nearby exact preimage'): raise
    LD = np.longdouble; g = np.column_stack([(q[2]-m.kappa*q[0])/m.eu, q[3]/m.nu, q[0]/m.bu, q[1]/m.su])
    for component, coordinate in [(2, 0), (3, 1), (1, 3), (0, 2)]:
        for attempt in range(8):
            actual = m.conserved(g)[coordinate]; mask = actual != q[coordinate]
            if not np.any(mask): break
            target = np.where(actual[mask] < q[coordinate, mask], LD('inf'), LD('-inf')); g[mask, component] = np.nextafter(g[mask, component], target)
    mismatch = float(np.max(abs(m.conserved(g) - q)/np.maximum(abs(q), LD('1e-290'))))
    inexact_inverses.append(mismatch); assert mismatch < 1e-16, ('Nearest coordinate preimage too far', mismatch)
    return g
rch.recovery.prior.restored_gas = restored_gas
def patched(fn, changes, **values):
    s = inspect.getsource(fn)
    for a, b in changes: assert s.count(a) == 1, (fn.__name__, a); s = s.replace(a, b)
    ns = dict(fn.__globals__, **values); exec(compile(s, fn.__code__.co_filename, 'exec'), ns); return ns[fn.__name__]
if stage == 'endpoints':
    saved_path, photon = Path(sys.argv[3]), sys.argv[4]; assert not out.exists()
    SRC = Path('native-true-momentum238-work')
    for p in list((SRC/'sweep-0').rglob('*.npz')) + [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
        dst = out/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
    for part in ['sweep-1/photons', 'sweep-1/material', 'gr']: (out/part).mkdir(parents=True, exist_ok=True)
    os.link(saved_path, out/'complete-64.npz'); z = dict(np.load(saved_path))
    (out/'plan.json').write_text(json.dumps(dict(original_return_horizon_seconds=float(z['actual_step_edges'][-1]), saved=str(saved_path), photon=photon), indent=1) + '\n')
    start = time.monotonic()
    if photon.startswith('captures:'):
        folder = Path(photon.split(':', 1)[1]); count = len(z['joint_stage_times'])
        captures = [dict(np.load(folder/f'captured-64-{i:03d}.npz')) for i in range(count)]
        c0 = captures[0]; empty = lambda v: np.zeros((0,) + np.shape(v), np.asarray(v).dtype)
        prefix = dict(times=empty(c0['time']), weights=empty(c0['weight']), photon_moments=empty(c0['photon_moments']),
                      collision_rates=empty(c0['collision_rates']), radial_ports=empty(c0['radial_ports']), angular=empty(z['accepted_angular_luminosity'][0]))
        d, row = cap.joined(z, prefix, captures); np.savez_compressed(out/'recovered-64.npz', **d); result['assembly'] = row
    else: os.link(photon, out/'recovered-64.npz')
    result['assembly_seconds'] = time.monotonic() - start; dump()
    endpoints = patched(rch.endpoints, [
        ("    initialize=FunctionType(",
         "    source=replace(source,'    for n in [64,128]:','    for n in [64]:')\n"
         "    source=replace(source,'    errors={k:aligned(*outputs,k) for k in keys}','    errors={}')\n"
         "    source=replace(source,'passed=max(errors.values())<.02','passed=None')\n"
         "    initialize=FunctionType("),
        ("    for n in [64,128]:\n        os.rename", "    for n in [64]:\n        os.rename")],
        OUT=out, saved=lambda n: out/'complete-64.npz', recovered=cap.bind(rch.recovered, INPUT=out))
    start = time.monotonic(); endpoints(); result['seconds'] = time.monotonic() - start
    result['check'] = {k: v for k, v in json.loads((out/'endpoint-64-check.json').read_text()).items() if k != 'local_material_ledger'}
elif stage == 'source':
    source = patched(rch.source, [
        ("    for a,b in changes:s=replace(s,a,b)", "    changes.append(('    for n in [64,128]:','    for n in [64]:'))\n    for a,b in changes:s=replace(s,a,b)"),
        ("    for n in [64,128]:\n        d=dict(np.load(OUT/'gr'/f'source-{n}.npz'))", "    for n in [64]:\n        d=dict(np.load(OUT/'gr'/f'source-{n}.npz'))")],
        OUT=out, INPUT=out, saved=lambda n: out/'complete-64.npz')
    start = time.monotonic(); source(); result['seconds'] = time.monotonic() - start
    result['check'] = {k: v for k, v in json.loads((out/'source-64-check.json').read_text()).items() if k not in ('dense_stage', 'polynomial', 'endpoint')}
    if len(sys.argv) > 3:
        raw, r = np.load(out/'gr/source-64.npz'), np.load(sys.argv[3])
        result['reference'] = sys.argv[3]
        result['reference_source_differing'] = [k for k in r.files if not (k in raw.files and raw[k].shape == r[k].shape and np.array_equal(raw[k], r[k]))]
elif stage == 'field':
    raw = dict(np.load(out/'gr/source-64.npz')); z = np.load(out/'complete-64.npz')
    candidates = [0.] + list(z['joint_stage_times']) + list(z['actual_step_edges']); clock = []
    for t in sorted(candidates):
        if not clock or t - clock[-1] > 1e-18: clock.append(t)
    clock = np.array(clock); knots, co = rch.coefficients(raw)
    d = dict(raw, t=clock, original_clock=np.array(64)); d.update({k: PPoly(np.asarray(v[::-1], float), knots)(clock) for k, v in co.items()})
    ids = np.array([np.argmin(abs(clock - t)) for t in raw['t']]); assert np.max(abs(clock[ids] - raw['t'])) < 1e-18
    for k in rch.KEYS: d[k][ids] = raw[k]
    np.savez_compressed(out/'gr/field-source-64.npz', **d)
    start = time.monotonic(); cap.bind(rch.endpoint.initialize, OUT=out)()
    m = ext.Response(); ext.bind(ext.previous.prior.Response.setup, INPUT=out)(m, d, 8)
    when = float(sys.argv[3]) if len(sys.argv) > 3 else float(d['t'][-1]); assert when <= d['t'][-1]  # retarded field: causal in the source
    free = ext.at(m, m.source, np.array([when])); M = float(d['M_cm'])
    result.update(seconds=time.monotonic() - start, M_cm=M, T=float(d['t'][-1]), evaluation_time=when, cells=len(d['radius']), output_times=len(clock),
                  endpoint_compact_charge=-float(free[0][-1, -1])/M)
else: raise ValueError(stage)
result.update(inexact_inverse_count=len(inexact_inverses), inexact_inverse_max=max(inexact_inverses, default=None))
dump(); print(json.dumps({k: v for k, v in result.items() if k not in ('assembly', 'check')}), flush=True)
