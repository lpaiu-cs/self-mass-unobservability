"""Phase267 diagnostic: does the final driven model's first t=0 macro step depend on the inputs judged unused?

Diagnostic only. Same construction and first coarse macro step as .phase267-audit.py, in its own scratch folder.
Variant 'base' changes nothing. Variant 'unused' doubles every read array of the inputs that the code reading judged
unused by the joint dynamics (Phase152 background template, Phase155 photon-lift runs, Phase157 self-GR metric,
Phase178 guides, sweep-0 histories). Variant 'metric' additionally doubles the Phase155 incident metric knots.
Substitutes are served through np.load / Path.read_bytes with hit counts; the step output npz is kept for comparison.
Usage: python3 .phase267-identity.py <variant> <result folder>
"""
import io, json, os, pathlib, shutil, sys, time
from pathlib import Path
sys.path.insert(0, 'verification')
import numpy as np
from types import FunctionType
variant, result = sys.argv[1], Path(sys.argv[2]); result.mkdir(parents=True, exist_ok=False)
root = os.path.abspath('.') + '/'
KEEP = {'t', 'radius_E', 'accepted_angular_times', 'RK_A', 'RK_b', 'RK_c', 'actual_step_edges', 'split_macro_steps', 'front_cells',
        'front_times', 'energy_offset_t', 'native_neutral_stage_times', 'native_neutral_stage_weights', 'accepted_angular_quadrature_weights',
        'time_integrator', 'ledger', 'escape', 'radial_ports', 'material_energy_units', 'material_neutral_units'}
UNUSED = ['retained-metric-return152-work/fields/source-128.npz',
          'native-incident-drive155-work/photons-lift/steps-64-reference-128.npz', 'native-incident-drive155-work/photons-lift/steps-128-reference-128.npz',
          'native-incident-self-gr157-work/metric/metric-128-g8.npz',
          'native-stage-collisions178-work/material-64.npz', 'native-stage-collisions178-work/sweep-1/photons/pilot-64.npz']
METRIC = ['native-incident-drive155-work/metric/corrected/metric-64-g8.npz', 'native-incident-drive155-work/metric/corrected/metric-128-g8.npz']
targets = {'base': [], 'unused': UNUSED, 'metric': UNUSED + METRIC}[variant]
def doubled(path):
    z = np.load(path); out = {}
    for k in z.files:
        v = z[k]; out[k] = v*2 if (k not in KEEP and v.dtype.kind == 'f') else v
    s = io.BytesIO(); np.savez(s, **out); return s.getvalue()
SUBST = {root + p: doubled(p) for p in targets}; hits = {p: 0 for p in SUBST}
original_read = pathlib.Path.read_bytes; original_load = np.load
def read_bytes(self):
    p = os.path.abspath(str(self))
    if p in SUBST: hits[p] += 1; return SUBST[p]
    return original_read(self)
pathlib.Path.read_bytes = read_bytes
def load(file, *a, **k):
    if isinstance(file, (str, Path)) and os.path.abspath(str(file)) in SUBST:
        p = os.path.abspath(str(file)); hits[p] += 1; file = io.BytesIO(SUBST[p])
    return original_load(file, *a, **k)
np.load = load
SRC = Path('native-true-momentum238-work'); scratch = Path(f'scratch-identity267-{variant}'); assert not scratch.exists()
for p in list((SRC/'sweep-0').rglob('*.npz')) + [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
    dst = scratch/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True)
    if variant != 'base' and '/sweep-0/' in str(p): dst.write_bytes(doubled(p)); hits[str(dst)] = 'rewritten'
    else: os.link(p, dst)
for part in ['sweep-1/photons', 'sweep-1/material']: (scratch/part).mkdir(parents=True, exist_ok=True)
import continue_true_momentum as ctm
ctm.OUT = scratch; row = dict(variant=variant)
try:
    t0 = time.monotonic(); FunctionType(ctm.initialize.__code__, dict(ctm.initialize.__globals__, OUT=scratch))(False)
    m = ctm.owner.Model(64); row['construct_seconds'] = time.monotonic() - t0
    t0 = time.monotonic(); out = m.run(64, f'identity-{variant}-64', 1); row['run_seconds'] = time.monotonic() - t0
    row['row'] = {k: v for k, v in out.items() if not isinstance(v, (list, dict))}
except BaseException as exc:
    row['error'] = repr(exc)[:2000]
finally:
    row['hits'] = {k.replace(root, ''): v for k, v in hits.items()}
    for f in (scratch/'sweep-1').rglob('*'):
        if f.is_file() and f.suffix in ('.npz', '.json'): shutil.copy2(f, result/f.name)
    (result/'identity.json').write_text(json.dumps(row, indent=1) + '\n')
    shutil.rmtree(scratch, ignore_errors=True)
