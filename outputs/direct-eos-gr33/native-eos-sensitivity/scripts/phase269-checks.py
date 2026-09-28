"""Phase269 stage 2 preparation checks (2x grid, PL-off EOS and level libraries; run inside the new runtime).

modules: every verification module equals the original runtime's except the two library constants, and those differ
         only in the intended line.
inputs : copy the grid-independent inputs recorded by phase 267 (verified against the original runtime) and the four
         extra files phase 267 placed by hand (all PL/MHD background inputs, unchanged by design).
bank   : 27 interior cells; geometry bitwise equal to the phase-267 2x bank (kept cells included); the EOS-derived keys
         of every cell (all from the PL-off native EOS, bank mode refined-all) are reported against the PL/MHD 2x bank.
thermal: 27-cell thermal/rate bank; differences of every cell from the PL/MHD 2x bank reported.
initial: 539-cell initial state and a native constitutive error below the original 0.2 percent.
ratio  : radius_E/re of the original grid (printed for the phase-150 zero metric).
exists : every existence check of the smoke step agrees with the original runtime (self-created phase-151 alias as in 268).
"""
import hashlib, json, os, shutil, sys
from pathlib import Path
import numpy as np
OLD = Path('/home/lpaiu/work/native-retained-tail-runtime'); R2 = Path('/home/lpaiu/work/native-refined267-runtime'); NEW = Path('.').resolve()
PATCHED = {'def_native_cold_population.py': ("CACHE=old.CACHE.parent/'native-cold-population'", "CACHE=old.CACHE.parent/'native-cold-population-ploff'"),
           'def_native_hydrogen_exchange.py': ("LEVELS=Path('/home/lpaiu/work/direct-eos-gr33/photon-eos-levels-repaired/levels.so')",
                                               "LEVELS=Path('/home/lpaiu/work/direct-eos-gr33/photon-eos-levels-ploff/levels.so')")}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
mode = sys.argv[1]
if mode == 'patch':
    for name, (a, b) in PATCHED.items():
        p = Path('verification')/name; t = p.read_text(); assert t.count(a) == 1, name; p.write_text(t.replace(a, b))
    print('patched', list(PATCHED))
elif mode == 'modules':
    old = {p.name: sha(p) for p in (OLD/'verification').glob('*.py')}; new = {p.name: sha(p) for p in Path('verification').glob('*.py')}
    assert set(old) == set(new), sorted(set(old) ^ set(new))[:5]
    differing = sorted(n for n in old if old[n] != new[n]); assert differing == sorted(PATCHED), differing
    for name, (a, b) in PATCHED.items():
        o = (OLD/'verification'/name).read_text().splitlines(); n = (Path('verification')/name).read_text().splitlines()
        changed = [(x, y) for x, y in zip(o, n) if x != y]; assert len(o) == len(n) and changed == [(a, b)], (name, changed)
    print('modules', len(new), 'identical except the two library constants:', differing)
elif mode == 'inputs':
    log = json.loads((R2/'.phase267-copied-inputs.json').read_text())
    extra = ['outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/lapse/corrected/metric-128-g8.npz',
             'native-incident-drive155-work/material-precision-plan.json', 'native-pressure-return173-work/sweep-1/expanded-front-run.py',
             'native-driver-radau212-work/expanded-source.py']
    for p in list(log) + extra:
        h = sha(OLD/p)
        if p in log: assert h == log[p], ('changed since phase 267', p)
        Path(p).parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(OLD/p, p); assert sha(p) == h; log[p] = h
    Path('.phase269-copied-inputs.json').write_text(json.dumps(log, indent=1, sort_keys=True) + '\n'); print('copied grid-independent inputs', len(log))
elif mode == 'bank':
    rel = 'outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz'; g, t = np.load(rel), np.load(R2/rel)
    print('cells', len(g['r'])); assert len(g['r']) == 27
    geometry = ['r', 'a', 'B', 'rho', 'T', 'phi', 'edges', 'face_a', 'face_B', 'face_T', 'face_L', 'Einf', 'num', 'cx']
    same = {k: bool(np.array_equal(g[k], t[k])) for k in geometry}; print('geometry bitwise vs PL/MHD 2x:', same); assert all(same.values())
    for k in ['y0', 'thermo', 'raw', 'coeff']:  # every cell's EOS-derived arrays now come from the PL-off native EOS
        a, b = np.asarray(t[k], float), np.asarray(g[k], float)
        rel = np.abs(b - a)/np.maximum(np.abs(a), 1e-300); rel = rel.reshape(len(rel), -1)
        print(k, 'cells 0-26 max relative change per cell', np.round(np.nanmax(np.where(np.abs(a.reshape(len(a), -1)) > 1e-280, rel, 0), axis=1), 4).tolist())
    print('native derivative errors (PL-off)', np.round(g['derivative_errors'], 10).tolist() if 'derivative_errors' in g.files else None)
elif mode == 'thermal':
    rel = 'outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/bank.npz'; a, b = np.load(rel), np.load(R2/rel)
    assert a['raw'].shape == b['raw'].shape and a['raw'].shape[0] == 27
    r = np.abs(a['rates'] - b['rates'])/np.maximum(np.abs(b['rates']), 1e-300)
    print('rates cells 0-26 max relative change per cell vs PL/MHD 2x', np.round(r.reshape(27, -1).max(1), 4).tolist())
elif mode == 'initial':
    F = 'outputs/direct-eos-gr33/def-native-initial-constraints/finite-volume'
    z = np.load(F + '/balanced-initial-state.npz'); a = json.load(open(F + '/audit.json'))
    print('initial state', z['radius_E'].shape, 'anchors', len(a['native_anchors']), 'native err', a['native_constitutive_relative'],
          'inventory', a['actual_cell_inventory_max_relative'], 'matched', a['constraint_input_matched'])
    assert z['radius_E'].shape[0] == 539 and a['native_constitutive_relative'] < 0.002
elif mode == 'ratio':
    early = 'outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/lapse/corrected/metric-128-g8.npz'
    z, t = np.load(OLD/early), np.load(OLD/'outputs/direct-eos-gr33/def-native-anisotropic-gr/source-128.npz')
    R = z['radius_E'][-1]/t['re'][-1]; assert np.array_equal(t['re']*R, z['radius_E']); print(repr(float(R)))
elif mode == 'exists':
    d = json.load(open('.phase269-exists-smoke.json')); bad = []
    for p, r in d.items():
        if not p.startswith(str(NEW) + '/'): continue
        rel = p[len(str(NEW)) + 1:]
        if rel.startswith(('primary269-', 'verification/', '.')): continue
        if r != (OLD/rel).exists(): bad.append((r, rel))
    alias = 'retained-motion-return151-work/immutable-coupled-128.npz'
    if (False, alias) in bad and os.readlink(alias) == str(NEW/os.path.relpath(os.readlink(OLD/alias), OLD)):
        bad.remove((False, alias)); print('self-created alias accepted:', alias, '->', os.readlink(alias))
    print('existence checks', len(d), 'disagreements', bad); assert not bad
else: raise ValueError(mode)
