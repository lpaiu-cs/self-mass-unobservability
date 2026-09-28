"""Phase268 preparation steps and checks (run inside the 4x runtime).

inputs : copy the grid-independent inputs recorded by phase 267 (verified against the original runtime) and the
         four extra grid-independent files phase 267 placed by hand.
bank   : 43 interior cells; the 2x edges are contained bitwise; cells 0-7 equal the 1x arrays bitwise.
thermal: 43-cell thermal/rate bank; cells 0-7 and the thin cells equal the 1x arrays.
initial: 555-cell initial state and a native constitutive error below the original 0.2 percent.
ratio  : radius_E/re of the original grid (printed for the phase-150 zero metric).
exists : every existence check of the smoke step agrees with the original runtime.
"""
import hashlib, json, os, shutil, sys
from pathlib import Path
import numpy as np
OLD = Path('/home/lpaiu/work/native-retained-tail-runtime'); R2 = Path('/home/lpaiu/work/native-refined267-runtime'); NEW = Path('.').resolve()
mode = sys.argv[1]
if mode == 'inputs':
    log = json.loads((R2/'.phase267-copied-inputs.json').read_text())
    extra = ['outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/lapse/corrected/metric-128-g8.npz',
             'native-incident-drive155-work/material-precision-plan.json', 'native-pressure-return173-work/sweep-1/expanded-front-run.py',
             'native-driver-radau212-work/expanded-source.py']
    for p in list(log) + extra:
        h = hashlib.sha256((OLD/p).read_bytes()).hexdigest()
        if p in log: assert h == log[p], ('changed since phase 267', p)
        Path(p).parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(OLD/p, p)
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == h; log[p] = h
    Path('.phase268-copied-inputs.json').write_text(json.dumps(log, indent=1, sort_keys=True) + '\n'); print('copied grid-independent inputs', len(log))
elif mode == 'bank':
    rel = 'outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz'
    g, o, t = np.load(rel), np.load(OLD/rel), np.load(R2/rel)
    print('cells', len(g['r']), 'edges', len(g['edges'])); assert len(g['r']) == 43
    print('widths km', np.round(np.diff(g['edges'])/1e5, 3).tolist())
    nested = set(t['edges'].tolist()) <= set(g['edges'].tolist()); print('2x edges contained bitwise:', nested); assert nested
    keys = ['r', 'a', 'B', 'rho', 'T', 'phi', 'raw', 'y0', 'thermo', 'target', 'coeff']
    kept = all(np.array_equal(g[k][:8], o[k][:8]) for k in keys); print('cells 0-7 bitwise vs 1x:', kept); assert kept
    print('thin cells 40-42 vs 1x 16-18:', {k: bool(np.array_equal(g[k][40:], o[k][16:])) for k in ['r', 'rho', 'T', 'raw', 'coeff', 'thermo']})
elif mode == 'thermal':
    rel = 'outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/bank.npz'; a, b = np.load(rel), np.load(OLD/rel)
    print('thermal raw', a['raw'].shape, 'cells 0-7 bitwise', bool(np.array_equal(a['raw'][:8], b['raw'][:8]) and np.array_equal(a['rates'][:8], b['rates'][:8])),
          'thin 40-42 vs 16-18', bool(np.array_equal(a['raw'][40:], b['raw'][16:]) and np.array_equal(a['rates'][40:], b['rates'][16:])))
    assert a['raw'].shape[0] == 43
elif mode == 'initial':
    F = 'outputs/direct-eos-gr33/def-native-initial-constraints/finite-volume'
    z = np.load(F + '/balanced-initial-state.npz'); a = json.load(open(F + '/audit.json'))
    print('initial state', z['radius_E'].shape, 'anchors', len(a['native_anchors']), 'native err', a['native_constitutive_relative'],
          'inventory', a['actual_cell_inventory_max_relative'], 'matched', a['constraint_input_matched'])
    assert z['radius_E'].shape[0] == 555 and a['native_constitutive_relative'] < 0.002
elif mode == 'ratio':
    early = 'outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/lapse/corrected/metric-128-g8.npz'
    z, t = np.load(OLD/early), np.load(OLD/'outputs/direct-eos-gr33/def-native-anisotropic-gr/source-128.npz')
    R = z['radius_E'][-1]/t['re'][-1]; assert np.array_equal(t['re']*R, z['radius_E']); print(repr(float(R)))
elif mode == 'exists':
    d = json.load(open('.phase268-exists-smoke.json')); bad = []
    for p, r in d.items():
        if not p.startswith(str(NEW) + '/'): continue
        rel = p[len(str(NEW)) + 1:]
        if rel.startswith(('primary268-', 'verification/', '.')): continue
        if r != (OLD/rel).exists(): bad.append((r, rel))
    # def_retained_motion_return.py:57-58 creates this alias of the immutable cache when it is missing (prior.initialize
    # then verifies the cache SHA); the original runtime created it on its first run. Accepted only if the alias now
    # exists and points at this runtime's copy of the same relative target.
    alias = 'retained-motion-return151-work/immutable-coupled-128.npz'
    if (False, alias) in bad and os.readlink(alias) == str(NEW/os.path.relpath(os.readlink(OLD/alias), OLD)):
        bad.remove((False, alias)); print('self-created alias accepted:', alias, '->', os.readlink(alias))
    print('existence checks', len(d), 'disagreements', bad); assert not bad
else: raise ValueError(mode)
