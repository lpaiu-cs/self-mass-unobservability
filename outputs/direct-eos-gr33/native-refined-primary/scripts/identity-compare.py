"""Compare the first-step outputs of the identity variants with the base variant (bitwise and relative)."""
import json, sys
from pathlib import Path
import numpy as np
R = Path(sys.argv[1]); out = {}
for v in ['unused', 'metric']:
    info = json.loads((R/v/'identity.json').read_text()); row = dict(error=info.get('error'), hits=info.get('hits'))
    for suffix in ['-64.npz', '-64-checkpoint.npz']:
        a, b = R/'base'/f'identity-base{suffix}', R/v/f'identity-{v}{suffix}'
        if not (a.exists() and b.exists()): row[suffix] = 'missing'; continue
        za, zb = np.load(a, allow_pickle=True), np.load(b, allow_pickle=True); diffs = {}
        for k in za.files:
            x, y = za[k], zb[k]
            if x.shape != y.shape: diffs[k] = f'shape {x.shape} {y.shape}'; continue
            if np.array_equal(x, y, equal_nan=x.dtype.kind in 'fc'): continue
            if x.dtype.kind in 'fc': diffs[k] = float(np.max(abs(x - y)) / max(np.max(abs(x)), 1e-300))
            else: diffs[k] = 'differs'
        row[suffix] = dict(arrays=len(za.files), differing=diffs)
    ra = json.loads((R/'base'/'identity.json').read_text()).get('row', {}); rb = info.get('row', {})
    row['row_differences'] = {k: (ra.get(k), rb.get(k)) for k in ra if k not in ('seconds', 'stepping_seconds', 'operator_point_seconds') and ra.get(k) != rb.get(k)}
    out[v] = row
print(json.dumps(out, indent=1))
