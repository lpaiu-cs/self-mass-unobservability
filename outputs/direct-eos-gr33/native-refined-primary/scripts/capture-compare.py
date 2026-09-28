"""Captured stage data of the t=0 validation run against the accepted recovered photon history (same stages)."""
import json
import numpy as np
root = '/home/lpaiu/work/native-retained-tail-runtime/'
rec = np.load(root + 'native-full-captured244-work/recovered-64.npz'); out = {}
for i in range(8):
    c = np.load(root + f'primary267-validate2-work/captures/captured-64-{i:03d}.npz')
    assert c['time'] == rec['times'][i] and c['weight'] == rec['weights'][i], i
    row = {}
    for key in ['photon_moments', 'radial_ports', 'collision_rates']:
        a, b = np.asarray(c[key], float), np.asarray(rec[key][i], float); s = np.max(abs(b))
        row[key] = float(np.max(abs(a - b))/s) if s else float(np.max(abs(a)))
    out[i] = row
print(json.dumps(dict(stages=out, maximum=max(max(r.values()) for r in out.values())), indent=1))
