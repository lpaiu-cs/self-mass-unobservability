"""Phase267 driver validation: t=0 run records against the accepted primary history at the same record times.

Usage: python3 validate-compare.py <new run npz> <accepted npz>
Reports, per record k shared by both, the max relative difference of each history array (normalised by the
accepted array's max at that record), and the overall maximum. Driver acceptance: <= 1e-6 (progress 3).
"""
import json, sys
import numpy as np
a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
ta, tb = a['t'], b['t']; rows = {}; worst = 0.
for k, t in enumerate(ta):
    j = int(np.argmin(abs(tb - t))); assert abs(tb[j] - t) <= 1e-18, (k, t, tb[j])
    row = {}
    for key in ['moments', 'photon_history_scaled_occupation', 'material_history', 'radial_ports', 'collision_transfer']:
        x, y = np.asarray(a[key][k], float), np.asarray(b[key][j], float)
        scale = np.max(abs(y)); row[key] = float(np.max(abs(x - y)) / scale) if scale else float(np.max(abs(x)))
    if k: worst = max(worst, max(row.values()))
    rows[f'{k}:{t:.6e}'] = row
print(json.dumps(dict(records=rows, maximum_relative_after_t0=worst, accepted=worst <= 1e-6), indent=1))
