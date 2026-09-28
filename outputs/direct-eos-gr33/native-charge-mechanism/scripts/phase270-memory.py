"""Phase270 diagnostic: how much of the 4x endpoint charge comes from layers the pulse has already left (frozen
displacement), and how frozen are they? Uses the 4x T readout source (baryon_g stage values) and the phase-268
per-sub-cell depth bands (charge per 4x cell at T).
"""
import json
import numpy as np
R = '/home/lpaiu/work/native-refined268-runtime/readout268-quad64-work/'
z = np.load(R + 'gr/source-64.npz'); t, M, x = z['t'], z['baryon_g'], z['drive_x']; D = float(z['drive_duration']); C = 2.99792458e10
bands = {r['band']: r['charge'] for r in json.load(open(R + 'depth-bands.json'))['rows']}
total = bands['all']; T = t[-1]
phase_T = (T + x/C)/D  # pulse phase at each cell centre at T (>1: the pulse has left the cell)
i_ref = int(np.argmin(abs(t - (T - 3.0e-4))))  # about 0.3 ms before T
passed, sweeping, rows = 0., 0., []
for c in range(8, 43):
    q = bands[f'cells:{c}-{c}']
    if phase_T[c] > 1: passed += q
    else: sweeping += q
    drift = (M[-1, c] - M[i_ref, c])/abs(M[-1, c]) if M[-1, c] else 0.
    rows.append((c, phase_T[c], q, drift))
print('T', T, 'total', total, 'reference time', t[i_ref])
print('charge from cells the pulse has left (phase>1): %.4e (%.1f%%)   still being swept: %.4e (%.1f%%)' % (passed, 100*passed/total, sweeping, 100*sweeping/total))
for c, p, q, dr in rows:
    if abs(q) > 1e-3*abs(total): print('cell %2d phase(T) %.3f  charge %+.3e (%5.1f%%)  dM drift over last 0.3 ms %+.2e' % (c, p, q, 100*q/total, dr))
