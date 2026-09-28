"""Phase270 diagnostic: is the charge-carrying baryon perturbation a mass redistribution (displacement) or a local
instantaneous change? Uses the 4x T readout source (stage values of baryon_g, grams per cell, amplitude-scaled).

Prints, over the stored stage times: the total sum over cells, the sum of |dM|, the ratio, and for the dominant 4x
sub-cells (original cells 10-12) the time series next to the incident pulse phase at the cell center.
"""
import numpy as np
z = np.load('/home/lpaiu/work/native-refined268-runtime/readout268-quad64-work/gr/source-64.npz')
t, M, r = z['t'], z['baryon_g'], z['radius']; x = z['drive_x']; D = float(z['drive_duration']); C = 2.99792458e10
print('stages', len(t), 'cells', M.shape[1], 'T', t[-1])
tot = M.sum(1); ab = np.abs(M).sum(1)
for i in list(range(0, len(t), 12)) + [len(t) - 1]:
    print('t %.4e  sum dM %+.3e  sum|dM| %.3e  ratio %+.2e' % (t[i], tot[i], ab[i], tot[i]/ab[i] if ab[i] else 0))
depth = (z['edges'][-1] - z['edges'])/1e5
cells = list(range(16, 32))  # 4x sub-cells of original cells 10-13
print('cells', cells, 'depth km', [round(float((depth[c] + depth[c+1])/2), 1) for c in cells[::4]])
for i in list(range(40, len(t), 10)) + [len(t) - 1]:
    phase = (t[i] + x[cells]/C)/D
    print('t %.4e dM %s | pulse phase %s' % (t[i], np.array2string(M[i, cells[::4]], precision=3), np.array2string(phase[::4], precision=3)))
