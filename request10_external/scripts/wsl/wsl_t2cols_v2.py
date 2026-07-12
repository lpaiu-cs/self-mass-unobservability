#!/usr/bin/env python3
"""10.8 T2 v2: empirical time-mapping calibration by spectral scan, then
production dynamic response columns.

v1's single-sinusoid convention test was defeated by the broadband/secular
envelope of the response (single-line fits captured ~2.6% of a 17 us
structure). v2: two probe drives at different frequencies; cubic detrend;
amplitude scan over a frequency grid; the response peak must MOVE in
proportion to the drive -> empirical mapping g = f_response/W_env, robust to
ANY time convention. Production uses W_env = Omega/g.
"""
import json
import os
import sys
import time

import numpy as np

SEPDYN_SRC = os.path.expanduser('~/work/nutimo_pilot/nutimo_sepdyn/src')
os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
sys.path.insert(0, SEPDYN_SRC)
import python_Fittriple_interface as pfi
assert 'nutimo_sepdyn' in pfi.__file__, pfi.__file__

base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)
errs = base['errs'].astype(np.float64)
t = base['toas'].astype(np.float64)
sw = 1.0/errs
N = res0.size

PARFILE = 'parfile-planetGR-max-bestfit'
TIMFILE = '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'
TIMEDAYS = 2.4077445558945513991e+01
OMS = {'in': 3.856136695791377, 'out': 0.019199654300698213, 'dif': 3.836937041490679}
ORBITAL_LINES = [OMS['in'], OMS['out'], OMS['dif'], 2*OMS['in'], 2*OMS['out'],
                 OMS['in']+OMS['out'], 2*OMS['dif']]
WRAP_GUARD = 1000.0

def run_env(A, W_env, ph):
    os.environ['SEPDYN_A'] = repr(A)
    os.environ['SEPDYN_W'] = repr(W_env)
    os.environ['SEPDYN_PH'] = repr(ph)
    t0 = time.perf_counter()
    fit = pfi.PyFittriple(PARFILE, TIMFILE)
    fit.Compute_lnposterior()
    r = fit.Get_time_residuals().copy()
    del fit
    print('  run A=%+.3e W_env=%.6g ph=%.3f: %.0fs' %
          (A, W_env, ph, time.perf_counter()-t0), flush=True)
    return r

# cubic detrend basis (weighted)
tc = (t - t.mean())/(t.max() - t.min())
POLY = sw[:, None]*np.column_stack([np.ones(N), tc, tc**2, tc**3])
QP, _ = np.linalg.qr(POLY)

def detrend(col):
    x = sw*col
    return x - QP @ (QP.T @ x)

def spectrum(xw, fgrid):
    amps = np.empty(len(fgrid))
    for i, w in enumerate(fgrid):
        X = np.column_stack([sw*np.cos(w*t), sw*np.sin(w*t)])
        X = X - QP @ (QP.T @ X)
        c, _, _, _ = np.linalg.lstsq(X, xw, rcond=None)
        amps[i] = np.hypot(*c)
    return amps

def offline_peak(fgrid, amps):
    mask = np.ones(len(fgrid), bool)
    for L in ORBITAL_LINES:
        mask &= np.abs(fgrid/L - 1.0) > 0.02
    i = int(np.argmax(np.where(mask, amps, 0.0)))
    return fgrid[i], amps[i]

out = {}
A0 = 1e-8
W1_days, W2_days = 2*np.pi/600.0, 2*np.pi/150.0
FGRID = np.geomspace(2*np.pi/30000.0, 3.0, 420)

# ---- probes ----
probes = {}
for tag, wd in (('p600', W1_days), ('p150', W2_days)):
    r = run_env(A0, wd*TIMEDAYS, 0.0)
    col = (r - res0)/A0
    xw = detrend(col)
    amps = spectrum(xw, FGRID)
    fpk, apk = offline_peak(FGRID, amps)
    probes[tag] = {'W_env': wd*TIMEDAYS, 'peak_rad_per_day': float(fpk),
                   'peak_amp': float(apk), 'maxdev_us': float(np.abs(r-res0).max()),
                   'g': float(fpk/(wd*TIMEDAYS))}
    print('%s: peak at %.6g rad/day (amp %.4g), g = f/W_env = %.6g' %
          (tag, fpk, apk, probes[tag]['g']), flush=True)
out['probes'] = probes

g1, g2 = probes['p600']['g'], probes['p150']['g']
gbar = float(np.sqrt(g1*g2))
consistent = abs(g1/g2 - 1.0) < 0.15
out['mapping'] = {'g1': g1, 'g2': g2, 'g_used': gbar, 'consistent': bool(consistent),
                  'note': 'g ~ 1/TIMEDAYS means integrator t is in days; g ~ 1 means timedays units'}
print('MAPPING: g1=%.6g g2=%.6g consistent=%s -> g=%.6g (1/TD=%.6g)' %
      (g1, g2, consistent, gbar, 1.0/TIMEDAYS), flush=True)
if not consistent:
    json.dump(out, open('sep_t2_results_v2.json', 'w'), indent=1)
    raise SystemExit('STOP RULE: probe peaks do not move proportionally with the drive')

# ---- production ----
cols = {}
meta = {}
for name, w in OMS.items():
    W_env = w/gbar
    for tag, ph in (('c', 0.0), ('s', np.pi/2.0)):
        A = 1e-8
        for attempt in range(3):
            rp = run_env(A, W_env, ph)
            maxdev = float(np.abs(rp - res0).max())
            if maxdev < 0.02 and attempt < 2:
                A *= 30.0
                print('%s_%s: maxdev %.3g too small -> A=%g' % (name, tag, maxdev, A), flush=True)
                continue
            if maxdev > 800.0 and attempt < 2:
                A /= 30.0
                print('%s_%s: maxdev %.3g too large -> A=%g' % (name, tag, maxdev, A), flush=True)
                continue
            break
        if maxdev > WRAP_GUARD:
            print('WARN %s_%s: maxdev %.3g near wrap' % (name, tag, maxdev), flush=True)
        rm = run_env(-A, W_env, ph)
        col = (rp - rm)/(2.0*A)
        cols['col_%s_%s' % (name, tag)] = col
        xw = detrend(col)
        X = np.column_stack([sw*np.cos(w*t), sw*np.sin(w*t)])
        X = X - QP @ (QP.T @ X)
        c, _, _, _ = np.linalg.lstsq(X, xw, rcond=None)
        meta['%s_%s' % (name, tag)] = {'A': A, 'W_env': W_env, 'maxdev_us': maxdev,
                                       'amp_at_drive_detrended': float(np.hypot(*c))}
        print('%s_%s: A=%g maxdev=%.3f us, detrended amp at drive = %.4g' %
              (name, tag, A, maxdev, meta['%s_%s' % (name, tag)]['amp_at_drive_detrended']), flush=True)

np.savez('sep_dynamic_columns.npz', **cols)
out['columns'] = meta
with open('sep_t2_results_v2.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('T2V2_DONE', flush=True)
