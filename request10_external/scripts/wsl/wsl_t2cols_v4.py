#!/usr/bin/env python3
"""10.8 T2 v4: quasi-static RAMP calibration (decisive), then production.

Calibration: drive at ultra-low frequency (period 30000 d), phase pi/2:
Delta(t) = A cos(W_env t_int + pi/2) = -A sin(W_env t_days / U), U = time unit
in days. Deep in the adiabatic regime the response factorizes on the MEASURED
static column j_Delta (gate F0): col(t) ~= -j_Delta(t) * sin(W_env t/U).
The two unit hypotheses U=1 (days) vs U=24.077 (timedays) give templates
differing by ~24x in amplitude and in shape; a weighted fit picks the winner
and the fitted amplitude a_hat ~= 1 is an end-to-end normalization check of
the patch against the static F0 measurement.
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
WRAP_GUARD = 1000.0

jD = np.load('jac_sep/col_SEP_D.npz')['dcol'].astype(np.float64)

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

out = {}

# ---- ramp calibration ----
A0 = 1e-8
W_ENV = 2.0*np.pi/30000.0
r = run_env(A0, W_ENV, np.pi/2.0)
col = (r - res0)/A0
np.savez('sep_ramp_probe.npz', col=col, A=A0, W_env=W_ENV)
fits = {}
for label, U in (('U_days', 1.0), ('U_timedays', TIMEDAYS)):
    tpl = -jD*np.sin(W_ENV*t/U)
    xw, yw = sw*tpl, sw*col
    a_hat = float(xw @ yw / (xw @ xw))
    resid = yw - a_hat*xw
    R2 = float(1.0 - (resid @ resid)/(yw @ yw))
    fits[label] = {'U': U, 'a_hat': a_hat, 'R2': R2,
                   'overlap_cos': float((xw @ yw)/np.linalg.norm(xw)/np.linalg.norm(yw))}
    print('%s: a_hat=%.4f  R2=%.5f  cos=%.5f' % (label, a_hat, R2, fits[label]['overlap_cos']),
          flush=True)
out['ramp'] = {'maxdev_us': float(np.abs(r-res0).max()), 'fits': fits}

best = max(fits, key=lambda k: fits[k]['R2'])
ok = (fits[best]['R2'] > 0.8 and fits[best]['overlap_cos'] > 0.9 and
      0.5 < abs(fits[best]['a_hat']) < 2.0)
out['ramp']['winner'] = best
out['ramp']['pass'] = bool(ok)
print('RAMP: winner %s (R2 %.4f vs %.4f), a_hat=%.3f -> %s' %
      (best, fits[best]['R2'], min(f['R2'] for f in fits.values()),
       fits[best]['a_hat'], 'CONFIRMED' if ok else 'AMBIGUOUS'), flush=True)
if not ok:
    json.dump(out, open('sep_t2_results_v4.json', 'w'), indent=1)
    raise SystemExit('STOP RULE: ramp calibration ambiguous')
U_UNIT = fits[best]['U']

# ---- production ----
cols = {}
meta = {}
for name, w in OMS.items():
    W_env = w*U_UNIT
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
        cols['col_%s_%s' % (name, tag)] = (rp - rm)/(2.0*A)
        meta['%s_%s' % (name, tag)] = {'A': A, 'W_env': W_env, 'maxdev_us': maxdev}
        print('%s_%s: A=%g maxdev=%.3f us' % (name, tag, A, maxdev), flush=True)

np.savez('sep_dynamic_columns.npz', **cols)
out['columns'] = meta
out['U_unit_days'] = U_UNIT
with open('sep_t2_results_v4.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('T2V4_DONE', flush=True)
