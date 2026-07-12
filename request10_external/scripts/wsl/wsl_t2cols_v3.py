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

# structure basis: cubic poly + growing orbital-line pairs {1,t}x{cos,sin}(W t)
# (the dominant response to an oscillating pulsar-pair G is orbital-phase
# modulation, i.e. t*sin(W_orb t)-type terms; remove them before line search)
tc = (t - t.mean())/(t.max() - t.min())
bases = [np.ones(N), tc, tc**2, tc**3]
for _w in ORBITAL_LINES[:3]:
    for trig in (np.cos, np.sin):
        bases += [trig(_w*t), tc*trig(_w*t)]
POLY = sw[:, None]*np.column_stack(bases)
QP, _ = np.linalg.qr(POLY)

def detrend(col):
    x = sw*col
    return x - QP @ (QP.T @ x)

def line_amp(xw, w):
    """norm of the projection of xw onto the (structure-projected) line basis
    at frequency w -- bounded by ||xw||, immune to ill-conditioning."""
    X = np.column_stack([sw*np.cos(w*t), sw*np.sin(w*t)])
    X = X - QP @ (QP.T @ X)
    Q, _ = np.linalg.qr(X)
    return float(np.linalg.norm(Q.T @ xw))

out = {}
A0 = 1e-8
W1_days, W2_days = 2*np.pi/600.0, 2*np.pi/150.0

# ---- probes: two physical hypotheses ----
#   H_timedays: integrator t in units of TIMEDAYS -> response at w_days when W_env = w*TD
#   H_days:     integrator t in days              -> response at w_days when W_env = w
probes = {}
probe_cols = {}
for tag, wd in (('p600', W1_days), ('p150', W2_days)):
    r = run_env(A0, wd*TIMEDAYS, 0.0)
    col = (r - res0)/A0
    probe_cols['probe_%s' % tag] = col
    xw = detrend(col)
    a_td = line_amp(xw, wd)             # response AT w_days -> supports H_timedays
    a_dy = line_amp(xw, wd*TIMEDAYS)    # response at w*TD    -> supports H_days
    probes[tag] = {'W_env': wd*TIMEDAYS, 'amp_H_timedays': a_td, 'amp_H_days': a_dy,
                   'maxdev_us': float(np.abs(r-res0).max()),
                   'norm_after_structure_removal': float(np.linalg.norm(xw))}
    print('%s: maxdev=%.3g us | amp@w(H_td)=%.4g vs amp@w*TD(H_days)=%.4g (struct-removed norm %.4g)' %
          (tag, probes[tag]['maxdev_us'], a_td, a_dy,
           probes[tag]['norm_after_structure_removal']), flush=True)
np.savez('sep_probe_columns.npz', **probe_cols)
out['probes'] = probes

r1 = probes['p600']['amp_H_timedays']/max(probes['p600']['amp_H_days'], 1e-30)
r2 = probes['p150']['amp_H_timedays']/max(probes['p150']['amp_H_days'], 1e-30)
if r1 >= 5 and r2 >= 5:
    WFAC = TIMEDAYS; hyp = 'H_timedays'
elif r1 <= 0.2 and r2 <= 0.2:
    WFAC = 1.0; hyp = 'H_days'
else:
    out['mapping'] = {'ratio_p600': r1, 'ratio_p150': r2, 'decision': 'ambiguous'}
    json.dump(out, open('sep_t2_results_v2.json', 'w'), indent=1)
    raise SystemExit('STOP RULE: hypothesis ratios ambiguous (%.3g, %.3g)' % (r1, r2))
out['mapping'] = {'ratio_p600': r1, 'ratio_p150': r2, 'hypothesis': hyp, 'wfac': WFAC}
gbar = 1.0/WFAC if WFAC != 0 else 0.0
print('MAPPING: %s confirmed (ratios %.3g, %.3g) -> W_env = w * %.6g' %
      (hyp, r1, r2, WFAC), flush=True)

# ---- production ----
cols = {}
meta = {}
for name, w in OMS.items():
    W_env = w*WFAC
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
