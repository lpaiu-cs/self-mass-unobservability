#!/usr/bin/env python3
"""10.8 T2 execution against the PATCHED build:
  Stage A: A=0 integrity gate (patched baseline must reproduce res0 to <1e-3 us)
  Stage B: time-convention verification (drive at 2pi/600 d, response frequency
           must sit at the drive frequency; auto-detects days vs timedays)
  Stage C: production dynamic response columns at {Omega_in, Omega_out,
           Omega_dif} x {cos, sin} quadratures, central difference in A.
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
print('patched module:', pfi.__file__, flush=True)

base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)
errs = base['errs'].astype(np.float64)
t = base['toas'].astype(np.float64)
sw = 1.0/errs
N = res0.size

PARFILE = 'parfile-planetGR-max-bestfit'
TIMFILE = '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'
TIMEDAYS = 2.4077445558945513991e+01     # integrator time scale (days), from construction dump
OMS = {'in': 3.856136695791377, 'out': 0.019199654300698213, 'dif': 3.836937041490679}
WRAP_GUARD = 1000.0                      # us; P_spin/2 = 1366 us

def run(A, w_days, ph, wfac):
    os.environ['SEPDYN_A'] = repr(A)
    os.environ['SEPDYN_W'] = repr(w_days * wfac)
    os.environ['SEPDYN_PH'] = repr(ph)
    t0 = time.perf_counter()
    fit = pfi.PyFittriple(PARFILE, TIMFILE)
    fit.Compute_lnposterior()
    r = fit.Get_time_residuals().copy()
    del fit
    print('  run A=%+.3e w=%.6g ph=%.3f (fac %.4g): %.0fs' %
          (A, w_days, ph, wfac, time.perf_counter()-t0), flush=True)
    return r

def clear_env():
    for k in ('SEPDYN_A', 'SEPDYN_W', 'SEPDYN_PH'):
        os.environ.pop(k, None)

def amp_at(col, w_days):
    X = np.column_stack([sw*np.cos(w_days*t), sw*np.sin(w_days*t)])
    c, _, _, _ = np.linalg.lstsq(X, sw*col, rcond=None)
    return float(np.hypot(*c))

out = {}

# ---- Stage A: A=0 integrity gate ----
clear_env()
rA0 = run(0.0, 0.0, 0.0, 1.0)  # env set to zeros (equivalent to unset)
gate = float(np.abs(rA0 - res0).max())
out['gate_A0'] = {'max_abs_diff_us': gate, 'pass': bool(gate < 1e-3)}
print('GATE A=0: max|diff| = %.3g us -> %s' % (gate, 'PASS' if gate < 1e-3 else 'FAIL'), flush=True)
if gate >= 1e-3:
    json.dump(out, open('sep_t2_results.json', 'w'), indent=1)
    raise SystemExit('STOP RULE: patched baseline does not reproduce')

# ---- Stage B: time-convention verification ----
WV = 2.0*np.pi/600.0                     # rad/day, clean band away from orbital lines
AV = 1e-8
rv = run(AV, WV, 0.0, TIMEDAYS)
colv = (rv - res0)/AV
cands = {'assumed_timedays': WV, 'if_t_in_days': WV*TIMEDAYS, 'if_double_scaled': WV/TIMEDAYS}
amps = {k: amp_at(colv, w) for k, w in cands.items()}
winner = max(amps, key=amps.get)
dom = amps[winner]/max(v for k, v in amps.items() if k != winner)
out['convention'] = {'candidate_amps': amps, 'winner': winner, 'dominance': dom,
                     'maxdev_us': float(np.abs(rv-res0).max()*1.0)}
print('CONVENTION: %s (dominance %.1fx); amps %s' % (winner, dom, amps), flush=True)
if winner == 'assumed_timedays' and dom >= 10:
    WFAC = TIMEDAYS
elif winner == 'if_t_in_days' and dom >= 10:
    WFAC = 1.0
    print('  -> switching to w-fac 1.0 (integrator time is in days)', flush=True)
else:
    json.dump(out, open('sep_t2_results.json', 'w'), indent=1)
    raise SystemExit('STOP RULE: ambiguous time convention (dominance %.1f)' % dom)
out['convention']['wfac_used'] = WFAC

# ---- Stage C: production columns ----
cols = {}
meta = {}
for name, w in OMS.items():
    for tag, ph in (('c', 0.0), ('s', np.pi/2.0)):
        A = 1e-8
        for attempt in range(3):
            rp = run(A, w, ph, WFAC)
            maxdev = float(np.abs(rp - res0).max())
            if maxdev < 0.02 and attempt < 2:
                print('%s_%s: A=%g maxdev=%.3g too small -> x30' % (name, tag, A, maxdev), flush=True)
                A *= 30.0
                continue
            if maxdev > 800.0 and attempt < 2:
                print('%s_%s: A=%g maxdev=%.3g too large -> /30' % (name, tag, A, maxdev), flush=True)
                A /= 30.0
                continue
            break
        if maxdev > WRAP_GUARD:
            print('WARN %s_%s: maxdev %.3g near wrap threshold' % (name, tag, maxdev), flush=True)
        rm = run(-A, w, ph, WFAC)
        col = (rp - rm)/(2.0*A)
        cols['col_%s_%s' % (name, tag)] = col
        meta['%s_%s' % (name, tag)] = {'A': A, 'maxdev_us': maxdev,
                                       'amp_at_drive': amp_at(col, w)}
        print('%s_%s: A=%g maxdev=%.3f us, response amp at drive = %.4g us' %
              (name, tag, A, maxdev, meta['%s_%s' % (name, tag)]['amp_at_drive']), flush=True)

np.savez('sep_dynamic_columns.npz', **cols)
out['columns'] = meta
with open('sep_t2_results.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('T2_DONE', flush=True)
