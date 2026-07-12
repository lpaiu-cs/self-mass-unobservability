#!/usr/bin/env python3
"""Request 10.8 Gate F0: static SEP_D response column j_Delta = d res/d SEP_D.

Central difference at h = 1e-6 (parfile scale), adapted x10 if maxdev < 0.02 us
or /10 if > 30 us; half-step linearity check < 5%. SEP_D is full-parameter
index 20 (Parameters.cpp: absparamnumber == 20), currently 0 and unfitted --
so we shift it via the FULL parameter vector interface if available, else via
a parfile edit. jacobian_runner used Set_fitted_parameter_relativeshifts over
the 28 FITTED parameters only; SEP_D is not among them, so this script edits
the parfile value directly (construction per value: 2 constructions per
difference is too slow) -- preferred: check for a full-parameter shift API
first.
"""
import json
import os
import sys
import time

import numpy as np

os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
sys.path.insert(0, os.path.expanduser('~/work/nutimo_pilot/nutimo/src'))
import python_Fittriple_interface as pfi

print('interface members with "param"/"shift":',
      [m for m in dir(pfi.PyFittriple) if 'aram' in m or 'hift' in m], flush=True)

base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)
errs = base['errs'].astype(np.float64)
w = 1.0 / errs

PARFILE = 'parfile-planetGR-max-bestfit'
TIMFILE = '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'

def read_parfile():
    return open(PARFILE).read().splitlines()

def write_parfile_sep(value, path):
    lines = read_parfile()
    out = []
    replaced = False
    for ln in lines:
        tok = ln.split()
        if tok and tok[0] == 'SEP_D':
            parts = tok[:]
            parts[1] = '%.16e' % value
            out.append(' '.join(parts))
            replaced = True
        else:
            out.append(ln)
    if not replaced:
        raise SystemExit('SEP_D line not found in parfile')
    with open(path, 'w') as fh:
        fh.write('\n'.join(out) + '\n')

def res_at_sep(value):
    path = 'parfile-sep-%s' % ('p' if value > 0 else ('m' if value < 0 else 'z'))
    write_parfile_sep(value, path)
    t0 = time.perf_counter()
    fit = pfi.PyFittriple(path, TIMFILE)
    fit.Compute_lnposterior()
    r = fit.Get_time_residuals().copy()
    print('  SEP_D=%+.3e: constructed+evaluated in %.0fs' % (value, time.perf_counter()-t0), flush=True)
    del fit
    return r

report = {}
h = 1e-6
for attempt in range(3):
    res_p = res_at_sep(+h)
    maxdev = float(np.abs(res_p - res0).max())
    if maxdev < 0.02 and attempt < 2:
        print('h=%g maxdev=%.3g us too small -> x10' % (h, maxdev), flush=True)
        h *= 10.0
        continue
    if maxdev > 30.0 and attempt < 2:
        print('h=%g maxdev=%.3g us too large -> /10' % (h, maxdev), flush=True)
        h /= 10.0
        continue
    break
res_m = res_at_sep(-h)
dcol = (res_p - res_m) / (2.0 * h)

# half-step linearity
res_h2 = res_at_sep(+0.5*h)
d_act = res_h2 - res0
d_lin = 0.5 * h * dcol
lin_dev = float(np.linalg.norm(w*(d_act - d_lin)) / np.linalg.norm(w*d_lin))

os.makedirs('jac_sep', exist_ok=True)
np.savez('jac_sep/col_SEP_D.npz', dcol=dcol, h_abs=h,
         maxdev_us=maxdev, halfstep_lin_dev=lin_dev)
report = {'h_abs': h, 'maxdev_us': maxdev, 'halfstep_lin_dev': lin_dev,
          'pass_F0': bool(lin_dev < 0.05),
          'col_wnorm': float(np.linalg.norm(w*dcol)),
          'note': 'j_Delta = dres/dSEP_D, us per unit Delta; via parfile edit + reconstruction'}
with open('jac_sep/sep_F0_report.json', 'w') as fh:
    json.dump(report, fh, indent=1)
print('F0: h=%g maxdev=%.3f us  halfstep_lin_dev=%.4f  -> %s'
      % (h, maxdev, lin_dev, 'PASS' if report['pass_F0'] else 'FAIL'), flush=True)
print('SEPF0_DONE')
