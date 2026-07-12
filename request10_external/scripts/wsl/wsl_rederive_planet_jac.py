#!/usr/bin/env python3
"""Request 10.7b: re-derive the 7 planet-block (*_extra1) Jacobian columns with
steps small enough to be true derivatives.

Why: the 10.7a steps were 0.3 * sigma_MCMC, but the planet elements are so
weakly constrained that those steps are huge fractions of the elements
(tasc_extra1 step = 8739 d = 2.8 planet periods; acosi step = 148% of value;
oman step = 1 rad). The resulting columns are secants over nonperturbative
arcs, not tangent vectors -- confirmed by jointfit_lindiag.json (only_planet
rel. deviation 7.9, oscillatory in amplitude; all other 21 params linear to
0.3%).

New steps: ~0.1% of the element scale (angle/time steps ~1e-3 of the planet
period), adapted once by x10 (up) if maxdev < 0.02 us or /10 (down) if
maxdev > 30 us. Central difference, checkpointed to jac_v2/.
Also validates linearity per column: eval at +h/2 must match 0.5*(column
prediction) to 5%.
"""
import json
import os
import sys
import time

import numpy as np

os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
sys.path.insert(0, os.path.expanduser('~/work/nutimo_pilot/nutimo/src'))
from python_Fittriple_interface import PyFittriple

base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)
errs = base['errs'].astype(np.float64)
scales = base['scales']
fmap = list(base['fmap'])
names = list(base['names'])
fitted_names = [names[p] for p in fmap]
nfit = len(fmap)
sc = np.array([float(scales[p]) for p in fmap])
w = 1.0 / errs

H_NEW = {
    'eta_extra1': 1e-3,
    'kappa_extra1': 1e-3,
    'asini_extra1': 1.5,
    'acosi_extra1': 3.2,
    'tasc_extra1': 3.0,
    'P_extra1': 0.3,
    'oman_extra1': 1e-3,
}

os.makedirs('jac_v2', exist_ok=True)
t0 = time.perf_counter()
fit = PyFittriple('parfile-planetGR-max-bestfit',
                  '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
lnp0 = fit.Compute_lnposterior()
print('constructed in %.0fs lnp0=%.6f' % (time.perf_counter() - t0, lnp0), flush=True)

def eval_shift(sh_vec):
    fit.Set_fitted_parameter_relativeshifts(sh_vec)
    fit.Compute_lnposterior()
    return fit.Get_time_residuals().copy()

report = {}
for pname, h0 in H_NEW.items():
    j = fitted_names.index(pname)
    h = float(h0)
    for attempt in range(3):
        sh = np.zeros(nfit)
        t1 = time.perf_counter()
        sh[j] = +h / sc[j]
        res_p = eval_shift(sh)
        sh[j] = -h / sc[j]
        res_m = eval_shift(sh)
        maxdev = float(np.abs(res_p - res0).max())
        if maxdev < 0.02 and attempt < 2:
            print('%s: h=%g maxdev=%.3g us too small -> x10' % (pname, h, maxdev), flush=True)
            h *= 10.0
            continue
        if maxdev > 30.0 and attempt < 2:
            print('%s: h=%g maxdev=%.3g us too large -> /10' % (pname, h, maxdev), flush=True)
            h /= 10.0
            continue
        break
    dcol = (res_p - res_m) / (2.0 * h)
    # linearity validation at h/2: actual vs 0.5*h*dcol prediction
    sh = np.zeros(nfit); sh[j] = +0.5 * h / sc[j]
    res_h2 = eval_shift(sh)
    d_act = res_h2 - res0
    d_lin = 0.5 * h * dcol
    lin_dev = float(np.linalg.norm(w * (d_act - d_lin)) / np.linalg.norm(w * d_lin))
    np.savez('jac_v2/col_%02d_%s.npz' % (j, pname), dcol=dcol, h_abs=h,
             h_rel=h / sc[j], maxdev_us=maxdev, halfstep_lin_dev=lin_dev)
    report[pname] = {'h_abs': h, 'maxdev_us': maxdev, 'halfstep_lin_dev': lin_dev,
                     'seconds': time.perf_counter() - t1}
    print('%-14s h=%-8g maxdev=%8.3f us  halfstep_lin_dev=%.4f  (%.0fs)'
          % (pname, h, maxdev, lin_dev, time.perf_counter() - t1), flush=True)

with open('jac_v2/rederive_report.json', 'w') as fh:
    json.dump(report, fh, indent=1)
print('REDERIVE_DONE')
