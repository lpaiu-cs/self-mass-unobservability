#!/usr/bin/env python3
"""Request 10.7b linearization check (pre-registered): one full Nutimo evaluation
at the fitted nuisance displacement of the tau=52 d joint-fit configuration.

Pass criterion: || res(p0+dtheta) - res0 - J dtheta ||_W / || J dtheta ||_W < 0.05
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
J = np.load('artifacts/finite_jacobian.npy')
dtheta = np.load('jointfit_dtheta52.npy')          # parameter units, len 28
sh = dtheta / np.array([float(scales[p]) for p in fmap])

t0 = time.perf_counter()
fit = PyFittriple('parfile-planetGR-max-bestfit',
                  '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
lnp0 = fit.Compute_lnposterior()
res_base = fit.Get_time_residuals().copy()
t1 = time.perf_counter()
base_repro = float(np.abs(res_base - res0).max())
print('constructed in %.0fs; baseline reproduction max|d| = %.3g us' % (t1 - t0, base_repro), flush=True)

fit.Set_fitted_parameter_relativeshifts(sh)
lnp1 = fit.Compute_lnposterior()
res1 = fit.Get_time_residuals().copy()
t2 = time.perf_counter()

w = 1.0 / errs
d_actual = res1 - res0
d_lin = J @ dtheta
num = float(np.linalg.norm(w * (d_actual - d_lin)))
den = float(np.linalg.norm(w * d_lin))
ratio = num / den
out = {
    'preregistration': 'notes/REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md (linearization check)',
    'configuration': 'tau = 52 d joint-fit nuisance displacement (jointfit_dtheta52.npy)',
    'dtheta_max_abs_over_step': float(np.max(np.abs(dtheta) / np.load('steps.npz', allow_pickle=True)['abs_steps'])),
    'baseline_reproduction_max_abs_us': base_repro,
    'lnp0': float(lnp0), 'lnp_shifted': float(lnp1),
    'W_norm_actual_minus_linear': num,
    'W_norm_linear_prediction': den,
    'relative_deviation': ratio,
    'pass_criterion': 'relative_deviation < 0.05',
    'pass': bool(ratio < 0.05),
    'eval_seconds': t2 - t1,
    'shift_convention': 'sh = dtheta_param / scales[fmap] via Set_fitted_parameter_relativeshifts (non-cumulative)',
}
with open('jointfit_linearization_check.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
print('LINCHECK_DONE pass=%s ratio=%.4f' % (out['pass'], ratio), flush=True)
