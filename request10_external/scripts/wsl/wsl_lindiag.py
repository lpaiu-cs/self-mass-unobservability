#!/usr/bin/env python3
"""Request 10.7b linearization-failure diagnosis: amplitude scaling + direction split.

One construction, six evaluations:
  f025, f050          -- dtheta scaled by 0.25, 0.50 (quadratic-scaling test)
  no_delta_i          -- dtheta with the delta_i component zeroed
  no_planet           -- dtheta with the 7 *_extra1 components zeroed
  only_delta_i        -- pure delta_i displacement
  only_planet         -- pure *_extra1 displacement
Deviation metric per case: ||res(p0+d) - res0 - J d||_W / ||J d||_W.
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
J = np.load('artifacts/finite_jacobian.npy')
dtheta = np.load('jointfit_dtheta52.npy')
sc = np.array([float(scales[p]) for p in fmap])
w = 1.0 / errs

idx_delta_i = fitted_names.index('delta_i')
idx_planet = [j for j, n in enumerate(fitted_names) if n.endswith('_extra1')]

cases = {}
d = dtheta.copy(); cases['f025'] = 0.25 * d
cases['f050'] = 0.50 * d
v = d.copy(); v[idx_delta_i] = 0.0; cases['no_delta_i'] = v
v = d.copy(); v[idx_planet] = 0.0; cases['no_planet'] = v
v = np.zeros_like(d); v[idx_delta_i] = d[idx_delta_i]; cases['only_delta_i'] = v
v = np.zeros_like(d); v[idx_planet] = d[idx_planet]; cases['only_planet'] = v

t0 = time.perf_counter()
fit = PyFittriple('parfile-planetGR-max-bestfit',
                  '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
lnp0 = fit.Compute_lnposterior()
print('constructed in %.0fs lnp0=%.6f' % (time.perf_counter() - t0, lnp0), flush=True)

out = {'reference_full_dtheta_deviation': 5.9315, 'cases': {}}
for label, dv in cases.items():
    t1 = time.perf_counter()
    fit.Set_fitted_parameter_relativeshifts(dv / sc)
    lnp = fit.Compute_lnposterior()
    res1 = fit.Get_time_residuals().copy()
    d_act = res1 - res0
    d_lin = J @ dv
    num = float(np.linalg.norm(w * (d_act - d_lin)))
    den = float(np.linalg.norm(w * d_lin))
    out['cases'][label] = {
        'relative_deviation': num / den,
        'W_norm_linear': den, 'W_norm_actual_minus_linear': num,
        'lnp': float(lnp), 'seconds': time.perf_counter() - t1,
    }
    print('%-14s dev=%.4f  |lin|_W=%.2f  lnp=%.4f  (%.0fs)'
          % (label, num / den, den, lnp, time.perf_counter() - t1), flush=True)

with open('jointfit_lindiag.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('LINDIAG_DONE')
