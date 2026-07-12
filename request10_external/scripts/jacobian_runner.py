#!/usr/bin/env python3
"""Request 10.7a external Nutimo pilot - finite Jacobian worker.

Computes central-difference residual derivatives d(res)/d(param) for the
standard fitted parameter set of parfile-planetGR-max-bestfit, one column
per fitted parameter, checkpointed to jac/col_*.npz.

Usage: jacobian_runner.py <worker_id> <n_workers>
Steps: absolute step per fitted parameter from steps.npz (FRAC * MCMC
sigma; the run of record used FRAC = 1.0 -- see prep_steps.py, 10.8f C11).
Column definition: J[:,j] = (res(p_j + h) - res(p_j - h)) / (2 h_abs)  [us / unit]
"""
import sys, os, time
import numpy as np

os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
sys.path.insert(0, os.path.expanduser('~/work/nutimo_pilot/nutimo/src'))
from python_Fittriple_interface import PyFittriple

worker, nworkers = int(sys.argv[1]), int(sys.argv[2])

base = np.load('baseline_planetGR.npz', allow_pickle=True)
names = list(base['names'])
fmap = list(base['fmap'])
scales = base['scales']
res0 = base['res']
steps = np.load('steps.npz', allow_pickle=True)
abs_steps = steps['abs_steps']  # aligned with fmap order
nfit = len(fmap)

os.makedirs('jac', exist_ok=True)
todo = [j for j in range(nfit) if j % nworkers == worker
        and not os.path.exists('jac/col_%02d_%s.npz' % (j, names[fmap[j]]))]
if not todo:
    print('WORKER %d: nothing to do' % worker); sys.exit(0)

t0 = time.perf_counter()
fit = PyFittriple('parfile-planetGR-max-bestfit',
                  '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
lnp0 = fit.Compute_lnposterior()
print('WORKER %d ready in %.0fs lnp0=%s todo=%s' % (worker, time.perf_counter()-t0, lnp0, todo), flush=True)

for j in todo:
    pname = names[fmap[j]]
    h_abs = float(abs_steps[j])
    h_rel = h_abs / float(scales[fmap[j]])
    sh = np.zeros(nfit)
    t1 = time.perf_counter()
    sh[j] = +h_rel
    fit.Set_fitted_parameter_relativeshifts(sh)
    lnp_p = fit.Compute_lnposterior()
    res_p = fit.Get_time_residuals().copy()
    sh[j] = -h_rel
    fit.Set_fitted_parameter_relativeshifts(sh)
    lnp_m = fit.Compute_lnposterior()
    res_m = fit.Get_time_residuals().copy()
    dcol = (res_p - res_m) / (2.0 * h_abs)
    dev = float(np.abs(res_p - res0).max())
    np.savez('jac/col_%02d_%s.npz' % (j, pname), dcol=dcol, h_abs=h_abs,
             h_rel=h_rel, lnp_p=lnp_p, lnp_m=lnp_m, maxdev_us=dev)
    print('WORKER %d col %02d %-18s %.0fs maxdev=%.3g us dlnp=(%.3g,%.3g)' %
          (worker, j, pname, time.perf_counter()-t1, dev, lnp_p-lnp0, lnp_m-lnp0), flush=True)
print('WORKER %d DONE' % worker, flush=True)
