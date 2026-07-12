#!/usr/bin/env python3
"""Request 10.7c gates L and A (pre-registered in
notes/REQUEST10_EXTERNAL_JOINT_FIT_AMENDMENT_10_7C.md).

Gate L: linearity of the truncated tau=52d displacement.
Gate A: can the REAL model (Gauss-Newton with live Nutimo residuals,
J_v2-truncated updates) absorb an injected chi signal beyond the
truncated-linear prediction?
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
sc = np.array([float(scales[p]) for p in fmap])
w = 1.0 / errs
N = res0.size

J = np.load('finite_jacobian_v2.npy')
dth_L = np.load('jointfit_dtheta52_trunc.npy')
xb52 = np.load('xb52.npy')
xc52 = np.load('xc52.npy')
gp = json.load(open('gateLA_params.json'))
SV_CUT = gp['sv_cut']
BETA_TEST = gp['beta_test_us']
SIGB = gp['sigma_beta_trunc_52']

# truncated spans (identical construction to the Windows pipeline)
A = w[:, None] * J
cn = np.linalg.norm(A, axis=0)
An = A / cn
cvec = w / np.linalg.norm(w)
B0 = np.column_stack([An, cvec])
U0, s0, V0t = np.linalg.svd(B0, full_matrices=False)
k0 = int(np.sum(s0 / s0[0] >= SV_CUT))
Qk = U0[:, :k0]
Ua, sa, Vat = np.linalg.svd(An, full_matrices=False)
ka = int(np.sum(sa / sa[0] >= SV_CUT))

def trunc_update(r_us):
    """truncated-pinv GN update (parameter units) from a us-residual vector."""
    c = Vat[:ka].T @ ((Ua[:, :ka].T @ (w * r_us)) / sa[:ka])
    return c / cn

t0 = time.perf_counter()
fit = PyFittriple('parfile-planetGR-max-bestfit',
                  '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
lnp0 = fit.Compute_lnposterior()
print('constructed in %.0fs lnp0=%.6f' % (time.perf_counter() - t0, lnp0), flush=True)

def eval_at(dtheta):
    fit.Set_fitted_parameter_relativeshifts(dtheta / sc)
    fit.Compute_lnposterior()
    return fit.Get_time_residuals().copy()

out = {'params': gp}

# ---- Gate L ----
t1 = time.perf_counter()
res_L = eval_at(dth_L)
d_act = res_L - res0
d_lin = J @ dth_L
num = float(np.linalg.norm(w * (d_act - d_lin)))
den = float(np.linalg.norm(w * d_lin))
out['gate_L'] = {'relative_deviation': num / den, 'W_norm_linear': den,
                 'pass': bool(num / den < 0.05), 'seconds': time.perf_counter() - t1}
print('GATE L: dev=%.4f (|lin|_W=%.2f) -> %s'
      % (num / den, den, 'PASS' if out['gate_L']['pass'] else 'FAIL'), flush=True)

# ---- Gate A ----
d_syn = res0 + BETA_TEST * xb52
dth = np.zeros(28)
r = d_syn - res0                      # iteration 0 residual (no eval needed)
hist = [float(np.linalg.norm(w * r))]
for it in range(1, 5):
    dth = dth + trunc_update(r)
    t1 = time.perf_counter()
    res_it = eval_at(dth)
    r = d_syn - res_it
    hist.append(float(np.linalg.norm(w * r)))
    print('GN iter %d: |w r| = %.3f (%.0fs)' % (it, hist[-1], time.perf_counter() - t1), flush=True)
    if it >= 2 and (hist[-2] - hist[-1]) / hist[-2] < 0.005:
        break

# measure surviving signal through the truncated joint pipeline on the remainder
Xs = np.column_stack([w * xc52, w * xb52])
X2 = Xs - Qk @ (Qk.T @ Xs)
Minv = np.linalg.inv(X2.T @ X2)
yr = w * r
yr_perp = yr - Qk @ (Qk.T @ yr)
th = Minv @ (X2.T @ yr_perp)
z_true = float(th[1] / SIGB)
out['gate_A'] = {
    'beta_test_us': BETA_TEST, 'gn_iterations': len(hist) - 1,
    'W_norm_history': hist, 'beta_hat_surviving': float(th[1]),
    'z_true': z_true, 'pred_z_trunc': gp['pred_z_trunc'],
    'pred_z_fullrank': gp['pred_z_fullrank'],
    'pass': bool(z_true >= gp['z_pass_threshold']),
    'max_dtheta_over_step_final': None,
}
print('GATE A: z_true=%.2f (trunc pred %.1f, full-rank pred %.2f, threshold %.1f) -> %s'
      % (z_true, gp['pred_z_trunc'], gp['pred_z_fullrank'], gp['z_pass_threshold'],
         'PASS' if out['gate_A']['pass'] else 'FAIL'), flush=True)

with open('jointfit_gateLA.json', 'w') as fh:
    json.dump(out, fh, indent=1)
np.save('gateA_final_dtheta.npy', dth)
print('GATES_DONE L=%s A=%s' % (out['gate_L']['pass'], out['gate_A']['pass']), flush=True)
