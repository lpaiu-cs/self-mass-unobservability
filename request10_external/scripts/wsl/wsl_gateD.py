#!/usr/bin/env python3
"""Request 10.7d gates D1 (live null calibration) and D2 (estimator
calibration at the limit), pre-registered in
notes/REQUEST10_EXTERNAL_JOINT_FIT_PROMOTION_10_7D.md.

GN protocol identical to 10.7c Gate A: truncated-pinv updates over the 28
parameter directions (rel SV cut 1e-3 on the weighted, column-normalized
Jacobian), <= 4 iterations, 0.5% convergence break, non-cumulative shifts.
Measurement: truncated joint pipeline (B0 span incl. offset column) on the
post-GN remainder.
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
cols = np.load('gateD_columns.npz')
gp = json.load(open('gateD_params.json'))
SV_CUT = gp['sv_cut']
ANCH = sorted(gp['anchors'], key=float)

A = w[:, None] * J
cn = np.linalg.norm(A, axis=0)
An = A / cn
cvec = w / np.linalg.norm(w)
B0 = np.column_stack([An, cvec])
U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
k0 = int(np.sum(s0 / s0[0] >= SV_CUT))
Qk = U0[:, :k0]
Ua, sa, Vat = np.linalg.svd(An, full_matrices=False)
ka = int(np.sum(sa / sa[0] >= SV_CUT))

def trunc_update(r_us):
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

def gn_absorb(beta_inj, xb, label):
    """Pre-registered semantics: the model fits the (real + injected) data.
    Minimize ||w * res_syn(p0+dtheta)||, res_syn(p) = res(p) + beta_inj*xb.
    Iteration 0 uses res_syn(p0) = res0 + beta_inj*xb (no eval needed)."""
    dth = np.zeros(28)
    r = res0 + beta_inj * xb
    hist = [float(np.linalg.norm(w * r))]
    for it in range(1, 5):
        dth = dth - trunc_update(r)
        r = eval_at(dth) + beta_inj * xb
        hist.append(float(np.linalg.norm(w * r)))
        print('  %s GN iter %d: |w r| = %.3f' % (label, it, hist[-1]), flush=True)
        if it >= 2 and abs(hist[-2] - hist[-1]) / max(hist[-2], 1e-12) < 0.005:
            break
    return r, hist

def measure(r_us, tau_key):
    xb = cols['xb_%s' % tau_key]
    Xs = np.column_stack([w * cols['xc'], w * xb])
    X2 = Xs - Qk @ (Qk.T @ Xs)
    Minv = np.linalg.inv(X2.T @ X2)
    yr = w * r_us
    th = Minv @ (X2.T @ (yr - Qk @ (Qk.T @ yr)))
    return float(th[1])

out = {'params': gp, 'D1': {}, 'D2': {}}

# ---- Gate D1: live null calibration (one GN run, measured at all anchors) ----
print('GATE D1 (null):', flush=True)
r_null, hist_null = gn_absorb(0.0, np.zeros(N), 'null')
out['D1']['W_norm_history'] = hist_null
d1_all = True
for a in ANCH:
    ap = gp['anchors'][a]
    bhat = measure(r_null, a)
    z_gn = bhat / ap['sigma_beta']
    ok = abs(z_gn - ap['z_lin']) <= gp['d1_tol_z']
    d1_all = d1_all and ok
    out['D1']['tau_%s' % a] = {'beta_hat_GN': bhat, 'z_GN': z_gn,
                               'z_lin': ap['z_lin'], 'pass': bool(ok)}
    print('  D1 tau=%s: z_GN=%.3f z_lin=%.3f -> %s'
          % (a, z_gn, ap['z_lin'], 'PASS' if ok else 'FAIL'), flush=True)
out['D1']['pass'] = bool(d1_all)

# ---- Gate D2: estimator calibration at the limit amplitude ----
print('GATE D2 (injections at u95):', flush=True)
d2_all = True
for a in ANCH:
    ap = gp['anchors'][a]
    b_inj = ap['u95']
    r_fin, hist = gn_absorb(b_inj, cols['xb_%s' % a], 'tau=%s' % a)
    bhat = measure(r_fin, a)
    dev_sigma = (bhat - ap['beta_hat_lin'] - b_inj) / ap['sigma_beta']
    ok = abs(dev_sigma) <= gp['d2_tol_sigma']
    d2_all = d2_all and ok
    out['D2']['tau_%s' % a] = {'beta_inj': b_inj, 'beta_hat_GN': bhat,
                               'beta_hat_lin_null': ap['beta_hat_lin'],
                               'deviation_sigma': dev_sigma,
                               'W_norm_history': hist, 'pass': bool(ok)}
    print('  D2 tau=%s: beta_inj=%.4f beta_hat_GN=%.4f dev=%.3f sigma -> %s'
          % (a, b_inj, bhat, dev_sigma, 'PASS' if ok else 'FAIL'), flush=True)
out['D2']['pass'] = bool(d2_all)

with open('jointfit_gateD.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('GATED_DONE D1=%s D2=%s' % (out['D1']['pass'], out['D2']['pass']), flush=True)
