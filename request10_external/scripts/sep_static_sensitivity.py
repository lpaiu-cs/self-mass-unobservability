#!/usr/bin/env python3
"""Request 10.8 Gate F1: static SEP lever-arm reproduction (data-independent).

sigma_Delta = marginal 1-sigma of a constant SEP_D amplitude fitted jointly
with the 10.7e headline nuisance (28 params + offset + 30 Fourier pairs,
truncated SVD 1e-3, noise scale s = 1.1148). PASS iff sigma_Delta <= 5e-6
(factor ~3 of the published static sensitivity ~1.8e-6, Voisin+20). An
over-tight sigma (< 0.5e-6) triggers the R5 nuisance-completeness review.

No real-data statistic is computed (design matrix only).
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

base = np.load('baseline_planetGR.npz', allow_pickle=True)
errs = base['errs'].astype(np.float64)
t = base['toas'].astype(np.float64)
N = errs.size
sw = 1.0 / errs
S_SCALE = 1.1148          # 10.7e rn-mode noise scale
SV_CUT = 1e-3
M_RN = 30

J = np.load('finite_jacobian_v2.npy')
d = np.load('sep_dynamic/col_SEP_D.npz')
jD = d['dcol'].astype(np.float64)          # us per unit Delta

A = sw[:, None] * J
An = A / np.linalg.norm(A, axis=0)
cvec = sw / np.linalg.norm(sw)
T_SPAN = float(t.max() - t.min())
fcols = []
for j in range(1, M_RN + 1):
    arg = 2.0*np.pi*j*(t - t.min())/T_SPAN
    fcols += [np.cos(arg), np.sin(arg)]
FB = sw[:, None]*np.column_stack(fcols)
FB = FB/np.linalg.norm(FB, axis=0, keepdims=True)

def sigma_and_survival(B0):
    U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
    k0 = int(np.sum(s0/s0[0] >= SV_CUT))
    Qk = U0[:, :k0]
    x = sw * jD
    x2 = x - Qk @ (Qk.T @ x)
    return (S_SCALE / np.linalg.norm(x2),
            float(np.linalg.norm(x2)/np.linalg.norm(x)), k0)

sig_rn, surv_rn, k_rn = sigma_and_survival(np.column_stack([An, cvec, FB]))
sig_wh, surv_wh, k_wh = sigma_and_survival(np.column_stack([An, cvec]))

out = {
    'column': 'sep_dynamic/col_SEP_D.npz',
    'h_abs': float(d['h_abs']), 'maxdev_us': float(d['maxdev_us']),
    'halfstep_lin_dev': float(d['halfstep_lin_dev']),
    'sigma_Delta_rn_marginalized': sig_rn,
    'sigma_Delta_white': sig_wh,
    'survival_fraction_rn': surv_rn, 'survival_fraction_white': surv_wh,
    'nuisance_rank_rn': k_rn, 'nuisance_rank_white': k_wh,
    'reference_published_sigma': 1.8e-6,
    'pass_F1': bool(sig_rn <= 5e-6),
    'R5_review_needed': bool(sig_rn < 0.5e-6),
}
with open('sep_dynamic/sep_static_sensitivity.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('F1: sigma_Delta (rn) = %.3e  (white %.3e); survival %.3f/%.3f; '
      'published ~1.8e-6 -> %s%s'
      % (sig_rn, sig_wh, surv_rn, surv_wh,
         'PASS' if out['pass_F1'] else 'FAIL',
         ' [R5 review: over-tight]' if out['R5_review_needed'] else ''))
