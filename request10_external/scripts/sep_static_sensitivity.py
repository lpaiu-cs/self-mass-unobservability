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

# R5 extension columns (fixed physically-degenerate parameters), if derived.
# Guard: several candidates turn out NULL in this parfile's active
# parametrization (dcol ~ 0 -> not an absorber here); keep only live columns.
import glob
r5_files = sorted(glob.glob('sep_dynamic/col_R5_*.npz'))
R5 = None
r5_kept, r5_null = [], []
if r5_files:
    cols = []
    for f in r5_files:
        name = os.path.basename(f)[7:-4]
        d5 = np.load(f)
        # A column is a usable nuisance DIRECTION iff the parameter actually
        # moves residuals (absolute floor, not relative to the huge SEP column).
        if float(d5['maxdev_us']) < 1e-4:
            r5_null.append(name)
            continue
        c = sw * d5['dcol'].astype(np.float64)
        cols.append(c / np.linalg.norm(c))
        r5_kept.append(name)
    if cols:
        R5 = np.column_stack(cols)
    print('R5 columns kept: %s | NULL (inert in this parametrization): %s'
          % (r5_kept or '(none)', r5_null or '(none)'))

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
}
if R5 is not None:
    sig_x, surv_x, k_x = sigma_and_survival(np.column_stack([An, cvec, FB, R5]))
    out['R5_extended'] = {'columns_kept': r5_kept, 'columns_null': r5_null,
                          'sigma_Delta_rn_R5': sig_x,
                          'survival_fraction_rn_R5': surv_x,
                          'nuisance_rank_rn_R5': k_x}
    sig_gate = sig_x
else:
    out['R5_extended'] = {'columns_kept': [], 'columns_null': r5_null,
                          'note': 'all candidate absorber columns null in this parametrization'}
    sig_gate = sig_rn
# Pre-registered decision: F1 validates the lever arm only if our Fisher sigma
# lands within a factor ~3 of published (<=5e-6) AND is not absurdly tight.
# Here it is ~1e3x too tight and R5 could not close the gap -> DOWNGRADE:
# anchor the channel's static sensitivity to the published marginal.
PUBLISHED = 1.8e-6
out['pass_F1_validation'] = bool(5e-7 <= sig_gate <= 5e-6)
out['decision'] = ('lever-arm validated; use Fisher sigma' if out['pass_F1_validation']
                   else 'Fisher sigma untrustworthy (%.2e, %.0fx vs published); '
                        'DOWNGRADE: anchor to published sigma_Delta = %.1e'
                        % (sig_gate, PUBLISHED / sig_gate, PUBLISHED))
out['static_anchor_used'] = PUBLISHED if not out['pass_F1_validation'] else sig_gate
with open('sep_dynamic/sep_static_sensitivity.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('F1: sigma_Delta (rn) = %.3e  (white %.3e); survival %.4f/%.4f'
      % (sig_rn, sig_wh, surv_rn, surv_wh))
if R5 is not None:
    print('F1+R5 (kept %s): sigma_Delta = %.3e; survival %.4f; rank %d'
          % (r5_kept, sig_x, surv_x, k_x))
print('DECISION: %s' % out['decision'])
