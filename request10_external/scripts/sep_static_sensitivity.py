#!/usr/bin/env python3
"""Request 10.8 Gate F1 (10.8f-corrected): static SEP lever-arm audit.

10.8f C10/C12 fixes vs the original: noise scales are COMPUTED per span
(the rn-mode 1.1148 was previously applied to the white-mode output too),
and the nuisance construction comes from sep_common. Data-independent.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

import sep_common as sc

inp = sc.load_inputs(HERE)
sw, res0, jD = inp['sw'], inp['res0'], inp['jD']

# For the STATIC column the guard direction IS the signal: build the spans
# without the jD guard column (keep_guard=False, and exclude jD from B0 by
# projecting the signal only off Jacobian+offset+Fourier).
def span_sigma(m_rn):
    A = sw[:, None]*inp['J']
    An = A/np.linalg.norm(A, axis=0)
    cvec = sw/np.linalg.norm(sw)
    blocks = [An, cvec[:, None]]
    if m_rn > 0:
        t = inp['t']
        T = float(t.max() - t.min())
        fc = []
        for j in range(1, m_rn + 1):
            arg = 2.0*np.pi*j*(t - t.min())/T
            fc += [np.cos(arg), np.sin(arg)]
        FB = sw[:, None]*np.column_stack(fc)
        blocks.append(FB/np.linalg.norm(FB, axis=0, keepdims=True))
    B0 = np.column_stack(blocks)
    U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
    k = int(np.sum(s0/s0[0] >= sc.SV_CUT_DEFAULT))
    Q = U0[:, :k]
    s_scale, chi2red = sc.noise_scale(Q, sw, res0)
    x = sw*jD
    x2 = sc.proj_out(Q, x)
    return (s_scale/np.linalg.norm(x2),
            float(np.linalg.norm(x2)/np.linalg.norm(x)), k, s_scale, chi2red)

sig_rn, surv_rn, k_rn, s_rn, c_rn = span_sigma(sc.M_RN_DEFAULT)
sig_wh, surv_wh, k_wh, s_wh, c_wh = span_sigma(0)

# R5 extension columns (live absorbers only)
import glob
r5_kept, r5_null = [], []
cols5 = []
for f in sorted(glob.glob('sep_dynamic/col_R5_*.npz')):
    name = os.path.basename(f)[7:-4]
    d5 = np.load(f)
    if float(d5['maxdev_us']) < 1e-4:
        r5_null.append(name)
        continue
    c = sw*d5['dcol'].astype(np.float64)
    cols5.append(c/np.linalg.norm(c))
    r5_kept.append(name)

sig_x = sig_rn
if cols5:
    A = sw[:, None]*inp['J']
    An = A/np.linalg.norm(A, axis=0)
    cvec = sw/np.linalg.norm(sw)
    t = inp['t']
    T = float(t.max() - t.min())
    fc = []
    for j in range(1, sc.M_RN_DEFAULT + 1):
        arg = 2.0*np.pi*j*(t - t.min())/T
        fc += [np.cos(arg), np.sin(arg)]
    FB = sw[:, None]*np.column_stack(fc)
    FB = FB/np.linalg.norm(FB, axis=0, keepdims=True)
    B0 = np.column_stack([An, cvec, FB] + [c[:, None] for c in cols5])
    U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
    Q = U0[:, :int(np.sum(s0/s0[0] >= sc.SV_CUT_DEFAULT))]
    s_x, _ = sc.noise_scale(Q, sw, res0)
    x2 = sc.proj_out(Q, sw*jD)
    sig_x = s_x/np.linalg.norm(x2)

PUBLISHED = 1.8e-6
out = {
    'column': 'sep_dynamic/col_SEP_D.npz',
    'sigma_Delta_rn_marginalized': sig_rn,
    'sigma_Delta_white': sig_wh,
    'noise_scale_rn': s_rn, 'noise_scale_white': s_wh,
    'chi2_red_rn': c_rn, 'chi2_red_white': c_wh,
    'survival_fraction_rn': surv_rn, 'survival_fraction_white': surv_wh,
    'nuisance_rank_rn': k_rn, 'nuisance_rank_white': k_wh,
    'R5_extended': {'columns_kept': r5_kept, 'columns_null': r5_null,
                    'sigma_Delta_rn_R5': sig_x},
    'reference_published_sigma': PUBLISHED,
    'K_static_ratio': PUBLISHED/sig_x,
    'pass_F1_validation': bool(5e-7 <= sig_x <= 5e-6),
    'decision': ('lever-arm validated; use Fisher sigma'
                 if 5e-7 <= sig_x <= 5e-6 else
                 'Fisher sigma untrustworthy (%.3e, %.0fx vs published); '
                 'DOWNGRADE: anchor to published sigma_Delta = %.1e'
                 % (sig_x, PUBLISHED/sig_x, PUBLISHED)),
    'static_anchor_used': (PUBLISHED if not (5e-7 <= sig_x <= 5e-6) else sig_x),
}
with open('sep_dynamic/sep_static_sensitivity.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('F1: sigma_Delta rn=%.4e (s=%.4f) | white=%.4e (s=%.4f) | R5-ext=%.4e'
      % (sig_rn, s_rn, sig_wh, s_wh, sig_x))
print('K_static = %.2f | DECISION: %s' % (out['K_static_ratio'], out['decision']))
