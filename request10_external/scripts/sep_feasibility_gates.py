#!/usr/bin/env python3
"""Request 10.8 gates F2/F3 (10.8f-corrected re-adjudication).

C1: templates come from sep_common with the CAUSAL sign (the original run
used the advanced-response combination). C4: the static guard is kept
explicitly. C10: the anchor ratio is read from the F1 artifact.
Data-independent (design matrices only).
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

import sep_common as sc

inp = sc.load_inputs(HERE)
nu = sc.build_nuisance(inp)
Q = nu['Q']
S_SCALE, CHI2RED = sc.noise_scale(Q, inp['sw'], inp['res0'])
K = sc.load_anchor_factors(HERE)['K_static']
sw = inp['sw']

taus = np.geomspace(5.0, 500.0, 41)
rows = []
for tau in taus:
    Tc, Tb = sc.templates(inp, tau)
    Xs = np.column_stack([sw*Tc, sw*Tb])
    X2 = sc.proj_out(Q, Xs)
    surv = float(np.linalg.norm(X2[:, 1])/np.linalg.norm(Xs[:, 1]))
    Minv = np.linalg.inv(X2.T @ X2)
    sig = S_SCALE*math.sqrt(Minv[1, 1])
    rows.append({'tau': float(tau), 'survival': surv,
                 'sigma_fisher': sig, 'sigma_anchored': sig*K})

best_surv = max(rows, key=lambda r: r['survival'])
best_sig = min(rows, key=lambda r: r['sigma_anchored'])
per_col = {}
for key, col in inp['cols'].items():
    x = sw*col
    per_col[key] = {'survival': float(np.linalg.norm(sc.proj_out(Q, x))/np.linalg.norm(x))}

F2 = best_surv['survival'] >= 0.03
F3 = best_sig['sigma_anchored'] <= 1e-5
out = {
    'preregistration': '../notes/REQUEST10_8_DYNAMIC_SEP_POLARIZATION_CHANNEL.md',
    'correction': '../notes/REQUEST10_8F_REVIEW_RESPONSE.md (C1 causal sign, C4 guard, C10 K from artifact)',
    'anchor_ratio_K': K, 'noise_scale': S_SCALE, 'chi2_red': CHI2RED,
    'nuisance_rank': nu['rank'], 'guard_residual_norm': nu['guard_residual_norm'],
    'per_column': per_col,
    'curve': rows,
    'F2': {'best_survival': best_surv['survival'], 'at_tau_d': best_surv['tau'],
           'threshold': 0.03, 'pass': bool(F2)},
    'F3': {'min_sigma_anchored': best_sig['sigma_anchored'],
           'min_sigma_fisher': min(r['sigma_fisher'] for r in rows),
           'at_tau_d': best_sig['tau'], 'threshold': 1e-5, 'pass': bool(F3)},
}
with open('sep_dynamic/sep_feasibility_gates.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('per-column survival:', {k: round(v['survival'], 4) for k, v in per_col.items()})
print('F2 (causal): best survival %.4f at tau=%.1f -> %s'
      % (best_surv['survival'], best_surv['tau'], 'PASS' if F2 else 'FAIL'))
print('F3 (causal): min anchored sigma = %.3e at tau=%.1f (Fisher %.3e) -> %s'
      % (best_sig['sigma_anchored'], best_sig['tau'],
         out['F3']['min_sigma_fisher'], 'PASS' if F3 else 'FAIL'))
