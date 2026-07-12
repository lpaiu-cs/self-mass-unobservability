#!/usr/bin/env python3
"""Request 10.8b real-data fit (10.8f-corrected re-run).

Corrections vs the original run (notes/REQUEST10_8F_REVIEW_RESPONSE.md):
C1 causal template sign (the original scanned the advanced-response branch
as 'causal' and vice versa); C4 active static-leakage guard; C10 noise scale
computed and K factors read from artifacts. Deterministic (seed 20260710).
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

NSIM = 500
STOP_S = 1.75
ANCHORS = [2.0, 5.0, 18.0, 52.0, 200.0]

inp = sc.load_inputs(HERE)
sw, res0, N = inp['sw'], inp['res0'], inp['N']
nu = sc.build_nuisance(inp)
Q0 = nu['Q']
S_SCALE, CHI2RED = sc.noise_scale(Q0, sw, res0)
KF = sc.load_anchor_factors(HERE)
print('nuisance rank %d/%d (guard residual %.3g kept); chi2_red = %.4f; s = %.4f; K_dyn=%g K_static=%.2f'
      % (nu['rank'], nu['ncols'], nu['guard_residual_norm'], CHI2RED, S_SCALE,
         KF['K_dyn'], KF['K_static']))
if S_SCALE > STOP_S:
    raise SystemExit('STOP RULE: s > %.2f' % STOP_S)

y_perp = sc.proj_out(Q0, sw*res0)

def scan(yv_perp, grid, anticausal=False, keep=False):
    rows = []
    for tau in grid:
        Tc, Tb = sc.templates(inp, tau, anticausal=anticausal)
        Xs = np.column_stack([sw*Tc, sw*Tb])
        X2 = sc.proj_out(Q0, Xs)
        M = X2.T @ X2
        th = np.linalg.solve(M, X2.T @ yv_perp)
        Minv = np.linalg.inv(M)
        sig = S_SCALE*math.sqrt(Minv[1, 1])
        r = {'tau': float(tau), 'c_hat': float(th[0]), 'beta_hat': float(th[1]),
             'sigma_fisher': sig, 'z': float(th[1]/sig),
             'u95_fisher': sc.u95_of(float(th[1]), sig),
             'u95_Kdyn': sc.u95_of(float(th[1]), KF['K_dyn']*sig),
             'u95_Kstatic': sc.u95_of(float(th[1]), KF['K_static']*sig),
             'survival': float(np.linalg.norm(X2[:, 1])/np.linalg.norm(Xs[:, 1]))}
        if keep:
            r['_X2'] = X2
            r['_Minv'] = Minv
        rows.append(r)
    return rows

grid = sorted(set(np.geomspace(2.0, 500.0, 61).tolist() + ANCHORS))
rows = scan(y_perp, grid, keep=True)
Z = max(abs(r['z']) for r in rows)
best = min(rows, key=lambda r: r['u95_Kdyn'])
print('real-data (causal): Z = %.3f; min u95_Kdyn = %.4g at tau = %.1f (Fisher %.3g)'
      % (Z, best['u95_Kdyn'], best['tau'], best['u95_fisher']))

rng = np.random.default_rng(sc.SEED)
P = np.vstack([(r['_Minv'] @ r['_X2'].T)[1] for r in rows])
sv = np.array([r['sigma_fisher'] for r in rows])
Zs = np.max(np.abs(P @ (rng.standard_normal((N, NSIM))*S_SCALE))/sv[:, None], axis=0)
p_global = float(np.mean(Zs >= Z))
print('null: median %.3f -> p_global = %.4f' % (float(np.median(Zs)), p_global))

rows_ac = scan(y_perp, grid, anticausal=True, keep=True)
Zac = max(abs(r['z']) for r in rows_ac)
Pa = np.vstack([(r['_Minv'] @ r['_X2'].T)[1] for r in rows_ac])
sva = np.array([r['sigma_fisher'] for r in rows_ac])
Zsa = np.max(np.abs(Pa @ (rng.standard_normal((N, NSIM))*S_SCALE))/sva[:, None], axis=0)
p_ac = float(np.mean(Zsa >= Zac))
print('anti-causal control: Z = %.3f, p = %.4f' % (Zac, p_ac))

detection = bool(p_global < 0.003 and p_ac > 0.05)

injections = {}
for a in ANCHORS:
    r0 = sc.row_at(rows, a)
    b_inj = r0['u95_Kdyn']
    _, Tb = sc.templates(inp, a)
    y_inj = sc.proj_out(Q0, sw*(res0 + b_inj*Tb))
    ri = scan(y_inj, [a])[0]
    dev = (ri['beta_hat'] - r0['beta_hat'] - b_inj)/r0['sigma_fisher']
    injections[sc.anchor_key(a)] = {'beta_inj': b_inj, 'beta_hat': ri['beta_hat'],
                                    'deviation_sigma_fisher': dev,
                                    'pass': bool(abs(dev) <= 0.5)}
    print('inject tau=%g: dev = %.3f sigma_F -> %s'
          % (a, dev, 'PASS' if abs(dev) <= 0.5 else 'FAIL'))

out = {
    'preregistration': '../notes/REQUEST10_8B_REALDATA_DYNAMIC_SEP_FIT.md',
    'correction': '../notes/REQUEST10_8F_REVIEW_RESPONSE.md',
    'seed': sc.SEED, 'nsim': NSIM,
    'K_dyn': KF['K_dyn'], 'K_static': KF['K_static'],
    'noise': {'chi2_red': CHI2RED, 's': S_SCALE, 'rank': nu['rank'],
              'guard_residual_norm': nu['guard_residual_norm']},
    'real_data': {'Z': Z, 'p_global': p_global, 'detection_candidate': detection,
                  'min_u95_Kdyn': best['u95_Kdyn'], 'tau_at_min': best['tau'],
                  'min_u95_fisher': min(r['u95_fisher'] for r in rows)},
    'anticausal_control': {'Z': Zac, 'p': p_ac, 'quiet': bool(p_ac > 0.05)},
    'injections': injections,
    'anchors': {sc.anchor_key(a): {k: v for k, v in sc.row_at(rows, a).items()
                                   if not k.startswith('_')}
                for a in ANCHORS},
    'quoted_window_d': [2.0, 500.0],
    'convention': ('CAUSAL template of record: response to the lagging pole; '
                   'col_s is the +pi/2 (-sin) drive response, causal = g1*col_c - g2*col_s. '
                   'beta_ff = dimensionless Delta-oscillation pole amplitude at unit drive.'),
}
with open('sep_dynamic/sep_joint_fit_10_8b.json', 'w') as fh:
    json.dump(out, fh, indent=1)

hdr = 'tau_d\tbeta_hat\tsigma_fisher\tu95_fisher\tu95_Kdyn\tu95_Kstatic\tz\tsurvival'
lines = ['%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.4f\t%.4f' %
         (r['tau'], r['beta_hat'], r['sigma_fisher'], r['u95_fisher'],
          r['u95_Kdyn'], r['u95_Kstatic'], r['z'], r['survival']) for r in rows]
with open('sep_dynamic/sep_beta_limit_curve_10_8b.tsv', 'w') as fh:
    fh.write(hdr + '\n' + '\n'.join(lines) + '\n')

# gate-G inputs (causal templates of record)
np.savez('sep_dynamic/gateG_inputs.npz',
         T18=sc.templates(inp, 18.0)[1], T52=sc.templates(inp, 52.0)[1],
         Tc=sc.templates(inp, 18.0)[0])
gg = {'anchors': {}, 'sv_cut': sc.SV_CUT_DEFAULT, 'convention': 'causal (10.8f C1)'}
for a in (18.0, 52.0):
    r0 = sc.row_at(rows, a)
    gg['anchors'][sc.anchor_key(a)] = {'u95_Kdyn': r0['u95_Kdyn'],
                                       'u95_Kstatic': r0['u95_Kstatic'],
                                       'sigma_fisher': r0['sigma_fisher'],
                                       'beta_hat_lin': r0['beta_hat'],
                                       'z_lin': r0['z']}
json.dump(gg, open('sep_dynamic/gateG_params.json', 'w'), indent=1)
print('10_8B_FIT_DONE detection=%s min_u95_Kdyn=%.4g @ tau=%.1f'
      % (detection, best['u95_Kdyn'], best['tau']))
