#!/usr/bin/env python3
"""Request 10.8 gates F2 (dynamic non-degeneracy) and F3 (projected
sensitivity) — data-independent, on the T2 exact response columns.

Template model (linear in Delta): a pole-lagged drive
Delta(t) = A * Re[e^{i w t}/(1 + i w tau)] = A [ g1 cos(w t) + g2 sin(w t) ],
g1 = 1/(1+w^2 tau^2), g2 = w tau/(1+w^2 tau^2), maps by linearity of the
orbit response onto  T_w(tau) = g1 * col_w_c + g2 * col_w_s.
Shared-tau template over the three drive carriers (unit drive, integrator
phase origin):  T_beta(tau) = sum_w T_w(tau);  T_cY = sum_w col_w_c.

F2: survival fraction of T_beta(tau) against the 10.7e headline nuisance
    (28 params + offset + 30 Fourier pairs, SV cut 1e-3) >= 3% somewhere in
    tau in [5, 500] d.
F3: sigma_beta_ff(tau) from the joint [T_cY, T_beta] fit; quoted BOTH as the
    raw Fisher value AND anchored by the pre-registered downgrade ratio
    K = published/Fisher static = 1.8e-6 / 1.922e-9 (the same
    nuisance-completeness deficiency is assumed to apply; conservative for
    oscillating templates). Gate: min anchored sigma <= 1e-5.
"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

base = np.load('baseline_planetGR.npz', allow_pickle=True)
errs = base['errs'].astype(np.float64)
t = base['toas'].astype(np.float64)
N = errs.size
sw = 1.0/errs
S_SCALE = 1.1148
SV_CUT = 1e-3
M_RN = 30
K_ANCHOR = 1.8e-6 / 1.922e-9        # published static marginal / our Fisher static

OMS = {'in': 3.856136695791377, 'out': 0.019199654300698213, 'dif': 3.836937041490679}
cols = np.load('sep_dynamic/sep_dynamic_columns.npz')

J = np.load('finite_jacobian_v2.npy')
A = sw[:, None]*J
An = A/np.linalg.norm(A, axis=0)
cvec = sw/np.linalg.norm(sw)
T_SPAN = float(t.max() - t.min())
fc = []
for j in range(1, M_RN+1):
    arg = 2.0*np.pi*j*(t - t.min())/T_SPAN
    fc += [np.cos(arg), np.sin(arg)]
FB = sw[:, None]*np.column_stack(fc)
FB = FB/np.linalg.norm(FB, axis=0, keepdims=True)
B0 = np.column_stack([An, cvec, FB])
U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
k0 = int(np.sum(s0/s0[0] >= SV_CUT))
Qk = U0[:, :k0]

def proj_out(x):
    return x - Qk @ (Qk.T @ x)

def templates(tau):
    Tb = np.zeros(N); Tc = np.zeros(N)
    for name, w in OMS.items():
        g1 = 1.0/(1.0 + (w*tau)**2)
        g2 = w*tau/(1.0 + (w*tau)**2)
        Tb += g1*cols['col_%s_c' % name] + g2*cols['col_%s_s' % name]
        Tc += cols['col_%s_c' % name]
    return Tc, Tb

taus = np.geomspace(5.0, 500.0, 41)
rows = []
for tau in taus:
    Tc, Tb = templates(tau)
    Xs = np.column_stack([sw*Tc, sw*Tb])
    X2 = proj_out(Xs)
    surv = float(np.linalg.norm(X2[:, 1])/np.linalg.norm(Xs[:, 1]))
    Minv = np.linalg.inv(X2.T @ X2)
    sig = S_SCALE*math.sqrt(Minv[1, 1])
    rows.append({'tau': float(tau), 'survival': surv,
                 'sigma_fisher': sig, 'sigma_anchored': sig*K_ANCHOR})

best_surv = max(rows, key=lambda r: r['survival'])
best_sig = min(rows, key=lambda r: r['sigma_anchored'])
per_col = {}
for key in cols.files:
    x = sw*cols[key]
    per_col[key] = {'wnorm_us': float(np.linalg.norm(x)),
                    'survival': float(np.linalg.norm(proj_out(x))/np.linalg.norm(x))}

F2 = best_surv['survival'] >= 0.03
F3 = best_sig['sigma_anchored'] <= 1e-5
out = {
    'preregistration': '../notes/REQUEST10_8_DYNAMIC_SEP_POLARIZATION_CHANNEL.md',
    'anchor_ratio_K': K_ANCHOR,
    'per_column': per_col,
    'curve': rows,
    'F2': {'best_survival': best_surv['survival'], 'at_tau_d': best_surv['tau'],
           'threshold': 0.03, 'pass': bool(F2)},
    'F3': {'min_sigma_anchored': best_sig['sigma_anchored'],
           'min_sigma_fisher': min(r['sigma_fisher'] for r in rows),
           'at_tau_d': best_sig['tau'], 'threshold': 1e-5, 'pass': bool(F3)},
    'clock_channel_comparison': {'clock_beta_phys_95': 0.4,
        'projected_gain_anchored': 0.4/(1.96*best_sig['sigma_anchored'])},
}
with open('sep_dynamic/sep_feasibility_gates.json', 'w') as fh:
    json.dump(out, fh, indent=1)

print('per-column survival:', {k: round(v['survival'], 4) for k, v in per_col.items()})
print('F2: best survival %.4f at tau=%.1f d (>=0.03) -> %s'
      % (best_surv['survival'], best_surv['tau'], 'PASS' if F2 else 'FAIL'))
print('F3: min anchored sigma_beta_ff = %.3e at tau=%.1f d (Fisher floor %.3e) (<=1e-5) -> %s'
      % (best_sig['sigma_anchored'], best_sig['tau'],
         out['F3']['min_sigma_fisher'], 'PASS' if F3 else 'FAIL'))
print('projected gain vs clock channel (anchored): %.3g x'
      % out['clock_channel_comparison']['projected_gain_anchored'])
