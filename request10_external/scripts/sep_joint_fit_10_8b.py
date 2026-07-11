#!/usr/bin/env python3
"""Request 10.8b: pre-registered real-data fit of the dynamic free-fall SEP
channel (notes/REQUEST10_8B_REALDATA_DYNAMIC_SEP_FIT.md).

Deterministic (seed 20260710), artifacts-only. FWL joint fit; nuisance =
28 Jacobian(v2) + offset + 30 Fourier pairs + static SEP column; signal =
[T_cY, T_beta(tau)] from the T2 measured response columns. Headline quote =
anchored (K=934) u95; Fisher floor reported alongside. Anti-causal control.
"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

SEED = 20260710
NSIM = 500
SV_CUT = 1e-3
M_RN = 30
STOP_S = 1.75
K_ANCHOR = 934.0
ANCHORS = [2.0, 5.0, 18.0, 52.0, 200.0]
OMS = {'in': 3.856136695791377, 'out': 0.019199654300698213, 'dif': 3.836937041490679}

base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)
errs = base['errs'].astype(np.float64)
t = base['toas'].astype(np.float64)
N = res0.size
sw = 1.0/errs

J = np.load('finite_jacobian_v2.npy')
jD = np.load('sep_dynamic/col_SEP_D.npz')['dcol'].astype(np.float64)
cols = np.load('sep_dynamic/sep_dynamic_columns.npz')

# ---- nuisance block ----
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
xjD = sw*jD
B0 = np.column_stack([An, cvec, FB, xjD/np.linalg.norm(xjD)])
U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
rel0 = s0/s0[0]
rank0 = int(np.sum(rel0 >= SV_CUT))
Q0 = U0[:, :rank0]

def proj_out(X):
    return X - Q0 @ (Q0.T @ X)

y = sw*res0
y_perp = proj_out(y)
chi2_nuis = float(y_perp @ y_perp)
dof = N - rank0
s_scale = math.sqrt(chi2_nuis/dof)
print('nuisance rank %d/%d; chi2_red = %.4f; s = %.4f' %
      (rank0, B0.shape[1], chi2_nuis/dof, s_scale))
if s_scale > STOP_S:
    raise SystemExit('STOP RULE: s > %.2f' % STOP_S)

# ---- templates ----
def templates(tau, anticausal=False):
    sgn = -1.0 if anticausal else 1.0
    Tb = np.zeros(N); Tc = np.zeros(N)
    for name, w in OMS.items():
        g1 = 1.0/(1.0 + (w*tau)**2)
        g2 = sgn*w*tau/(1.0 + (w*tau)**2)
        Tb += g1*cols['col_%s_c' % name] + g2*cols['col_%s_s' % name]
        Tc += cols['col_%s_c' % name]
    return Tc, Tb

def phi(x):
    return 0.5*(1.0 + math.erf(x/math.sqrt(2.0)))

def u95_of(b, sig):
    lo, hi = 0.0, abs(b) + 10.0*sig
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if phi((mid-b)/sig) - phi((-mid-b)/sig) < 0.95:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)

def scan(yv_perp, grid, anticausal=False, keep=False):
    rows = []
    for tau in grid:
        Tc, Tb = templates(tau, anticausal)
        Xs = np.column_stack([sw*Tc, sw*Tb])
        X2 = proj_out(Xs)
        M = X2.T @ X2
        th = np.linalg.solve(M, X2.T @ yv_perp)
        Minv = np.linalg.inv(M)
        sig = s_scale*math.sqrt(Minv[1, 1])
        r = {'tau': float(tau), 'c_hat': float(th[0]), 'beta_hat': float(th[1]),
             'sigma_fisher': sig, 'sigma_anchored': sig*K_ANCHOR,
             'z': float(th[1]/sig),
             'u95_fisher': u95_of(float(th[1]), sig),
             'u95_anchored': u95_of(float(th[1]), sig*K_ANCHOR),
             'survival': float(np.linalg.norm(X2[:, 1])/np.linalg.norm(Xs[:, 1]))}
        if keep:
            r['_X2'] = X2; r['_Minv'] = Minv
        rows.append(r)
    return rows

grid = sorted(set(np.geomspace(2.0, 500.0, 61).tolist() + ANCHORS))
rows = scan(y_perp, grid, keep=True)
Z = max(abs(r['z']) for r in rows)
best = min(rows, key=lambda r: r['u95_anchored'])
print('real-data: Z = %.3f; min u95_anchored = %.4g at tau = %.1f d (Fisher floor %.3g)'
      % (Z, best['u95_anchored'], best['tau'], best['u95_fisher']))

# ---- null calibration ----
rng = np.random.default_rng(SEED)
P = np.vstack([(r['_Minv'] @ r['_X2'].T)[1] for r in rows])
sv = np.array([r['sigma_fisher'] for r in rows])
Nn = rng.standard_normal((N, NSIM))*s_scale
Zs = np.max(np.abs(P @ Nn)/sv[:, None], axis=0)
p_global = float(np.mean(Zs >= Z))
print('null: median %.3f, 99.7%% %.3f -> p_global = %.4f'
      % (float(np.median(Zs)), float(np.quantile(Zs, 0.997)), p_global))

# ---- anti-causal control ----
rows_ac = scan(y_perp, grid, anticausal=True, keep=True)
Zac = max(abs(r['z']) for r in rows_ac)
Pa = np.vstack([(r['_Minv'] @ r['_X2'].T)[1] for r in rows_ac])
sva = np.array([r['sigma_fisher'] for r in rows_ac])
Na = rng.standard_normal((N, NSIM))*s_scale
Zsa = np.max(np.abs(Pa @ Na)/sva[:, None], axis=0)
p_ac = float(np.mean(Zsa >= Zac))
print('anti-causal control: Z = %.3f, p = %.4f' % (Zac, p_ac))

detection = bool(p_global < 0.003 and p_ac > 0.05)

# ---- linear injection-recovery at anchors (anchored u95 amplitudes) ----
injections = {}
for a in ANCHORS:
    r0 = next(r for r in rows if abs(r['tau']-a) < 1e-9)
    b_inj = r0['u95_anchored']
    _, Tb = templates(a)
    y_inj = proj_out(sw*(res0 + b_inj*Tb))
    ri = scan(y_inj, [a])[0]
    dev = (ri['beta_hat'] - r0['beta_hat'] - b_inj)/r0['sigma_fisher']
    injections['tau_%g' % a] = {'beta_inj': b_inj, 'beta_hat': ri['beta_hat'],
                                'deviation_sigma_fisher': dev,
                                'pass': bool(abs(dev) <= 0.5)}
    print('inject tau=%g: beta_inj=%.4g -> dev = %.3f sigma_F -> %s'
          % (a, b_inj, dev, 'PASS' if abs(dev) <= 0.5 else 'FAIL'))

def strip(rr):
    return [{k: v for k, v in r.items() if not k.startswith('_')} for r in rr]

out = {
    'preregistration': '../notes/REQUEST10_8B_REALDATA_DYNAMIC_SEP_FIT.md',
    'seed': SEED, 'nsim': NSIM, 'K_anchor': K_ANCHOR,
    'noise': {'chi2_red': chi2_nuis/dof, 's': s_scale, 'rank': rank0},
    'real_data': {'Z': Z, 'p_global': p_global, 'detection_candidate': detection,
                  'min_u95_anchored': best['u95_anchored'], 'tau_at_min': best['tau'],
                  'min_u95_fisher': min(r['u95_fisher'] for r in rows)},
    'anticausal_control': {'Z': Zac, 'p': p_ac, 'quiet': bool(p_ac > 0.05)},
    'injections': injections,
    'anchors': {('tau_%g' % a): next({'beta_hat': r['beta_hat'],
                                      'sigma_fisher': r['sigma_fisher'],
                                      'u95_anchored': r['u95_anchored'],
                                      'u95_fisher': r['u95_fisher'],
                                      'z': r['z'], 'survival': r['survival']}
                                     for r in rows if abs(r['tau']-a) < 1e-9)
                for a in ANCHORS},
    'quoted_window_d': [2.0, 500.0],
    'convention': ('beta_ff = dimensionless amplitude of the Delta-oscillation pole '
                   'residue at unit drive per carrier, zero drive phase; templates from '
                   'T2 measured integrator response columns'),
}
with open('sep_dynamic/sep_joint_fit_10_8b.json', 'w') as fh:
    json.dump(out, fh, indent=1)

hdr = 'tau_d\tbeta_hat\tsigma_fisher\tu95_fisher\tu95_anchored\tz\tsurvival'
lines = ['%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.4f\t%.4f' %
         (r['tau'], r['beta_hat'], r['sigma_fisher'], r['u95_fisher'],
          r['u95_anchored'], r['z'], r['survival']) for r in strip(rows)]
with open('sep_dynamic/sep_beta_limit_curve_10_8b.tsv', 'w') as fh:
    fh.write(hdr + '\n' + '\n'.join(lines) + '\n')

# save the tau=18 and 52 linear info + templates for the WSL G-gates
np.savez('sep_dynamic/gateG_inputs.npz',
         T18=templates(18.0)[1], T52=templates(52.0)[1],
         Tc=templates(18.0)[0])
gg = {'anchors': {}}
for a in (18.0, 52.0):
    r0 = next(r for r in rows if abs(r['tau']-a) < 1e-9)
    gg['anchors']['%g' % a] = {'u95_anchored': r0['u95_anchored'],
                               'sigma_fisher': r0['sigma_fisher'],
                               'beta_hat_lin': r0['beta_hat'], 'z_lin': r0['z']}
gg['sv_cut'] = SV_CUT
json.dump(gg, open('sep_dynamic/gateG_params.json', 'w'), indent=1)

print('10_8B_FIT_DONE detection=%s min_u95_anchored=%.4g @ tau=%.1f d'
      % (detection, best['u95_anchored'], best['tau']))
