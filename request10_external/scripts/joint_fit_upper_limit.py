#!/usr/bin/env python3
"""Request 10.7b: pre-registered joint-fit shared-tau_chi upper limit (real data).

Implements notes/REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md exactly, on the
returned 10.7a artifacts (no live Nutimo needed):
  baseline_planetGR.npz, finite_jacobian.npy, carrier_projection_rank.json

Joint fitting per C2 via Frisch-Waugh-Lovell block elimination: both the data
and the signal columns are projected off the 29-column nuisance span before the
2x2 (or 4x4 comparator) solve -- algebraically identical to the simultaneous
weighted least squares of all 31 (33) amplitudes, and NOT the sequential
project-then-fit that biased tau-hat in 10.7a Stage 3.

Stdlib + numpy only. Deterministic (fixed RNG seed 20260710).
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # request10_external/
os.chdir(HERE)

# optional: argv[1] = jacobian file (default v1), argv[2] = output suffix,
#           argv[3] = relative-SV cut for the nuisance span (default 1e-10 = full rank),
#           argv[4] = quoted-window minimum tau in days (default 10.0; 10.7d uses 1.0),
#           argv[5] = number of red-noise Fourier pairs appended to the nuisance (10.7e E3),
#           argv[6] = signal template: 'unit' (default) or 'phys' (10.7e E2 dictionary drive)
JAC_PATH = sys.argv[1] if len(sys.argv) > 1 else 'finite_jacobian.npy'
SUFFIX = sys.argv[2] if len(sys.argv) > 2 else ''
SV_CUT = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
WIN_MIN = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
M_RN = int(sys.argv[5]) if len(sys.argv) > 5 else 0
TEMPLATE = sys.argv[6] if len(sys.argv) > 6 else 'unit'

SEED = 20260710
NSIM = 500
STOP_S = 1.75
CONTROL_FACTORS = [0.87, 0.93, 1.07, 1.13]
ANCHORS = [26.0, 52.0, 104.0] if WIN_MIN >= 10.0 else [1.0, 3.0, 26.0, 52.0, 104.0]

# ---- load artifacts of record ----
base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)          # us
errs = base['errs'].astype(np.float64)         # us
toas = base['toas'].astype(np.float64)         # days
names = list(base['names']); fmap = list(base['fmap'])
fitted_names = [names[p] for p in fmap]
N = res0.size

J = np.load(JAC_PATH)                          # us per parameter unit, N x 28
rank_json = json.load(open('carrier_projection_rank.json'))
OMS = [rank_json['carriers_rad_per_day'][k] for k in ('Omega_in', 'Omega_out', 'Omega_dif')]

# stop rule: carrier consistency (recompute from periods, match to 1e-9 rad/day)
P_i = rank_json['periods_days']['P_inner']; P_o = rank_json['periods_days']['P_outer']
W_re = [2*np.pi/P_i, 2*np.pi/P_o, abs(2*np.pi/P_i - 2*np.pi/P_o)]
assert all(abs(a-b) < 1e-9 for a, b in zip(OMS, W_re)), 'carrier mismatch vs artifacts'

sw = 1.0/errs
t = toas

# ---- 10.7e E2: Request-8 dictionary drive template (unit template = 1) ----
if TEMPLATE == 'phys':
    GMSUN = 1.32712440018e20; C2 = 8.9875517873681764e16
    M_A, M_B, M_C = 1.43781441, 0.19753639, 0.41010271          # Msun, Nutimo construction
    P_in_s, P_out_s = P_i*86400.0, P_o*86400.0
    a_in = (GMSUN*(M_A+M_B)*P_in_s**2/(4*np.pi**2))**(1.0/3.0)
    a_out = (GMSUN*(M_A+M_B+M_C)*P_out_s**2/(4*np.pi**2))**(1.0/3.0)
    E_IN = math.hypot(6.9373227475173604e-4, -8.5950920080896147e-5)   # eta_p, kappa_p
    E_OUT = math.hypot(3.5114692182346233e-2, -3.5249206779379461e-3)  # eta_b, kappa_b
    BETA_A = M_B/(M_A+M_B); EPS = a_in/a_out
    OM_S = [w/86400.0 for w in OMS]                              # rad/s
    D_K = [GMSUN*M_B/(a_in*C2)*E_IN/OM_S[0]*1e6,
           GMSUN*M_C/(a_out*C2)*E_OUT/OM_S[1]*1e6,
           GMSUN*M_C/(a_out*C2)*BETA_A*EPS/OM_S[2]*1e6]          # us
    TASC_P, TASC_B = -575.13824374291314, -262.10186433770934    # internal days
    PH = [-OMS[0]*TASC_P - np.pi/2.0, -OMS[1]*TASC_B - np.pi/2.0]
    PH.append(PH[0] - PH[1])                                     # 10.5 phase lock
    AMPW = [D_K[k]*np.exp(1j*PH[k]) for k in range(3)]
    print('phys template: D_k = (%.4g, %.4g, %.4g) us (in, out, dif); phases locked to tasc'
          % tuple(D_K))
else:
    AMPW = [1.0+0j, 1.0+0j, 1.0+0j]

# ---- nuisance block: column-normalized weighted Jacobian + constant column ----
A = sw[:, None] * J
colnorm = np.linalg.norm(A, axis=0)
An = A / colnorm[None, :]
cvec = sw / np.linalg.norm(sw)
B0 = np.column_stack([An, cvec])               # N x 29

# ---- 10.7e E3: red-noise Fourier nuisance (exact sinusoids, linear by construction) ----
if M_RN > 0:
    T_SPAN = float(t.max() - t.min())
    fcols = []
    for j in range(1, M_RN + 1):
        arg = 2.0*np.pi*j*(t - t.min())/T_SPAN
        fcols += [np.cos(arg), np.sin(arg)]
    FB = sw[:, None]*np.column_stack(fcols)
    FB = FB/np.linalg.norm(FB, axis=0, keepdims=True)
    B0 = np.column_stack([B0, FB])             # N x (29 + 2M)
    print('red-noise marginalization: %d Fourier pairs appended (f = 1/T .. %d/T, T = %.2f d)'
          % (M_RN, M_RN, T_SPAN))

U0, s0, V0t = np.linalg.svd(B0, full_matrices=False)
rel0 = s0 / s0[0]
rank0 = int(np.sum(rel0 > SV_CUT))
Q0 = U0[:, :rank0]
if SV_CUT > 1e-10:
    print('nuisance span truncated at rel SV %g -> rank %d/%d (10.7c protocol)'
          % (SV_CUT, rank0, B0.shape[1]))

def proj_out(X):
    return X - Q0 @ (Q0.T @ X)

y = sw * res0
y_perp = proj_out(y)
chi2_nuis = float(y_perp @ y_perp)
dof_nuis = N - rank0
s_scale = math.sqrt(chi2_nuis / dof_nuis)
print('nuisance block rank %d/29; chi2_red(29-par) = %.4f; noise scale s = %.4f'
      % (rank0, chi2_nuis/dof_nuis, s_scale))
if s_scale > STOP_S:
    raise SystemExit('STOP RULE: s = %.3f > %.2f -- red-noise modeling required' % (s_scale, STOP_S))

# ---- signal columns ----
def sig_cols_raw(tau, oms):
    """(x_cY, x_beta) raw us columns; template weights AMPW (unit or phys)."""
    x_c = np.zeros(N); x_b = np.zeros(N)
    for w, aw in zip(oms, AMPW):
        x_c += aw.real*np.cos(w*t) - aw.imag*np.sin(w*t)
        a = aw/(1.0 + 1j*w*tau)
        x_b += a.real*np.cos(w*t) - a.imag*np.sin(w*t)
    return x_c, x_b

def comparator_cols_raw(oms):
    """complex degree-1 per-carrier: A_k = c0 + c1 (i w_k) -> 4 real columns."""
    c0r = np.zeros(N); c0i = np.zeros(N); c1r = np.zeros(N); c1i = np.zeros(N)
    for w in oms:
        c0r += np.cos(w*t)
        c0i += -np.sin(w*t)
        c1r += -w*np.sin(w*t)
        c1i += -w*np.cos(w*t)
    return np.column_stack([c0r, c0i, c1r, c1i])

def joint_fit_at(tau, yv_perp, oms):
    """FWL joint fit at fixed tau. Returns dict incl. beta_hat, sigma_beta, chi2."""
    x_c, x_b = sig_cols_raw(tau, oms)
    Xs = np.column_stack([sw*x_c, sw*x_b])
    X2 = proj_out(Xs)
    M = X2.T @ X2
    v = X2.T @ yv_perp
    theta = np.linalg.solve(M, v)
    r = yv_perp - X2 @ theta
    Minv = np.linalg.inv(M)
    return {
        'c_Y_hat': float(theta[0]), 'beta_hat': float(theta[1]),
        'sigma_c_Y': s_scale*math.sqrt(Minv[0, 0]), 'sigma_beta': s_scale*math.sqrt(Minv[1, 1]),
        'chi2': float(r @ r), 'X2': X2, 'Minv': Minv,
        'survival_beta_col': float(np.linalg.norm(X2[:, 1])/np.linalg.norm(Xs[:, 1])),
    }

def phi(x):
    return 0.5*(1.0 + math.erf(x/math.sqrt(2.0)))

def u95_of(beta_hat, sigma):
    """flat-prior Gaussian posterior on beta: P(|beta| <= u) = 0.95."""
    def mass(u):
        return phi((u-beta_hat)/sigma) - phi((-u-beta_hat)/sigma)
    lo, hi = 0.0, abs(beta_hat) + 10.0*sigma
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if mass(mid) < 0.95:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)

# ---- grids (pre-registered; 10.7d extends the window to [1, 327] d) ----
npts = 49 if WIN_MIN >= 10.0 else 61
quoted_grid = sorted(set(np.geomspace(WIN_MIN, 327.0, npts).tolist() + ANCHORS))
diag_grid = [x for x in np.geomspace(min(1.0, WIN_MIN), 1000.0, 61).tolist()
             if not (WIN_MIN <= x <= 327.0)]

# ---- real-data scan ----
def scan(yv_perp, grid, oms, keep_X2=False):
    out = []
    for tau in grid:
        f = joint_fit_at(tau, yv_perp, oms)
        row = {k: f[k] for k in ('c_Y_hat', 'beta_hat', 'sigma_c_Y', 'sigma_beta',
                                 'chi2', 'survival_beta_col')}
        row['tau'] = float(tau)
        row['z'] = row['beta_hat']/row['sigma_beta']
        row['u95'] = u95_of(row['beta_hat'], row['sigma_beta'])
        if keep_X2:
            row['_X2'] = f['X2']; row['_Minv'] = f['Minv']
        out.append(row)
    return out

rows_q = scan(y_perp, quoted_grid, OMS, keep_X2=True)
rows_d = scan(y_perp, diag_grid, OMS)
Z_obs = max(abs(r['z']) for r in rows_q)
best = min(rows_q, key=lambda r: r['u95'])
print('real-data scan: Z = max|z| = %.3f; min u95 = %.4f us at tau = %.1f d'
      % (Z_obs, best['u95'], best['tau']))

# ---- null calibration: NSIM white-noise sims through the identical pipeline ----
# beta_hat_sim(tau) = [Minv @ X2^T n]_beta with n ~ N(0, s^2 I) (whitened space);
# X2 already _|_ Q0, so projection of n is a no-op for the signal solve.
rng = np.random.default_rng(SEED)
P_rows = np.vstack([ (r['_Minv'] @ r['_X2'].T)[1] for r in rows_q ])   # ntau x N
sig_vec = np.array([r['sigma_beta'] for r in rows_q])
Nn = rng.standard_normal((N, NSIM)) * s_scale
Zsim = np.max(np.abs(P_rows @ Nn) / sig_vec[:, None], axis=0)          # NSIM
p_global = float(np.mean(Zsim >= Z_obs))
z_thresh_997 = float(np.quantile(Zsim, 0.997))
print('null: median max|z| = %.3f, 99.7%% = %.3f; observed Z = %.3f -> p_global = %.4f'
      % (float(np.median(Zsim)), z_thresh_997, Z_obs, p_global))

# ---- off-carrier controls (each with its own null) ----
controls = {}
for fscale in CONTROL_FACTORS:
    oms_f = [w*fscale for w in OMS]
    rows_f = scan(y_perp, quoted_grid, oms_f, keep_X2=True)
    Zf = max(abs(r['z']) for r in rows_f)
    Pf = np.vstack([(r['_Minv'] @ r['_X2'].T)[1] for r in rows_f])
    sf = np.array([r['sigma_beta'] for r in rows_f])
    Nf = rng.standard_normal((N, NSIM)) * s_scale
    Zsim_f = np.max(np.abs(Pf @ Nf) / sf[:, None], axis=0)
    controls[str(fscale)] = {'Z': float(Zf), 'p': float(np.mean(Zsim_f >= Zf))}
    print('control f=%.2f: Z = %.3f, p = %.4f' % (fscale, Zf, controls[str(fscale)]['p']))

detection = bool(p_global < 0.003 and all(c['p'] > 0.05 for c in controls.values()))

# ---- pre-registered comparator: complex degree-1 per-carrier (4 real params) ----
Xc = proj_out(sw[:, None] * comparator_cols_raw(OMS))
Mc = Xc.T @ Xc
tc = np.linalg.solve(Mc, Xc.T @ y_perp)
rc = y_perp - Xc @ tc
chi2_c1 = float(rc @ rc)
chi2_shared_min = min(r['chi2'] for r in rows_q)
tau_at_min = min(rows_q, key=lambda r: r['chi2'])['tau']
print('comparator: chi2(shared-tau,3par) = %.4f at tau=%.1f vs chi2(complex-deg1,4par) = %.4f; chi2(nuis only) = %.4f'
      % (chi2_shared_min, tau_at_min, chi2_c1, chi2_nuis))

# ---- injection-recovery through the FULL joint pipeline ----
injections = {}
for tau_star in ANCHORS:
    row0 = min(rows_q, key=lambda r: abs(r['tau']-tau_star))
    b_inj = row0['u95']
    _, x_b = sig_cols_raw(tau_star, OMS)
    y_inj_perp = proj_out(sw*(res0 + b_inj*x_b))
    rows_inj = scan(y_inj_perp, quoted_grid, OMS)
    row_at = min(rows_inj, key=lambda r: abs(r['tau']-tau_star))
    row_tau_hat = min(rows_inj, key=lambda r: r['chi2'])
    injections['tau_%g' % tau_star] = {
        'beta_inj_us': b_inj,
        'beta_hat_at_tau_star': row_at['beta_hat'],
        'beta_hat_minus_beta_inj_over_sigma':
            (row_at['beta_hat'] - b_inj)/row_at['sigma_beta'],
        'beta_hat_baseline_at_tau_star': row0['beta_hat'],
        'differenced_recovery_exact':
            abs((row_at['beta_hat'] - row0['beta_hat']) - b_inj) < 1e-8*max(1.0, b_inj),
        'tau_hat_chi2min': row_tau_hat['tau'],
        'coverage_pass': bool(abs(row_at['beta_hat'] - b_inj) < 2*row_at['sigma_beta']),
    }
    print('inject tau*=%g: beta_inj=%.4f -> beta_hat=%.4f (dz=%.2f sigma), tau_hat=%.1f'
          % (tau_star, b_inj, row_at['beta_hat'],
             injections['tau_%g' % tau_star]['beta_hat_minus_beta_inj_over_sigma'],
             row_tau_hat['tau']))

# ---- sequential-refit bias reproduction (diagnostic; 1x and 5x amplitude) ----
# 1x = the pre-registered injection amplitude (u95 ~ 2 sigma: tau is not
# localizable by EITHER method there); 5x (~10 sigma) is an additional pure
# linear-algebra diagnostic where tau IS localizable, so joint-vs-sequential
# tau recovery can be compared meaningfully. Quoted limits are unaffected.
S6 = np.column_stack(sum([[np.cos(w*t), np.sin(w*t)] for w in OMS], []))
B6 = sw[:, None] * S6

def seq_tau_hat(y_inj_weighted):
    """10.7a-style project-then-fit: nuisance lstsq, carrier amps, shared-tau grid."""
    dth, _, _, _ = np.linalg.lstsq(B0, y_inj_weighted, rcond=None)
    r_post = y_inj_weighted - B0 @ dth
    amp6 = np.linalg.solve(B6.T @ B6, B6.T @ r_post)
    Aks = [complex(amp6[0], -amp6[1]), complex(amp6[2], -amp6[3]), complex(amp6[4], -amp6[5])]
    yv = np.array(sum([[a.real, a.imag] for a in Aks], []))
    c = np.array(sum([[1.0, 0.0] for _ in OMS], []))
    best_seq = None
    for tau in np.geomspace(0.01, 2000, 4001):
        m = np.array(sum([[(1.0/(1.0+1j*w*tau)).real, (1.0/(1.0+1j*w*tau)).imag]
                          for w in OMS], []))
        Mm = np.column_stack([c, m])
        coef, _, _, _ = np.linalg.lstsq(Mm, yv, rcond=None)
        c2 = float(np.sum((yv - Mm@coef)**2))
        if best_seq is None or c2 < best_seq[0]:
            best_seq = (c2, float(tau))
    return best_seq[1]

_, x_b52 = sig_cols_raw(52.0, OMS)
b1 = injections['tau_52']['beta_inj_us']
seq_bias = {}
for label, fac in (('amp_1x_u95', 1.0), ('amp_5x_u95', 5.0)):
    b_amp = fac*b1
    yw = sw*(res0 + b_amp*x_b52)
    tau_seq = seq_tau_hat(yw)
    rows_amp = scan(proj_out(yw), quoted_grid, OMS)
    row_min = min(rows_amp, key=lambda r: r['chi2'])
    row_at52 = min(rows_amp, key=lambda r: abs(r['tau']-52.0))
    seq_bias[label] = {
        'beta_inj_us': b_amp,
        'tau_hat_sequential_d': tau_seq,
        'tau_hat_joint_d': row_min['tau'],
        'joint_z_at_52d': row_at52['z'],
        'joint_beta_hat_at_52d': row_at52['beta_hat'],
    }
    print('seq-bias diag %s (beta=%.1f us): tau_hat_seq=%.2f d, tau_hat_joint=%.1f d, joint z(52d)=%.1f'
          % (label, b_amp, tau_seq, row_min['tau'], row_at52['z']))

# ---- nuisance displacement of the tau=52 configuration (for WSL linearization check) ----
row52 = min(rows_q, key=lambda r: abs(r['tau']-52.0))
x_c, x_b = sig_cols_raw(52.0, OMS)
y_minus_sig = y - row52['c_Y_hat']*(sw*x_c) - row52['beta_hat']*(sw*x_b)
coef_n, _, _, _ = np.linalg.lstsq(B0, y_minus_sig, rcond=None)
dtheta_param = (coef_n[:28] / colnorm)          # parameter units
np.save('jointfit_dtheta52%s.npy' % SUFFIX, dtheta_param)

# ---- outputs ----
def strip(rows):
    return [{k: v for k, v in r.items() if not k.startswith('_')} for r in rows]

out = {
    'preregistration': '../notes/REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md',
    'jacobian': JAC_PATH,
    'template': TEMPLATE, 'rn_fourier_pairs': M_RN,
    'seed': SEED, 'nsim': NSIM,
    'noise': {'chi2_red_29par': chi2_nuis/dof_nuis, 's_scale': s_scale,
              'stop_rule_s_max': STOP_S, 'nuisance_rank': rank0},
    'real_data': {'Z_max_abs_z': Z_obs, 'p_global': p_global,
                  'null_z_median': float(np.median(Zsim)), 'null_z_997': z_thresh_997,
                  'detection_criterion_met': detection,
                  'min_u95_us': best['u95'], 'tau_at_min_u95_d': best['tau']},
    'controls_off_carrier': controls,
    'comparator': {'chi2_shared_tau_3par_min': chi2_shared_min,
                   'tau_at_chi2_min_d': tau_at_min,
                   'chi2_complex_deg1_4par': chi2_c1,
                   'chi2_nuisance_only': chi2_nuis,
                   'shared_fits_at_least_as_well': bool(chi2_shared_min <= chi2_c1)},
    'injections': injections,
    'sequential_bias_diagnostic': {'tau_star_d': 52.0, 'note':
        'beyond-preregistration pure-linear-algebra diagnostic; quoted limits unaffected',
        **seq_bias},
    'anchors': {('tau_%g' % a): next({'beta_hat': r['beta_hat'], 'sigma_beta': r['sigma_beta'],
                                      'u95': r['u95'], 'z': r['z']}
                                     for r in rows_q if abs(r['tau']-a) < 1e-9)
                for a in ANCHORS},
    'quoted_window_d': [WIN_MIN, 327.0],
    'convention': ('unit drive Lambda_k F_k = 1; beta, c_Y, u95 in us of common pre-transfer drive amplitude'
                   if TEMPLATE == 'unit' else
                   'Request-8 dictionary drive template (10.7e E2): beta, c_Y, u95 DIMENSIONLESS '
                   '(amplitude of the dynamic alpha_A response); D_k and phases recorded here'),
}
if TEMPLATE == 'phys':
    out['phys_template'] = {'D_k_us': D_K, 'phases_rad': PH,
                            'masses_msun': [M_A, M_B, M_C], 'e_in': E_IN, 'e_out': E_OUT,
                            'a_in_m': a_in, 'a_out_m': a_out,
                            'note': 'O(1) geometric projection factors set to 1; see 10.7e note'}
with open('joint_fit_upper_limit%s.json' % SUFFIX, 'w') as fh:
    json.dump(out, fh, indent=1)

hdr = 'tau_chi_d\tbeta_hat_us\tsigma_beta_us\tu95_us\tz\tsurvival\twindow'
lines = []
for r in strip(rows_q):
    lines.append('%.6g\t%.6g\t%.6g\t%.6g\t%.4f\t%.4f\tquoted' %
                 (r['tau'], r['beta_hat'], r['sigma_beta'], r['u95'], r['z'], r['survival_beta_col']))
for r in strip(rows_d):
    lines.append('%.6g\t%.6g\t%.6g\t%.6g\t%.4f\t%.4f\tdiagnostic' %
                 (r['tau'], r['beta_hat'], r['sigma_beta'], r['u95'], r['z'], r['survival_beta_col']))
with open('beta_limit_curve%s.tsv' % SUFFIX, 'w') as fh:
    fh.write(hdr + '\n' + '\n'.join(lines) + '\n')

print('JOINT_FIT_DONE detection=%s min_u95=%.4f us @ tau=%.1f d' % (detection, best['u95'], best['tau']))
