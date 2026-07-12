#!/usr/bin/env python3
"""Request 10.7a external Nutimo pilot - Stage 2 (rank gate) + Stage 3 (injection).

Pure linear algebra on the finite Jacobian computed by jacobian_runner.py.

Conventions (documented in configuration_manifest.json):
- t = internal TOA epochs in days (timeshift-subtracted MJD), as returned by Get_toas().
- Carriers: Omega_in = 2pi/P_i, Omega_out = 2pi/P_o, Omega_c = |Omega_in - Omega_out| [rad/day].
- Carrier basis S (ntoa x 6): [cos(Wk t), sin(Wk t)] for k = in, out, c.
- Complex amplitude convention: signal = Re[A_k exp(i Wk t)] = Re(A_k) cos - Im(A_k) sin.
- Weighted metric: W = diag(1/err^2), err from tim file (us).
- Transfer: G(z) = c_Y + beta/(1 + tau_chi z), z_k = i Omega_k. Test columns use c_Y=0, beta=1,
  unit drive amplitude and zero drive phase per carrier (Lambda_k F_k = 1), several tau_chi.
"""
import os, json, glob
import numpy as np

os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
OUT = 'artifacts'
os.makedirs(OUT, exist_ok=True)

base = np.load('baseline_planetGR.npz', allow_pickle=True)
names = list(base['names']); fmap = list(base['fmap'])
fitted_names = [names[p] for p in fmap]
res0, errs, toas = base['res'], base['errs'], base['toas']
N = res0.size; nfit = len(fmap)

# ---- assemble Jacobian ----
cols, hs, missing = {}, {}, []
for j in range(nfit):
    f = 'jac/col_%02d_%s.npz' % (j, fitted_names[j])
    if not os.path.exists(f): missing.append(j); continue
    d = np.load(f); cols[j] = d['dcol']; hs[j] = float(d['h_abs'])
if missing:
    print('MISSING columns:', [(j, fitted_names[j]) for j in missing]); raise SystemExit(1)
J = np.column_stack([cols[j] for j in range(nfit)])  # us per param-unit
np.save(os.path.join(OUT, 'finite_jacobian.npy'), J)
meta = {'columns': fitted_names, 'abs_steps': [hs[j] for j in range(nfit)],
        'units': 'us per parameter unit', 'ntoa': int(N),
        'definition': 'central difference (res(p+h)-res(p-h))/(2h) around parfile-planetGR-max-bestfit'}
with open(os.path.join(OUT, 'finite_jacobian_meta.json'), 'w') as fh: json.dump(meta, fh, indent=1)
print('J assembled', J.shape)

# ---- carriers ----
P_i, P_o = None, None
for line in open('parfile-planetGR-max-bestfit'):
    t_ = line.split()
    if len(t_) >= 2 and t_[0] == 'period_i': P_i = float(t_[1])
    if len(t_) >= 2 and t_[0] == 'period_o': P_o = float(t_[1])
W_in, W_out = 2*np.pi/P_i, 2*np.pi/P_o
W_c = abs(W_in - W_out)
OMS = [W_in, W_out, W_c]; OMNAMES = ['Omega_in', 'Omega_out', 'Omega_dif']
t = toas.astype(np.float64)
S = np.column_stack(sum([[np.cos(w*t), np.sin(w*t)] for w in OMS], []))  # N x 6

sw = 1.0/errs
A = sw[:, None] * J          # weighted Jacobian
B = sw[:, None] * S          # weighted carrier basis

# ---- nuisance-span basis via SVD (rank-aware) ----
# Column-normalize first: raw columns are 'us per parameter unit' with unit scales
# spanning ~16 decades, which would make an SVD cutoff drop genuinely independent
# directions. Normalization changes column scaling only, not the span.
colnorm = np.linalg.norm(A, axis=0)
An = A / colnorm[None, :]
Uj, sj, Vjt = np.linalg.svd(An, full_matrices=False)
rel = sj / sj[0]
rank_J = int(np.sum(rel > 1e-10))
QJ = Uj[:, :rank_J]
rank_sensitivity = {('1e-%d' % k): int(np.sum(rel > 10.0**(-k))) for k in (6, 8, 10, 12)}
print('J (column-normalized) singular values: %.3e .. %.3e; numerical rank %d/%d; sensitivity %s'
      % (sj[0], sj[-1], rank_J, nfit, rank_sensitivity))

# ---- principal cosines between span(J) and 6D carrier space ----
QS, _ = np.linalg.qr(B)
M = QJ.T @ QS
cosines = np.linalg.svd(M, compute_uv=False)  # 6 principal cosines, descending
# per-carrier-coordinate absorption
absorb_cols = []
for c in range(6):
    b = B[:, c]; b = b/np.linalg.norm(b)
    absorb_cols.append(float(np.linalg.norm(QJ.T @ b)))
# carrier amplitudes of each Jacobian column (6 x nfit), columns normalized in W-metric
G6 = np.linalg.solve(B.T @ B, B.T @ A)                  # 6 x nfit carrier amps per unit param
Cn = G6 / np.linalg.norm(A, axis=0, keepdims=True)
sv_C = np.linalg.svd(Cn, compute_uv=False)

th_list = [0.9, 0.99, 0.999, 0.9999]
rank_eff = {str(th): int(np.sum(cosines >= th)) for th in th_list}
rank_json = {
  'carriers_rad_per_day': dict(zip(OMNAMES, OMS)),
  'periods_days': {'P_inner': P_i, 'P_outer': P_o},
  'nonresonance': {'Omega_in/Omega_out': W_in/W_out,
                   'excluded_cases': ['1:1', '2:1', '1:2'], 'pass': bool(abs(W_in/W_out-1)>0.05 and abs(W_in/W_out-2)>0.05 and abs(W_in/W_out-0.5)>0.05)},
  'principal_cosines_spanJ_vs_carrier6D': [float(x) for x in cosines],
  'per_coordinate_absorption': dict(zip(['cos_in','sin_in','cos_out','sin_out','cos_dif','sin_dif'], absorb_cols)),
  'effective_rank_at_cosine_threshold': rank_eff,
  'svd_of_normalized_carrier_projection_of_J': [float(x) for x in sv_C],
  'jacobian_numerical_rank': rank_J, 'rank_at_relative_sv_cutoff': rank_sensitivity,
  'nfit': nfit, 'ntoa': int(N),
  'gate_rule': 'collapse if effective rank = 6/6 (all six carrier coordinates absorbable by finite nuisance span)',
}
with open(os.path.join(OUT, 'carrier_projection_rank.json'), 'w') as fh: json.dump(rank_json, fh, indent=1)
print('principal cosines:', np.round(cosines, 6))
print('rank_eff:', rank_eff)

# ---- dynamic-chi test columns ----
def chi_col(tau, c_Y=0.0, beta=1.0):
    amps = [c_Y + beta/(1.0 + 1j*w*tau) for w in OMS]
    s = np.zeros(N)
    for A_k, w in zip(amps, OMS):
        s += A_k.real*np.cos(w*t) - A_k.imag*np.sin(w*t)
    return s, amps

TAUS = [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 26.0, 52.0, 100.0, 200.0, 327.0]
rows, span_res = [], {}
for tau in TAUS:
    s, amps = chi_col(tau)
    b = sw*s; b = b/np.linalg.norm(b)
    ab = float(np.linalg.norm(QJ.T @ b))
    span_res[str(tau)] = {'absorption': ab, 'survival': float(np.sqrt(max(0.0, 1-ab*ab))),
                          'amps': [[a.real, a.imag] for a in amps]}
    rows.append((tau, s))
    print('tau=%7.2f d  absorption=%.6f  survival=%.4f' % (tau, ab, np.sqrt(max(0.0,1-ab*ab))))
hdr = ['toa_internal_days'] + ['tau_%g' % tau for tau, _ in rows]
tab = np.column_stack([t] + [s for _, s in rows])
np.savetxt(os.path.join(OUT, 'dynamic_chi_test_column.tsv'), tab, delimiter='\t',
           header='\t'.join(hdr), comments='')
with open(os.path.join(OUT, 'dynamic_chi_span_test.json'), 'w') as fh:
    json.dump({'convention': 'unit drive, zero phase, c_Y=0, beta=1; absorption = |Q_J^T b| of unit W-normalized column',
               'results': span_res}, fh, indent=1)

# ---- Stage 3: minimal synthetic injection with linearized refit ----
def carrier_amps_of(r):
    x = np.linalg.solve(B.T @ B, B.T @ (sw*r))
    return [complex(x[0], -x[1]), complex(x[2], -x[3]), complex(x[4], -x[5])]
    # note: r = a cos + b sin = Re[(a - i b) e^{iwt}] -> A = a - i b

def fit_shared_tau(Aks):
    y = np.array(sum([[a.real, a.imag] for a in Aks], []))
    def model(cY, beta, tau):
        m = [cY + beta/(1.0+1j*w*tau) for w in OMS]
        return np.array(sum([[a.real, a.imag] for a in m], []))
    best = None
    for tau in np.geomspace(0.01, 2000, 4001):
        Mm = np.column_stack([model(1,0,tau), model(0,1,tau)])  # linear in (cY, beta)
        coef, _, _, _ = np.linalg.lstsq(Mm, y, rcond=None)
        r = y - Mm@coef; c2 = float(r@r)
        if best is None or c2 < best[0]: best = (c2, tau, float(coef[0]), float(coef[1]))
    return {'chi2': best[0], 'tau_hat': best[1], 'c_Y_hat': best[2], 'beta_hat': best[3], 'npar': 3}

def fit_real_deriv(Aks, deg):
    y = np.array(sum([[a.real, a.imag] for a in Aks], []))
    colsm = []
    for p in range(deg+1):
        m = [(1j*w)**p for w in OMS]
        colsm.append(np.array(sum([[a.real, a.imag] for a in m], [])))
    Mm = np.column_stack(colsm)
    coef, _, _, _ = np.linalg.lstsq(Mm, y, rcond=None)
    r = y - Mm@coef
    return {'chi2': float(r@r), 'coef': [float(c) for c in coef], 'npar': deg+1}

def fit_complex_poly(Aks, deg):
    yc = np.array(Aks)
    Mc = np.column_stack([(1j*np.array(OMS))**p for p in range(deg+1)])
    coef, _, _, _ = np.linalg.lstsq(Mc, yc, rcond=None)
    r = yc - Mc@coef
    return {'chi2': float(np.sum(np.abs(r)**2)), 'coef': [[c.real, c.imag] for c in coef], 'npar': 2*(deg+1)}

INJ_TAUS = [2.0, 5.0, 10.0, 26.0, 52.0, 100.0, 200.0, 327.0]
AMP = 5.0  # us notional scaling of the unit-drive chi column
inj_report = {'amplitude_us': AMP, 'convention': 'r_inj = AMP * chi_col(tau); linearized refit on finite Jacobian span (W-metric)'}
for tau in INJ_TAUS:
    s, amps0 = chi_col(tau)
    r_inj = AMP*s
    dthn, _, _, _ = np.linalg.lstsq(An, sw*r_inj, rcond=None)
    fitted = An@dthn
    r_post = r_inj - fitted/sw
    ab = np.linalg.norm(fitted)/np.linalg.norm(sw*r_inj)
    A_pre = carrier_amps_of(r_inj); A_post = carrier_amps_of(r_post)
    fs = fit_shared_tau(A_post)
    f1 = fit_real_deriv(A_post, 1); f2 = fit_real_deriv(A_post, 2); fc1 = fit_complex_poly(A_post, 1)
    inj_report['tau_%g' % tau] = {
        'injected_amps': [[a.real, a.imag] for a in A_pre],
        'postrefit_amps': [[a.real, a.imag] for a in A_post],
        'absorbed_fraction_W': float(ab),
        'postfit_over_prefit_norm': float(np.linalg.norm(sw*r_post)/np.linalg.norm(sw*r_inj)),
        'shared_tau_fit': fs,
        'real_deriv_deg1': f1, 'real_deriv_deg2': f2, 'complex_poly_deg1': fc1,
        'tau_recovered_ok': bool(0.5 < fs['tau_hat']/tau < 2.0),
        'shared_beats_deg2': bool(fs['chi2'] < f2['chi2']),
        'shared_beats_complex1': bool(fs['chi2'] < fc1['chi2']),
    }
    print('INJ tau=%g: absorbed=%.4f tau_hat=%.3f chi2 shared/deg2/cplx1 = %.3e/%.3e/%.3e' %
          (tau, ab, fs['tau_hat'], fs['chi2'], f2['chi2'], fc1['chi2']))
with open(os.path.join(OUT, 'synthetic_injection_recovery.json'), 'w') as fh:
    json.dump(inj_report, fh, indent=1)
print('ANALYSIS_DONE')
