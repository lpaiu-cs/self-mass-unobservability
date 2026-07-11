#!/usr/bin/env python3
"""Request 10.8e: phase-marginalized estimator (E1) and dictionary-anchored
physical limit (E2). Pre-registered in
notes/REQUEST10_8E_PHASE_MARG_PHYSICAL_ANCHOR.md. Deterministic; no live runs.

Template rotation: the response to a drive cos(w t + phi) is
cos(phi) col_w_c + sin(phi) col_w_s (measured quadrature columns). The pole
transfer adds (g1, g2); a common time-origin t_off adds w*t_off per carrier.
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
K10, K934 = 10.0, 934.0
P_IN = 1.6293990080893948
OMS = {'in': 3.856136695791377, 'out': 0.019199654300698213, 'dif': 3.836937041490679}
ANCHORS = [2.0, 5.0, 18.0, 52.0, 200.0]

base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)
errs = base['errs'].astype(np.float64)
t = base['toas'].astype(np.float64)
N = res0.size
sw = 1.0/errs

J = np.load('finite_jacobian_v2.npy')
jD = np.load('sep_dynamic/col_SEP_D.npz')['dcol'].astype(np.float64)
cols = np.load('sep_dynamic/sep_dynamic_columns.npz')

A = sw[:, None]*J
An = A/np.linalg.norm(A, axis=0)
cvec = sw/np.linalg.norm(sw)
T_SPAN = float(t.max() - t.min())
fc = []
for j in range(1, 31):
    arg = 2.0*np.pi*j*(t - t.min())/T_SPAN
    fc += [np.cos(arg), np.sin(arg)]
FB = sw[:, None]*np.column_stack(fc)
FB = FB/np.linalg.norm(FB, axis=0, keepdims=True)
xjD = sw*jD
B0 = np.column_stack([An, cvec, FB, xjD/np.linalg.norm(xjD)])
U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
Q0 = U0[:, :int(np.sum(s0/s0[0] >= SV_CUT))]

def proj_out(X):
    return X - Q0 @ (Q0.T @ X)

y_perp = proj_out(sw*res0)
S_SCALE = 1.1148

# drive models: unit (E1) and dictionary (E2)
F_UNIT = {'in': (1.0, 0.0), 'out': (1.0, 0.0), 'dif': (1.0, 0.0)}   # (amp, phase0)
TASC_P, TASC_B = -575.13824374291314, -262.10186433770934
PH_IN = -OMS['in']*TASC_P - np.pi/2.0
PH_OUT = -OMS['out']*TASC_B - np.pi/2.0
F_PHYS = {'in': (4.27e-11, PH_IN), 'out': (1.21e-10, PH_OUT),
          'dif': (1.12e-11, PH_IN - PH_OUT)}

def templates(tau, t_off, drive):
    Tb = np.zeros(N); Tc = np.zeros(N)
    for name, w in OMS.items():
        amp, ph0 = drive[name]
        ph = ph0 + w*t_off
        g1 = 1.0/(1.0 + (w*tau)**2)
        g2 = w*tau/(1.0 + (w*tau)**2)
        # pole part: Re[e^{i ph}/(1+iw tau) e^{iwt}] response
        cb = amp*(g1*math.cos(ph) + g2*math.sin(ph))
        sb = amp*(g1*math.sin(ph) - g2*math.cos(ph))
        # response to cos(wt+ph_eff): cos-comp*col_c + sin-comp*col_s with
        # cos(wt+p) = cos p cos wt - sin p sin wt; col_s is the (+pi/2) drive
        # response, i.e. the -sin drive: response = cos p col_c + sin p col_s
        Tb += cb*cols['col_%s_c' % name] + sb*cols['col_%s_s' % name]
        Tc += amp*(math.cos(ph)*cols['col_%s_c' % name] +
                   math.sin(ph)*cols['col_%s_s' % name])
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

TAUS = sorted(set(np.geomspace(2.0, 500.0, 61).tolist() + ANCHORS))
TOFFS = np.linspace(0.0, P_IN, 24, endpoint=False)
rng = np.random.default_rng(SEED)

def full_scan(drive, label, anticausal=False):
    """scan (tau, t_off); return per-tau worst-phase stats + P-rows for nulls."""
    sgn = -1.0 if anticausal else 1.0
    rows = []
    P_rows, sig_all = [], []
    for tau in TAUS:
        worst = None
        zmax_tau = 0.0
        for toff in TOFFS:
            Tc, Tb = templates(tau, toff, drive)
            if anticausal:
                Tc2, Tb2 = templates(-tau, toff, drive)   # g2 -> -g2 equivalent
                Tb = Tb2
            Xs = np.column_stack([sw*Tc, sw*Tb])
            X2 = proj_out(Xs)
            M = X2.T @ X2
            Minv = np.linalg.inv(M)
            th = Minv @ (X2.T @ y_perp)
            sig = S_SCALE*math.sqrt(Minv[1, 1])
            z = float(th[1]/sig)
            zmax_tau = max(zmax_tau, abs(z))
            u = u95_of(float(th[1]), sig)
            if worst is None or u > worst['u95_fisher']:
                worst = {'tau': float(tau), 't_off': float(toff),
                         'beta_hat': float(th[1]), 'sigma_fisher': sig,
                         'u95_fisher': u, 'z': z,
                         'survival': float(np.linalg.norm(X2[:, 1])/np.linalg.norm(Xs[:, 1]))}
            P_rows.append((Minv @ X2.T)[1])
            sig_all.append(sig)
        worst['zmax_over_phase'] = zmax_tau
        rows.append(worst)
    return rows, np.vstack(P_rows), np.array(sig_all)

out = {'preregistration': '../notes/REQUEST10_8E_PHASE_MARG_PHYSICAL_ANCHOR.md',
       'seed': SEED, 'nsim': NSIM, 'n_toff': len(TOFFS)}

# ---- E1: unit drive, phase-marginalized ----
rows_u, P_u, sig_u = full_scan(F_UNIT, 'unit')
Z_u = max(r['zmax_over_phase'] for r in rows_u)
Nn = rng.standard_normal((N, NSIM))*S_SCALE
Zs_u = np.max(np.abs(P_u @ Nn)/sig_u[:, None], axis=0)
p_u = float(np.mean(Zs_u >= Z_u))
rows_ac, P_a, sig_a = full_scan(F_UNIT, 'unit-ac', anticausal=True)
Z_a = max(r['zmax_over_phase'] for r in rows_ac)
Na = rng.standard_normal((N, NSIM))*S_SCALE
p_a = float(np.mean(np.max(np.abs(P_a @ Na)/sig_a[:, None], axis=0) >= Z_a))
print('E1 unit phase-marg: Z = %.3f (p = %.3f); anti-causal Z = %.3f (p = %.3f)'
      % (Z_u, p_u, Z_a, p_a))

# ---- E2: dictionary drive, phase-marginalized ----
rows_p, P_p, sig_p = full_scan(F_PHYS, 'phys')
Z_p = max(r['zmax_over_phase'] for r in rows_p)
Np_ = rng.standard_normal((N, NSIM))*S_SCALE
p_p = float(np.mean(np.max(np.abs(P_p @ Np_)/sig_p[:, None], axis=0) >= Z_p))
print('E2 phys phase-marg: Z = %.3f (p = %.3f)' % (Z_p, p_p))

detection = bool((p_u < 0.003 and p_a > 0.05) or (p_p < 0.003))
out['E1'] = {'Z': Z_u, 'p_global': p_u, 'anticausal_Z': Z_a, 'anticausal_p': p_a}
out['E2'] = {'Z': Z_p, 'p_global': p_p}
out['detection_candidate'] = detection

hdr = ('tau_d\tu95pm_fisher\tu95pm_K10\tu95pm_K934\tworst_toff_d\tsurvival'
       '\tu95phys_fisher\tu95phys_K10')
lines = []
anch = {}
for ru, rp in zip(rows_u, rows_p):
    u_f = ru['u95_fisher']
    u10 = u95_of(ru['beta_hat'], K10*ru['sigma_fisher'])
    u934 = u95_of(ru['beta_hat'], K934*ru['sigma_fisher'])
    up_f = rp['u95_fisher']
    up10 = u95_of(rp['beta_hat'], K10*rp['sigma_fisher'])
    lines.append('%.6g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4f\t%.4g\t%.4g' %
                 (ru['tau'], u_f, u10, u934, ru['t_off'], ru['survival'], up_f, up10))
    if any(abs(ru['tau']-a) < 1e-9 for a in ANCHORS):
        anch['tau_%g' % ru['tau']] = {'u95pm_fisher': u_f, 'u95pm_K10': u10,
                                      'u95pm_K934': u934,
                                      'u95phys_K10': up10,
                                      'worst_toff': ru['t_off'],
                                      'survival_worst_phase': ru['survival']}
with open('sep_dynamic/sep_limit_curve_10_8e.tsv', 'w') as fh:
    fh.write(hdr + '\n' + '\n'.join(lines) + '\n')
out['anchors'] = anch
mn = min((json.loads(json.dumps(a))['u95pm_K10'], k) for k, a in anch.items())
out['min_u95pm_K10'] = {'value': mn[0], 'at': mn[1]}
with open('sep_dynamic/sep_phase_marg_10_8e.json', 'w') as fh:
    json.dump(out, fh, indent=1)

print('anchors (phase-marginalized, K10 | phys K10):')
for a in ANCHORS:
    v = anch['tau_%g' % a]
    print('  tau=%5g: u95pm_K10 = %.4g (Fisher %.4g, K934 %.4g) | beta_phys < %.4g'
          % (a, v['u95pm_K10'], v['u95pm_fisher'], v['u95pm_K934'], v['u95phys_K10']))
print('10_8E_DONE detection=%s' % detection)
