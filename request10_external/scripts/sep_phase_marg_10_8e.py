#!/usr/bin/env python3
"""Request 10.8e phase-marginalized limits (10.8f-corrected re-run).

Corrections (notes/REQUEST10_8F_REVIEW_RESPONSE.md):
C3  per-K worst-phase maxima (the original evaluated K-scaled limits at the
    Fisher-argmax phase, understating 4/5 anchors);
C4  active static-leakage guard;
C5  COMPLETE common-origin scan: t_off over [0, P_out) with inner-phase step
    <= 15 deg (the original scanned [0, P_in), moving the outer phase 1.8 deg),
    made affordable by a 6-projected-column factorization;
C6  span-truncation bracket: the full-rank curve is computed alongside;
C9  the pre-registered detection rule (E1 with anti-causal control; E2 is
    limit-only) replaces the drifted implementation;
C10 noise scale computed; K factors from artifacts.
Deterministic (seed 20260710).
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
ANCHORS = [2.0, 5.0, 18.0, 52.0, 200.0]

inp = sc.load_inputs(HERE)
sw, res0, N, OMS = inp['sw'], inp['res0'], inp['N'], inp['OMS']
KF = sc.load_anchor_factors(HERE)
K10, K934 = KF['K_dyn'], KF['K_static']

TASC_P, TASC_B = -575.13824374291314, -262.10186433770934
PH_IN = -OMS['in']*TASC_P - np.pi/2.0
PH_OUT = -OMS['out']*TASC_B - np.pi/2.0
# free-fall U/c^2 drive amplitudes (leading order, O(1) geometry; 10.8e note)
F_PHYS = {'in': (4.27e-11, PH_IN), 'out': (1.21e-10, PH_OUT),
          'dif': (1.12e-11, PH_IN - PH_OUT)}

TAUS = sorted(set(np.geomspace(2.0, 500.0, 61).tolist() + ANCHORS))
# complete common-origin domain: full outer period, inner-phase step <= 15 deg
DT = inp['P_in']/24.0
TOFFS = np.arange(0.0, inp['P_out'], DT)
print('t_off grid: %d points over [0, %.2f) d (step %.4f d)'
      % (len(TOFFS), inp['P_out'], DT))

rng = np.random.default_rng(sc.SEED)


def factor_scan(Q, label, drive=None, anticausal=False, do_null=True):
    """Factorized (tau x t_off) scan on span Q. Returns per-tau worst-phase
    rows (per-K maxima), the global max|z|, and its null p-value."""
    s_scale, _ = sc.noise_scale(Q, sw, res0)
    C6p, order = sc.projected_sixcols(inp, Q)
    G = C6p.T @ C6p                       # 6x6 Gram
    y_perp = sc.proj_out(Q, sw*res0)
    b6 = C6p.T @ y_perp
    # raw (unprojected, weighted) Gram for survival denominators
    C6w = np.column_stack([sw*inp['cols'][k] for k in order])
    Graw = C6w.T @ C6w
    rows = []
    z_abs_max = 0.0
    p_cells = []                          # per-cell 6-vector for null z
    sig_cells = []
    for tau in TAUS:
        worst = {'fisher': None, 'K10': None, 'K934': None}
        zmax_tau = 0.0
        for toff in TOFFS:
            wc, wb = sc.template_stencil(tau, toff, OMS, drive, anticausal)
            W = np.column_stack([wc, wb])
            M = W.T @ G @ W
            Minv = np.linalg.inv(M)
            th = Minv @ (W.T @ b6)
            sig = s_scale*math.sqrt(Minv[1, 1])
            b = float(th[1])
            z = b/sig
            zmax_tau = max(zmax_tau, abs(z))
            pvec = W @ Minv[:, 1]         # beta-row projector in 6-space
            p_cells.append(pvec)
            sig_cells.append(sig)
            for key, kk in (('fisher', 1.0), ('K10', K10), ('K934', K934)):
                u = sc.u95_of(b, kk*sig)
                if worst[key] is None or u > worst[key]['u95']:
                    surv = math.sqrt(max(0.0, float(wb @ G @ wb))) / \
                        math.sqrt(max(1e-300, float(wb @ Graw @ wb)))
                    worst[key] = {'u95': u, 't_off': float(toff),
                                  'beta_hat': b, 'sigma_fisher': sig,
                                  'survival': surv}
        rows.append({'tau': float(tau), 'zmax': zmax_tau, 'worst': worst})
        z_abs_max = max(z_abs_max, zmax_tau)
    p_global = None
    if do_null:
        Pm = np.vstack(p_cells)           # ncells x 6
        sv = np.array(sig_cells)
        Nn = rng.standard_normal((N, NSIM))*s_scale
        Z6 = C6p.T @ Nn                   # 6 x NSIM
        Zsim = np.max(np.abs(Pm @ Z6)/sv[:, None], axis=0)
        p_global = float(np.mean(Zsim >= z_abs_max))
    print('%s: Z = %.3f%s' % (label, z_abs_max,
          '' if p_global is None else ' (p = %.4f)' % p_global), flush=True)
    return rows, z_abs_max, p_global, s_scale


# ---- spans: truncated+guard (of record) and full-rank (bracket, C6) ----
nu_t = sc.build_nuisance(inp)
nu_f = sc.build_nuisance(inp, sv_cut=0.0, keep_guard=False)
print('spans: truncated rank %d (guard residual %.3g kept) | full rank %d'
      % (nu_t['rank'], nu_t['guard_residual_norm'], nu_f['rank']))

# E1 unit drive, causal, truncated span (headline)
rows_u, Z_u, p_u, s_t = factor_scan(nu_t['Q'], 'E1 unit causal (trunc)')
# anti-causal control on the same span
rows_a, Z_a, p_a, _ = factor_scan(nu_t['Q'], 'E1 anti-causal control',
                                  anticausal=True)
# E2 dictionary drive (limit-only per pre-registration, C9)
rows_p, Z_p, p_p, _ = factor_scan(nu_t['Q'], 'E2 phys causal (trunc)',
                                  drive=F_PHYS)
# full-rank bracket (limits only; no detection role)
rows_f, Z_f, _, s_f = factor_scan(nu_f['Q'], 'E1 unit causal (FULL-rank bracket)',
                                  do_null=False)

# pre-registered detection rule (E1 with control; E2 limit-only)
detection = bool(p_u < 0.003 and p_a > 0.05)

out = {
    'preregistration': '../notes/REQUEST10_8E_PHASE_MARG_PHYSICAL_ANCHOR.md',
    'correction': '../notes/REQUEST10_8F_REVIEW_RESPONSE.md (C3/C4/C5/C6/C9/C10)',
    'seed': sc.SEED, 'nsim': NSIM,
    'n_toff': len(TOFFS), 'toff_domain_d': [0.0, inp['P_out']],
    'K_dyn': K10, 'K_static': K934,
    'noise_scale_trunc': s_t, 'noise_scale_full': s_f,
    'nuisance_rank_trunc': nu_t['rank'], 'nuisance_rank_full': nu_f['rank'],
    'E1': {'Z': Z_u, 'p_global': p_u},
    'anticausal_control': {'Z': Z_a, 'p': p_a, 'quiet': bool(p_a > 0.05)},
    'E2': {'Z': Z_p, 'p_global_reported_not_gating': p_p},
    'detection_candidate': detection,
}

hdr = ('tau_d\tu95pm_fisher\tu95pm_K10\tu95pm_K934\tu95pm_K10_fullrank'
       '\tworst_toff_K10_d\tsurvival_K10\tu95phys_K10')
lines = []
anch = {}
for ru, rf, rp in zip(rows_u, rows_f, rows_p):
    wF, w10, w934 = ru['worst']['fisher'], ru['worst']['K10'], ru['worst']['K934']
    f10 = rf['worst']['K10']
    p10 = rp['worst']['K10']
    lines.append('%.6g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4f\t%.4g' %
                 (ru['tau'], wF['u95'], w10['u95'], w934['u95'], f10['u95'],
                  w10['t_off'], w10['survival'], p10['u95']))
    if any(abs(ru['tau'] - a) < 1e-9 for a in ANCHORS):
        anch[sc.anchor_key(ru['tau'])] = {
            'u95pm_fisher': wF['u95'], 'u95pm_K10': w10['u95'],
            'u95pm_K934': w934['u95'], 'u95pm_K10_fullrank': f10['u95'],
            'u95phys_K10': p10['u95'], 'worst_toff_K10': w10['t_off']}
with open('sep_dynamic/sep_limit_curve_10_8e.tsv', 'w') as fh:
    fh.write(hdr + '\n' + '\n'.join(lines) + '\n')
out['anchors'] = anch
mins = min((a['u95pm_K10'], k) for k, a in anch.items())
curve_min = min((float(l.split('\t')[2]), float(l.split('\t')[0])) for l in lines)
out['min_u95pm_K10'] = {'value': curve_min[0], 'at_tau_d': curve_min[1],
                        'min_over_anchors': {'value': mins[0], 'at': mins[1]}}
with open('sep_dynamic/sep_phase_marg_10_8e.json', 'w') as fh:
    json.dump(out, fh, indent=1)

print('anchors (per-K worst-phase | full-rank bracket | phys):')
for a in ANCHORS:
    v = anch[sc.anchor_key(a)]
    print('  tau=%5g: K10 %.4g (Fisher %.4g, K934 %.4g) | full-rank %.4g | beta_phys < %.4g'
          % (a, v['u95pm_K10'], v['u95pm_fisher'], v['u95pm_K934'],
             v['u95pm_K10_fullrank'], v['u95phys_K10']))
print('curve min u95pm_K10 = %.4g at tau = %g' % curve_min)
print('10_8E_DONE detection=%s' % detection)
