#!/usr/bin/env python3
"""Request 10.8h: red-noise robustness of the free-fall limits.

Pre-registered in notes/REQUEST10_8H_REDNOISE_ROBUSTNESS.md (committed with
this harness, before the run). Mirrors the 10.8e factorized worst-phase
machinery exactly (truncated+guard span, per-span noise scale, 6-projected-
column scan over the full common-origin domain, per-K worst-phase maxima),
varying ONLY the red-noise Fourier order m_rn.

Registered grid:   m_rn in {0, 5, 10, 15, 20, 30, 45, 60}   (30 = of record)
Registered anchors: tau in {2, 5, 18, 52, 200, 500} d
Statistic of record: worst-phase u95 per tier per anchor per m_rn, and
R(tau, m) = u95pm_K10(tau; m) / u95pm_K10(tau; 30).
Stability criterion: max over m in {45, 60} of R(tau, m) <= 1.5.
Integrity stop-rule: the m_rn = 30 row must reproduce the committed 10.8e
values (JSON anchors to rel < 1e-6; tau = 500 tsv row to < 0.5%).

No null simulations; no detection rescoring; deterministic.
Writes: sep_dynamic/sep_rn_robustness_10_8h.json
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

import sep_common as sc  # noqa: E402

M_GRID = (0, 5, 10, 15, 20, 30, 45, 60)
M_RICH = (45, 60)                      # criterion set (richer-than-record)
TAUS = (2.0, 5.0, 18.0, 52.0, 200.0, 500.0)
R_FLAG = 1.5

inp = sc.load_inputs(HERE)
sw, res0, OMS = inp['sw'], inp['res0'], inp['OMS']
KF = sc.load_anchor_factors(HERE)
K10, K934 = KF['K_dyn'], KF['K_static']

DT = inp['P_in'] / 24.0
TOFFS = np.arange(0.0, inp['P_out'], DT)
print('t_off grid: %d points over [0, %.2f) d (step %.4f d)'
      % (len(TOFFS), inp['P_out'], DT), flush=True)
T_SPAN = float(inp['t'].max() - inp['t'].min())
print('span T = %.1f d; outer carrier at j* = T/P_out = %.2f'
      % (T_SPAN, T_SPAN / inp['P_out']), flush=True)


def scan_span(m_rn: int) -> dict:
    """10.8e-identical worst-phase scan on the truncated+guard span with
    m_rn Fourier pairs. Unit drive, causal convention of record."""
    nu = sc.build_nuisance(inp, m_rn=m_rn)
    Q = nu['Q']
    s_scale, chi2red = sc.noise_scale(Q, sw, res0)
    C6p, order = sc.projected_sixcols(inp, Q)
    G = C6p.T @ C6p
    y_perp = sc.proj_out(Q, sw * res0)
    b6 = C6p.T @ y_perp
    C6w = np.column_stack([sw * inp['cols'][k] for k in order])
    Graw = C6w.T @ C6w
    out = {'m_rn': m_rn, 'rank': int(nu['rank']), 'ncols': int(nu['ncols']),
           'noise_scale': float(s_scale), 'chi2_red': float(chi2red),
           'anchors': {}}
    for tau in TAUS:
        worst = {'fisher': None, 'K10': None, 'K934': None}
        zmax = 0.0
        for toff in TOFFS:
            wc, wb = sc.template_stencil(tau, float(toff), OMS)
            W = np.column_stack([wc, wb])
            M = W.T @ G @ W
            Minv = np.linalg.inv(M)
            th = Minv @ (W.T @ b6)
            sig = s_scale * math.sqrt(Minv[1, 1])
            b = float(th[1])
            zmax = max(zmax, abs(b / sig))
            for key, kk in (('fisher', 1.0), ('K10', K10), ('K934', K934)):
                u = sc.u95_of(b, kk * sig)
                if worst[key] is None or u > worst[key]['u95']:
                    surv = math.sqrt(max(0.0, float(wb @ G @ wb))) / \
                        math.sqrt(max(1e-300, float(wb @ Graw @ wb)))
                    worst[key] = {'u95': u, 't_off': float(toff),
                                  'beta_hat': b, 'sigma_fisher': sig,
                                  'survival': surv}
        out['anchors'][sc.anchor_key(tau)] = {
            'u95pm_fisher': worst['fisher']['u95'],
            'u95pm_K10': worst['K10']['u95'],
            'u95pm_K934': worst['K934']['u95'],
            'worst_toff_K10': worst['K10']['t_off'],
            'survival_K10': worst['K10']['survival'],
            'zmax': zmax,
        }
        print('  m_rn=%2d tau=%6g: u95pm_K10 = %.4e (fisher %.4e)  '
              'zmax = %.3f  surv = %.3f'
              % (m_rn, tau, worst['K10']['u95'], worst['fisher']['u95'],
                 zmax, worst['K10']['survival']), flush=True)
    return out


def integrity_check(row30: dict) -> dict:
    """Registered stop-rule: m_rn = 30 must reproduce the committed 10.8e
    values (JSON anchors rel < 1e-6; tau = 500 tsv row < 0.5%)."""
    ref = json.load(open(os.path.join(
        HERE, 'sep_dynamic', 'sep_phase_marg_10_8e.json')))
    checks = {}
    ok = True
    for key, a in ref['anchors'].items():
        mine = row30['anchors'][key]['u95pm_K10']
        rel = abs(mine - a['u95pm_K10']) / a['u95pm_K10']
        checks[key] = {'committed': a['u95pm_K10'], 'reproduced': mine,
                       'rel_dev': rel, 'pass': rel < 1e-6}
        ok = ok and rel < 1e-6
    tsv = open(os.path.join(HERE, 'sep_dynamic',
                            'sep_limit_curve_10_8e.tsv')).read().splitlines()
    hdr = tsv[0].split('\t')
    row500 = None
    for line in tsv[1:]:
        f = line.split('\t')
        if abs(float(f[0]) - 500.0) < 1e-9:
            row500 = dict(zip(hdr, f))
    assert row500 is not None, 'tau=500 row missing from committed tsv'
    committed = float(row500['u95pm_K10'])
    mine = row30['anchors']['tau_500']['u95pm_K10']
    rel = abs(mine - committed) / committed
    checks['tau_500_tsv'] = {'committed': committed, 'reproduced': mine,
                             'rel_dev': rel, 'pass': rel < 5e-3}
    ok = ok and rel < 5e-3
    return {'pass': ok, 'checks': checks}


def main() -> None:
    rows = {}
    for m in M_GRID:
        print('scanning m_rn = %d ...' % m, flush=True)
        rows['m_%d' % m] = scan_span(m)

    integ = integrity_check(rows['m_30'])
    print('integrity (m_rn = 30 reproduces 10.8e): %s'
          % ('PASS' if integ['pass'] else '*** FAIL -- STOP ***'), flush=True)

    ratios = {}
    flags = {}
    for tau in TAUS:
        key = sc.anchor_key(tau)
        base = rows['m_30']['anchors'][key]['u95pm_K10']
        ratios[key] = {'m_%d' % m: rows['m_%d' % m]['anchors'][key]['u95pm_K10'] / base
                       for m in M_GRID}
        rmax_rich = max(ratios[key]['m_%d' % m] for m in M_RICH)
        flags[key] = {'R_max_rich': rmax_rich,
                      'stable': bool(rmax_rich <= R_FLAG)}
        print('%s: R over grid = %s  | R_max{45,60} = %.3f -> %s'
              % (key,
                 {k: round(v, 3) for k, v in ratios[key].items()},
                 rmax_rich,
                 'STABLE' if flags[key]['stable'] else '*** FLAGGED ***'),
              flush=True)

    out = {
        'preregistration': 'notes/REQUEST10_8H_REDNOISE_ROBUSTNESS.md',
        'purpose': ('red-noise robustness of the free-fall worst-phase '
                    'limits: vary the Fourier order m_rn only; 10.8e '
                    'machinery otherwise unchanged'),
        'm_grid': list(M_GRID),
        'criterion': {'set': list(M_RICH), 'R_flag': R_FLAG},
        'n_toff': int(len(TOFFS)),
        'span_T_d': T_SPAN,
        'outer_carrier_j': T_SPAN / inp['P_out'],
        'integrity': integ,
        'spans': rows,
        'ratios_K10_vs_m30': ratios,
        'stability_flags': flags,
    }
    dst = os.path.join(HERE, 'sep_dynamic', 'sep_rn_robustness_10_8h.json')
    with open(dst, 'w') as fh:
        json.dump(out, fh, indent=1)
    print('wrote', dst, flush=True)
    if not integ['pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
