#!/usr/bin/env python3
"""10.8g D4: causal vs advanced quadrature discriminability (template geometry).

Measures, in the estimator's weighted nuisance-projected metric, the cosine
between the CAUSAL lag template T_beta(tau) and the ADVANCED (anticausal
control) template at the same tau and phase -- both before and after
marginalizing the co-fitted instantaneous direction T_cY.  This quantifies
how much leverage the data geometry has on the lag/pole itself (as opposed
to the amplitude of an oscillating Delta): a cosine near +1 means the two
quadrature conventions are nearly indistinguishable for the estimator and
any common excess loads both, which is exactly the behaviour seen in the
phase-marginalized E1/anticausal-control statistics of 10.8e.

Design-only measurement: consumes the committed template columns, TOA times
and errors (weights), and the nuisance span -- never the residuals res0.
No detection rule or limit is touched.

Writes: sep_dynamic/sep_quadrature_overlap_10_8g.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sep_common as sc  # noqa: E402

HERE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

WORST_TOFF_HEADLINE = 104.34942814305832  # d; 10.8e tau=2 worst phase (K10)
TAUS = (2.0, 5.0, 18.0, 52.0, 200.0, 500.0)


def unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def main() -> None:
    inp = sc.load_inputs(HERE)
    nu = sc.build_nuisance(inp)
    Q, sw = nu['Q'], inp['sw']

    out = {
        'purpose': ('cosine overlap of causal vs advanced beta templates in the '
                    'weighted nuisance-projected estimator metric; design-only '
                    '(no residuals consumed)'),
        'nuisance_rank': int(nu['rank']),
        'worst_toff_headline_d': WORST_TOFF_HEADLINE,
        'taus': {},
    }
    for tau in TAUS:
        entry = {}
        for label, toff in (('toff_0', 0.0),
                            ('toff_worst_headline', WORST_TOFF_HEADLINE)):
            Tc, Tb_c = sc.templates(inp, tau, t_off=toff, anticausal=False)
            _, Tb_a = sc.templates(inp, tau, t_off=toff, anticausal=True)
            wc = sc.proj_out(Q, sw * Tc)
            bc = sc.proj_out(Q, sw * Tb_c)
            ba = sc.proj_out(Q, sw * Tb_a)
            cos_full = float(bc @ ba / (np.linalg.norm(bc) * np.linalg.norm(ba)))
            wcn = unit(wc)
            bc2 = bc - wcn * (wcn @ bc)
            ba2 = ba - wcn * (wcn @ ba)
            cos_marg = float(bc2 @ ba2 /
                             (np.linalg.norm(bc2) * np.linalg.norm(ba2)))
            entry[label] = {
                'cos_beta_causal_vs_advanced': cos_full,
                'cos_after_cY_marginalization': cos_marg,
            }
        out['taus'][sc.anchor_key(tau)] = entry
        print('tau %6g d: cos(causal, advanced) = %+0.4f (toff 0), '
              '%+0.4f (worst headline phase); after cY marginalization '
              '%+0.4f / %+0.4f'
              % (tau,
                 entry['toff_0']['cos_beta_causal_vs_advanced'],
                 entry['toff_worst_headline']['cos_beta_causal_vs_advanced'],
                 entry['toff_0']['cos_after_cY_marginalization'],
                 entry['toff_worst_headline']['cos_after_cY_marginalization']),
              flush=True)

    dst = os.path.join(HERE, 'sep_dynamic', 'sep_quadrature_overlap_10_8g.json')
    with open(dst, 'w') as fh:
        json.dump(out, fh, indent=1)
    print('wrote', dst)


if __name__ == '__main__':
    main()
