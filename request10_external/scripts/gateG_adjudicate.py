#!/usr/bin/env python3
"""Adjudicate the 10.8f R6 causal live-gate run (sep_gateG.json).

The gates were re-run in WSL (scripts/wsl/wsl_gateG.py) on the corrected
CAUSAL templates with amplitudes from the corrected 10.8b artifact. This
script does NOT rescore the registered criteria -- those verdicts stand as
recorded by the gate run itself:

    G1  (|z_GN - z_lin| <= 0.3):        FAIL at both anchors
    G2a (|dev| <= 0.5 sigma_F):         FAIL at tau=52 d
    G2b (|dev_rel| <= 0.02):            PASS at both anchors

It computes the post-hoc DIAGNOSTIC decomposition: the registered G2a/G2b
deviations reference the LINEAR null point (beta_hat_lin), so the GN-vs-
linear null offset measured by G1 enters them as a common mode. The
differential response

    diff = beta_hat_GN(injected) - beta_hat_GN(null) - beta_inj

isolates the live pipeline's response to the injected signal itself.

Interpretation of record (see notes/REQUEST10_8F_REVIEW_RESPONSE.md, R6):
the live truncated-GN refit displaces beta_hat by the null offset
(+0.34 / +0.60 sigma_F at tau = 18 / 52 d) relative to the frozen linear
pipeline, identically with and without injection; the injected-signal
response is accurate to <= 0.02 sigma_F (detection scale) and <= 0.1%
(limit scale). The offset is carried as a live-refit systematic on
Fisher-scale numbers and is covered by the K_dyn = 10 inflation.

Reference: the superseded pre-C1 (anticausal-template) gate run is in git
history at 52153b3:request10_external/sep_dynamic/sep_gateG.json; its G1
measured the SAME GN null residual with the anticausal quadrature
(offsets +0.012 / +0.019 sigma_F), which is why the swap masked the
displacement.
"""
import json
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

d = json.load(open('sep_dynamic/sep_gateG.json'))
an = d['params']['anchors']

out = {
    'registered_verdicts': {
        'G1': d['G1']['pass'], 'G2a': d['G2a']['pass'], 'G2b': d['G2b']['pass'],
        'criteria': {'G1': '|z_GN - z_lin| <= 0.3',
                     'G2a': '|(beta_hat_GN - beta_hat_lin - beta_inj)| <= 0.5 sigma_F',
                     'G2b': '|(beta_hat_GN - beta_hat_lin - beta_inj)| <= 0.02 * beta_inj'},
        'note': 'verdicts as recorded by the gate run; NOT rescored here',
    },
    'diagnostics': {},
}

for a in ('tau_18', 'tau_52'):
    sF = an[a]['sigma_fisher']
    z_gn = d['G1'][a]['z_GN']
    b_gn_null = z_gn * sF
    null_off = z_gn - d['G1'][a]['z_lin']
    g2a = d['G2a'][a]
    g2b = d['G2b'][a]
    diff_a = (g2a['beta_hat_GN'] - b_gn_null - g2a['beta_inj']) / sF
    diff_b = (g2b['beta_hat_GN'] - b_gn_null - g2b['beta_inj']) / g2b['beta_inj']
    out['diagnostics'][a] = {
        'sigma_fisher': sF,
        'null_offset_sigma_F': null_off,
        'G2a_registered_dev_sigma_F': g2a['deviation_sigma_fisher'],
        'G2a_differential_dev_sigma_F': diff_a,
        'G2b_registered_dev_rel': g2b['deviation_relative'],
        'G2b_differential_dev_rel': diff_b,
        'decomposition_check_G2a': {
            'null_offset_plus_differential': null_off + diff_a,
            'registered': g2a['deviation_sigma_fisher'],
        },
    }
    print('%s: null offset %+.4f sigma_F | G2a diff %+.4f sigma_F '
          '(registered %+.4f) | G2b diff %+.5f rel (registered %+.5f)'
          % (a, null_off, diff_a, g2a['deviation_sigma_fisher'],
             diff_b, g2b['deviation_relative']))

offs = [abs(out['diagnostics'][a]['null_offset_sigma_F']) for a in ('tau_18', 'tau_52')]
out['adjudication'] = {
    'finding': ('single common-mode explanation: the GN truncated refit of the '
                'null data displaces beta_hat by the G1 offset; G2a/G2b excesses '
                'over their differential values equal that offset to <= 0.02 '
                'sigma_F at both anchors'),
    'live_refit_systematic_sigma_F_max': max(offs),
    'impact': {
        'detection_tests': ('z-scores carry a +/- %.2f sigma_F live-refit '
                            'systematic; the registered E1/E2 non-detections '
                            '(min p = 0.26 with an equally loud anticausal '
                            'control) are unaffected' % max(offs)),
        'upper_limits': ('worst-case linear addition inflates a Fisher-only '
                         'u95 by <= %.0f%%; the K_dyn = 10 headline scaling '
                         'covers this with an order of magnitude to spare; '
                         'Fisher-only numbers are quoted as statistical-only '
                         'for this reason' % (100.0 * max(offs) / 2.0)),
    },
    'status': ('G1/G2a recorded as FAIL under the registered absolute criteria; '
               'no criterion is redefined post hoc. The template-shape '
               'validation (G2b, and the differential responses) PASSES; the '
               'absolute-null agreement does not.'),
}

with open('sep_dynamic/sep_gateG_adjudication.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('wrote sep_dynamic/sep_gateG_adjudication.json | max null offset %.3f sigma_F' % max(offs))
