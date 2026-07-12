#!/usr/bin/env python3
"""Adjudicate the live-gate record: the 10.8f R6 causal run (sep_gateG.json)
and the amendment G-2 re-run with the C4-complete harness (sep_gateG2.json).

R6 verdicts under the registered criteria (recorded, never rescored):
G1 FAIL at both anchors, G2a FAIL at tau = 52 d, G2b PASS.

This script computes the two-layer diagnosis:

1. DIFFERENTIAL decomposition (as in the first adjudication): the
   registered G2a/G2b deviations reference the LINEAR null point
   (beta_hat_lin), so the GN-vs-linear null offset measured by G1 enters
   them as a common mode; subtracting the GN null isolates the injected-
   signal response.

2. SPAN decomposition (amendment G-2; this SUPERSEDES the first
   adjudication's attribution): the R6 harness measured beta_hat in a
   truncation-only span (rank 70) -- correction C4 (keep the static-guard
   residual direction) was applied to every pipeline script but never
   propagated into the WSL gate harness. Measuring the SAME baseline
   residuals res0 with NO GN refit in that span reproduces the entire "G1
   offset"; the guard-kept span reproduces z_lin on res0 to machine
   precision. The first adjudication's LS-vs-MAP-walk attribution was
   therefore wrong as stated: the walk exists but its projection onto the
   estimator is at the 0.01-sigma level.

Record of record: sep_gateG2.json (guard-kept harness, anchors tau =
2/18/52 d including the headline anchor, criteria unchanged) -- G1/G2a/G2b
PASS 9/9. The R6 artifact remains recorded with this diagnosis.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

d = json.load(open('sep_dynamic/sep_gateG.json'))
an = d['params']['anchors']

out = {
    'registered_verdicts_R6': {
        'G1': d['G1']['pass'], 'G2a': d['G2a']['pass'], 'G2b': d['G2b']['pass'],
        'criteria': {'G1': '|z_GN - z_lin| <= 0.3',
                     'G2a': '|(beta_hat_GN - beta_hat_lin - beta_inj)| <= 0.5 sigma_F',
                     'G2b': '|(beta_hat_GN - beta_hat_lin - beta_inj)| <= 0.02 * beta_inj'},
        'note': 'verdicts as recorded by the R6 run; NOT rescored here',
    },
    'differential_diagnostics': {},
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
    out['differential_diagnostics'][a] = {
        'sigma_fisher': sF,
        'null_offset_sigma_F': null_off,
        'G2a_registered_dev_sigma_F': g2a['deviation_sigma_fisher'],
        'G2a_differential_dev_sigma_F': diff_a,
        'G2b_registered_dev_rel': g2b['deviation_relative'],
        'G2b_differential_dev_rel': diff_b,
    }
    print('%s: R6 null offset %+.4f sigma_F | G2a diff %+.4f | G2b diff %+.5f'
          % (a, null_off, diff_a, diff_b))

# ---- span decomposition (amendment G-2): rebuild both harness spans and
# measure res0 (no GN) against the pipeline z_lin ----
import sys
sys.path.insert(0, os.path.join(HERE, 'scripts'))
import sep_common as sc

inp = sc.load_inputs(HERE)
sw, res0, t = inp['sw'], inp['res0'], inp['t']
J = np.load('finite_jacobian_v2.npy')
jD = np.load('sep_dynamic/col_SEP_D.npz')['dcol'].astype(np.float64)
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
xn = sw*jD
xn = xn/np.linalg.norm(xn)
B0 = np.column_stack([An, cvec, FB, xn])
U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
k0 = int(np.sum(s0/s0[0] >= 1e-3))
Q_r6 = U0[:, :k0]                                   # R6 harness span
g = xn - Q_r6 @ (Q_r6.T @ xn)
Q_c4 = np.column_stack([Q_r6, g/np.linalg.norm(g)])  # C4-complete span

gi = np.load('sep_dynamic/gateG_inputs.npz')
TC = gi['Tc'].astype(np.float64)


def measure(Q, r_us, Tb):
    Xs = np.column_stack([sw*TC, sw*Tb])
    X2 = Xs - Q @ (Q.T @ Xs)
    Minv = np.linalg.inv(X2.T @ X2)
    yr = sw*r_us
    th = Minv @ (X2.T @ (yr - Q @ (Q.T @ yr)))
    return float(th[1])


out['span_decomposition'] = {'note': ('z(res0) measured with NO GN refit in '
                                      'each harness span, vs the pipeline z_lin')}
for a, key in (('tau_18', 'T18'), ('tau_52', 'T52')):
    ap = an[a]
    Tb = gi[key].astype(np.float64)
    z_r6 = measure(Q_r6, res0, Tb)/ap['sigma_fisher']
    z_c4 = measure(Q_c4, res0, Tb)/ap['sigma_fisher']
    out['span_decomposition'][a] = {
        'z_res0_R6_span': z_r6, 'z_res0_C4_span': z_c4, 'z_lin': ap['z_lin'],
        'span_only_delta_sigma_F': z_r6 - ap['z_lin'],
        'C4_span_reproduces_z_lin_to': z_c4 - ap['z_lin'],
        'R6_G1_offset_recorded': d['G1'][a]['z_GN'] - ap['z_lin'],
        'true_refit_displacement_in_R6_span': d['G1'][a]['z_GN'] - z_r6,
    }
    print('%s: span-only delta %+.4f (recorded offset %+.4f) | '
          'C4 span residual %+.2g | true refit in-span %+.4f'
          % (a, z_r6 - ap['z_lin'], d['G1'][a]['z_GN'] - ap['z_lin'],
             z_c4 - ap['z_lin'], d['G1'][a]['z_GN'] - z_r6))

# ---- gateG2 record of record ----
try:
    g2 = json.load(open('sep_dynamic/sep_gateG2.json'))
    out['gateG2_record_of_record'] = {
        'G1': {a: {'offset_sigma_F': g2['G1'][a]['offset'],
                   'pass': g2['G1'][a]['pass']}
               for a in ('tau_2', 'tau_18', 'tau_52')},
        'G2a': {a: {'dev_sigma_F': g2['G2a'][a]['deviation_sigma_fisher'],
                    'pass': g2['G2a'][a]['pass']}
                for a in ('tau_2', 'tau_18', 'tau_52')},
        'G2b': {a: {'dev_rel': g2['G2b'][a]['deviation_relative'],
                    'pass': g2['G2b'][a]['pass']}
                for a in ('tau_2', 'tau_18', 'tau_52')},
        'all_pass': bool(g2['G1']['pass'] and g2['G2a']['pass'] and g2['G2b']['pass']),
    }
    print('gateG2: all_pass =', out['gateG2_record_of_record']['all_pass'])
except FileNotFoundError:
    out['gateG2_record_of_record'] = 'sep_gateG2.json not present'

# verdict text is DERIVED from the gateG2 output, never hardcoded: a
# missing or failing re-run must show up in this artifact as such
g2rec = out['gateG2_record_of_record']
if isinstance(g2rec, dict):
    max_g1 = max(abs(v['offset_sigma_F']) for v in g2rec['G1'].values())
    max_g2a = max(abs(v['dev_sigma_F']) for v in g2rec['G2a'].values())
    max_g2b = max(abs(v['dev_rel']) for v in g2rec['G2b'].values())
    n_pass = sum(bool(v['pass']) for blk in ('G1', 'G2a', 'G2b')
                 for v in g2rec[blk].values())
    n_tot = sum(len(g2rec[blk]) for blk in ('G1', 'G2a', 'G2b'))
    record = ('sep_gateG2.json: C4-complete harness, anchors tau = 2/18/52 d '
              '(headline anchor directly gated), registered criteria '
              'unchanged -- %d/%d %s; max |G1 offset| %.3f sigma_F, max '
              '|G2a dev| %.3f sigma_F, max |G2b dev| %.3f%%.'
              % (n_pass, n_tot, 'PASS' if g2rec['all_pass'] else
                 'with FAILURES', max_g1, max_g2a, 100.0*max_g2b))
    if g2rec['all_pass']:
        systematic = ('live-refit displacement <= %.3f sigma_F as measured '
                      '(replaces the earlier 0.6 sigma_F carry, which was '
                      'the harness artifact)' % max_g1)
    else:
        systematic = ('gateG2 recorded FAILURES: carry each failed '
                      'deviation at its measured size; the supersession '
                      'claim does NOT hold until adjudicated')
else:
    record = ('sep_gateG2.json NOT PRESENT -- there is no superseding '
              'live-gate record; the R6 verdicts (G1/G2a FAIL) stand')
    systematic = ('R6-era carry applies: the span decomposition bounds the '
                  'harness artifact, but no corrected-harness measurement '
                  'exists')
out['adjudication'] = {
    'finding': ('the R6 G1/G2a failures were a HARNESS defect, not a live-'
                'model property: the R6 measurement span omitted the C4 '
                'guard-residual direction; the span-only delta on res0 (no '
                'GN) reproduces the recorded offsets, and the C4-complete '
                'span reproduces z_lin to machine precision. The first '
                'adjudication attributed the offset to the truncated-GN '
                'LS-vs-MAP walk; that attribution was wrong as stated -- '
                'the walk projects onto the estimator at only the '
                '0.01-sigma level.'),
    'record_of_record': record,
    'systematic_carried': systematic,
    'status': ('R6 verdicts stand as recorded (FAIL under its harness); '
               'amendment G-2 (pre-committed) governs the supersession; '
               'no criterion was redefined.'),
}

with open('sep_dynamic/sep_gateG_adjudication.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('wrote sep_dynamic/sep_gateG_adjudication.json')
