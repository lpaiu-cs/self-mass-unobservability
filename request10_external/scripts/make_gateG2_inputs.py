#!/usr/bin/env python3
"""Build the gateG2 inputs (amendment G-2): causal templates and anchor
parameters at tau = 2 / 18 / 52 d for the C4-completed live-gate harness.

Templates use the identical recipe as the R6 inputs (sep_common.templates,
zero phase); anchor statistics come from the committed corrected 10.8b
artifact. tau = 2 d is the headline anchor, newly gated."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

import sep_common as sc

inp = sc.load_inputs(HERE)
fit = json.load(open('sep_dynamic/sep_joint_fit_10_8b.json'))

TAUS = (2.0, 18.0, 52.0)
tpl = {('T%d' % t): sc.templates(inp, float(t))[1] for t in (2, 18, 52)}
tpl['Tc'] = sc.templates(inp, 18.0)[0]
# stage the committed static-guard column too, so the WSL harness builds
# its C4 span from the archived inputs alone (no scratch jac_sep/ needed)
tpl['jD'] = inp['jD']

gg = {'amendment': '../notes/REQUEST10_8F_REVIEW_RESPONSE.md (G-2)',
      'sv_cut': sc.SV_CUT_DEFAULT,
      'convention': 'causal (10.8f C1); measurement span guard-kept (C4)',
      'anchors': {}}
for a in TAUS:
    r0 = fit['anchors'][sc.anchor_key(a)]
    gg['anchors'][sc.anchor_key(a)] = {
        'u95_Kdyn': r0['u95_Kdyn'], 'u95_Kstatic': r0['u95_Kstatic'],
        'sigma_fisher': r0['sigma_fisher'], 'beta_hat_lin': r0['beta_hat'],
        'z_lin': r0['z']}

# ---- G-2 addendum: the WORST-PHASE headline template pair ----
# The quoted headline is the 10.8e phase-marginalized worst-phase row at
# tau = 2 d; gate that exact template pair, not only the zero-phase one.
pm = json.load(open('sep_dynamic/sep_phase_marg_10_8e.json'))
TOFF_WP = pm['anchors']['tau_2']['worst_toff_K10']
Tc_wp, Tb_wp = sc.templates(inp, 2.0, t_off=TOFF_WP)
tpl['T2wp'] = Tb_wp
tpl['Tc2wp'] = Tc_wp

nu = sc.build_nuisance(inp)
s, _ = sc.noise_scale(nu['Q'], inp['sw'], inp['res0'])
X = np.column_stack([inp['sw']*Tc_wp, inp['sw']*Tb_wp])
X2 = sc.proj_out(nu['Q'], X)
Minv = np.linalg.inv(X2.T @ X2)
y2 = sc.proj_out(nu['Q'], inp['sw']*inp['res0'])
th = Minv @ (X2.T @ y2)
b_wp = float(th[1])
sig_wp = s*float(np.sqrt(Minv[1, 1]))
u95_K10_wp = sc.u95_of(b_wp, 10.0*sig_wp)
K_static = sc.load_anchor_factors(HERE)['K_static']
# consistency: must reproduce the committed 10.8e headline anchor
ref = pm['anchors']['tau_2']['u95pm_K10']
assert abs(u95_K10_wp - ref)/ref < 1e-6, (u95_K10_wp, ref)
gg['anchors']['tau_2_wp'] = {
    't_off_d': TOFF_WP,
    'u95_Kdyn': u95_K10_wp,
    'u95_Kstatic': sc.u95_of(b_wp, K_static*sig_wp),
    'sigma_fisher': sig_wp, 'beta_hat_lin': b_wp, 'z_lin': b_wp/sig_wp,
    'note': ('headline worst-phase row of sep_phase_marg_10_8e.json; '
             'params recomputed on the pipeline span and asserted against '
             'the committed u95pm_K10')}

np.savez('sep_dynamic/gateG2_inputs.npz', **tpl)
json.dump(gg, open('sep_dynamic/gateG2_params.json', 'w'), indent=1)
print('worst-phase anchor: t_off=%.4f d, u95_K10=%.6g (committed %.6g), '
      'z_lin=%.4f' % (TOFF_WP, u95_K10_wp, ref, b_wp/sig_wp))
print('wrote gateG2_inputs.npz', {k: v.shape for k, v in tpl.items()},
      'and gateG2_params.json (anchors: %s)' % list(gg['anchors']))
