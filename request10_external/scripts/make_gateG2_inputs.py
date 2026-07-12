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
np.savez('sep_dynamic/gateG2_inputs.npz', **tpl)

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
json.dump(gg, open('sep_dynamic/gateG2_params.json', 'w'), indent=1)
print('wrote gateG2_inputs.npz', {k: v.shape for k, v in tpl.items()},
      'and gateG2_params.json (anchors: %s)' % list(gg['anchors']))
