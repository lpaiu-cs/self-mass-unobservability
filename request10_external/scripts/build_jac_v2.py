#!/usr/bin/env python3
"""Assemble finite_jacobian_v2.npy: the 10.7a Jacobian with the 7 planet-block
(*_extra1) columns replaced by the re-derived small-step columns in jac_v2/.

Rationale recorded in jac_v2/rederive_report.json and joint_fit_summary.md:
the 10.7a planet steps (0.3 sigma_MCMC) were nonperturbative (up to 2.8 planet
periods), so those columns were secants, not derivatives.
"""
import glob
import json
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

J = np.load('finite_jacobian.npy').copy()
meta = json.load(open('finite_jacobian_meta.json'))
cols = list(meta['columns'])
steps = list(meta['abs_steps'])

replaced = {}
for f in sorted(glob.glob('jac_v2/col_*.npz')):
    m = re.match(r'col_(\d+)_(.+)\.npz$', os.path.basename(f))
    j, pname = int(m.group(1)), m.group(2)
    assert cols[j] == pname, (j, pname, cols[j])
    d = np.load(f)
    J[:, j] = d['dcol']
    steps[j] = float(d['h_abs'])
    replaced[pname] = {'h_abs': float(d['h_abs']), 'maxdev_us': float(d['maxdev_us']),
                       'halfstep_lin_dev': float(d['halfstep_lin_dev'])}
assert len(replaced) == 7, list(replaced)

np.save('finite_jacobian_v2.npy', J)
meta_v2 = dict(meta)
meta_v2['abs_steps'] = steps
meta_v2['definition'] = (meta['definition'] +
    '; v2: the 7 *_extra1 columns re-derived with small (perturbative) steps, '
    'see jac_v2/rederive_report.json')
meta_v2['v2_replaced_columns'] = replaced
with open('finite_jacobian_v2_meta.json', 'w') as fh:
    json.dump(meta_v2, fh, indent=1)
print('finite_jacobian_v2.npy written; replaced:', ', '.join(sorted(replaced)))
for k, v in sorted(replaced.items()):
    print('  %-14s h=%-10g maxdev=%8.3f us  halfstep_lin_dev=%.4f'
          % (k, v['h_abs'], v['maxdev_us'], v['halfstep_lin_dev']))
