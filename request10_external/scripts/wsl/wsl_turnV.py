#!/usr/bin/env python3
"""10.8d Stage V: measure the live wrap response to turn-alias spinfreq shifts.

delta_f in {0.5, 1.0}/T_span cycles/day via the fitted-parameter interface
(spinfreq = fitted index 4). Exports the measured residual-change columns for
offline comparison against the analytic wrapped-ramp model.
"""
import json
import os
import sys
import time

import numpy as np

os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
sys.path.insert(0, os.path.expanduser('~/work/nutimo_pilot/nutimo/src'))
from python_Fittriple_interface import PyFittriple

base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)
t = base['toas'].astype(np.float64)
scales = base['scales']
fmap = list(base['fmap'])
names = list(base['names'])
fitted = [names[p] for p in fmap]
j_f = fitted.index('spinfreq')
sc_f = float(scales[fmap[j_f]])
T_SPAN = float(t.max() - t.min())

fit = PyFittriple('parfile-planetGR-max-bestfit',
                  '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
fit.Compute_lnposterior()
print('constructed; spinfreq fitted idx %d scale %.4g; T=%.2f d' % (j_f, sc_f, T_SPAN), flush=True)

out = {}
cols = {}
for n in (0.5, 1.0):
    df = n/T_SPAN
    sh = np.zeros(len(fmap))
    sh[j_f] = df/sc_f
    t0 = time.perf_counter()
    fit.Set_fitted_parameter_relativeshifts(sh)
    fit.Compute_lnposterior()
    r = fit.Get_time_residuals().copy()
    d = r - res0
    cols['dres_n%g' % n] = d
    out['n%g' % n] = {'delta_f_cpd': df, 'maxdev_us': float(np.abs(d).max()),
                      'rms_us': float(np.sqrt(np.mean(d**2))),
                      'seconds': time.perf_counter()-t0}
    print('n=%g: df=%.4g c/d  maxdev=%.1f us  rms=%.1f us (%.0fs)' %
          (n, df, out['n%g' % n]['maxdev_us'], out['n%g' % n]['rms_us'],
           out['n%g' % n]['seconds']), flush=True)

np.savez('turnV_columns.npz', **cols)
with open('turnV_meta.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('TURNV_DONE', flush=True)
