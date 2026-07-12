#!/usr/bin/env python3
"""Request 10.8 R5 resolution: central-difference columns for the FIXED
parameters that Delta physically mixes with (GM/geometry absorbers):
masspar_i, mass_o, apcosi_i, distance.

Parfile-edit + reconstruction route (these are outside the fitted map).
Adaptive step targeting maxdev in [0.02, 30] us (3 attempts, accept last).
No half-step check: these are nuisance-span columns, not signal templates.
"""
import json
import os
import sys
import time

import numpy as np

os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
sys.path.insert(0, os.path.expanduser('~/work/nutimo_pilot/nutimo/src'))
import python_Fittriple_interface as pfi

base = np.load('baseline_planetGR.npz', allow_pickle=True)
res0 = base['res'].astype(np.float64)

PARFILE = 'parfile-planetGR-max-bestfit'
TIMFILE = '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'
LINES = open(PARFILE).read().splitlines()

PARAMS = {                       # name -> initial absolute step
    'masspar_i': 3e-9,
    'mass_o': 3e-9,
    'apcosi_i': 1e-5,
    'distance': 1.0,
}

def base_value(name):
    for ln in LINES:
        tok = ln.split()
        if tok and tok[0] == name:
            return float(tok[1])
    raise SystemExit('%s not in parfile' % name)

def res_at(name, value):
    out = []
    for ln in LINES:
        tok = ln.split()
        if tok and tok[0] == name:
            parts = tok[:]
            parts[1] = '%.16e' % value
            out.append(' '.join(parts))
        else:
            out.append(ln)
    path = 'parfile-r5-tmp'
    with open(path, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    t0 = time.perf_counter()
    fit = pfi.PyFittriple(path, TIMFILE)
    fit.Compute_lnposterior()
    r = fit.Get_time_residuals().copy()
    print('  %s=%+.10e: %.0fs' % (name, value, time.perf_counter()-t0), flush=True)
    del fit
    return r

os.makedirs('jac_sep', exist_ok=True)
report = {}
for name, h0 in PARAMS.items():
    v0 = base_value(name)
    h = float(h0)
    for attempt in range(3):
        res_p = res_at(name, v0 + h)
        maxdev = float(np.abs(res_p - res0).max())
        if maxdev < 0.02 and attempt < 2:
            print('%s: h=%g maxdev=%.3g too small -> x30' % (name, h, maxdev), flush=True)
            h *= 30.0
            continue
        if maxdev > 30.0 and attempt < 2:
            print('%s: h=%g maxdev=%.3g too large -> /30' % (name, h, maxdev), flush=True)
            h /= 30.0
            continue
        break
    res_m = res_at(name, v0 - h)
    dcol = (res_p - res_m) / (2.0*h)
    np.savez('jac_sep/col_R5_%s.npz' % name, dcol=dcol, h_abs=h, maxdev_us=maxdev)
    report[name] = {'h_abs': h, 'maxdev_us': maxdev, 'base_value': v0}
    print('%-10s h=%-10g maxdev=%8.3f us' % (name, h, maxdev), flush=True)

with open('jac_sep/sep_R5_report.json', 'w') as fh:
    json.dump(report, fh, indent=1)
print('SEPR5_DONE')
