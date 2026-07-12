#!/usr/bin/env python3
"""Prepare per-parameter absolute finite-difference steps (FRAC * MCMC sigma).

PROVENANCE (10.8f C11): the steps of record used FRAC = 1.0 (below); older
docstrings said 0.3 * sigma, which is wrong -- FRAC as stored in steps.npz
is authoritative. The source chain lives in the runtime tree at
Analysis/Synthesis/MCMC-planetGR.npz (not in this repository)."""
import os
import numpy as np

os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
base = np.load('baseline_planetGR.npz', allow_pickle=True)
names = list(base['names'])
fmap = list(base['fmap'])
params = base['params']

chain = np.load(os.path.expanduser(
    '~/work/nutimo_pilot/Analysis/Synthesis/MCMC-planetGR.npz'), allow_pickle=True)

FRAC = 1.0
abs_steps, srcs = [], []
for j, p in enumerate(fmap):
    nm = names[p]
    if nm in chain.files:
        sig = float(np.std(chain[nm]))
        if sig > 0:
            abs_steps.append(FRAC * sig); srcs.append('mcmc'); continue
    val = abs(float(params[p]))
    abs_steps.append(max(1e-8 * val, 1e-12)); srcs.append('fallback')

abs_steps = np.array(abs_steps)
for j, p in enumerate(fmap):
    print('%02d %-18s value=%+.6e step=%.3e (%s)' % (j, names[p], float(params[p]), abs_steps[j], srcs[j]))
np.savez('steps.npz', abs_steps=abs_steps, frac=FRAC,
         fitted_names=np.array([names[p] for p in fmap], dtype=object))
print('STEPS_SAVED nfit=%d' % len(fmap))
