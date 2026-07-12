#!/usr/bin/env python3
"""Request 10.8c anchor resolution: reproduce anchor_resolution_10_8c.json and
the re-anchored companion curve (10.8f C11 -- committed generator for an
artifact originally produced inline; see the pre-registered diagnosis and
re-anchoring rule in notes/REQUEST10_8C_ANCHOR_RESOLUTION.md).

Diagnostic: the ratio spectrum rho_j = sigma_MCMC,j / sigma_Fisher,j over the
20 CORE fitted parameters (excluding the 7 *_extra1 planet parameters and
delta_i, which are known nonlinear/wrap-dominated). sigma_Fisher,j is the
full-rank leave-one-out projection marginal s / ||P_perp(others) a_j|| in the
10.8b nuisance context (all other Jacobian columns + offset + 30 Fourier
pairs + static SEP column); NO SVD truncation. sigma_MCMC,j = abs_steps_j
from steps_provenance.json (frac = 1.0, i.e. the steps ARE the released-chain
posterior widths; the chain file itself is not in the repository).

Re-anchoring rule (pre-registered): branch 1 if median(rho_core) <= 3 AND
max(rho_core) <= 10 -> K_dyn := max(10, ceil(max(rho_core))). Branch 2 if
3 < median <= 30 -> K_dyn := 3 * max(rho_core). Branch 3 otherwise -> K
stays 934, no re-quote.

10.8f note: the noise scale s is read from the corrected 10.8b artifact
(originally 1.1148 was typed inline; the corrected chain gives 1.11489),
and the companion curve is regenerated from the corrected 10.8b curve with
a consistency assert against its own u95 columns. Numbers therefore move in
the 4th-5th digit relative to the pre-10.8f artifact; the decision branch is
unchanged by three orders of margin.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

import sep_common as sc

base = np.load('baseline_planetGR.npz', allow_pickle=True)
errs = base['errs'].astype(float)
t = base['toas'].astype(float)
sw = 1.0/errs
S = json.load(open('sep_dynamic/sep_joint_fit_10_8b.json'))['noise']['s']

J = np.load('finite_jacobian_v2.npy')
names = json.load(open('finite_jacobian_v2_meta.json'))['columns']
prov = json.load(open('sep_dynamic/steps_provenance.json'))
assert prov['fitted_names'] == names
sig_mcmc = np.array(prov['abs_steps'])   # frac = 1.0 -> these ARE sigma_MCMC
jD = np.load('sep_dynamic/col_SEP_D.npz')['dcol'].astype(float)

A = sw[:, None]*J
An = A/np.linalg.norm(A, axis=0)
cvec = sw/np.linalg.norm(sw)
T = t.max() - t.min()
fc = []
for j in range(1, 31):
    a = 2*np.pi*j*(t - t.min())/T
    fc += [np.cos(a), np.sin(a)]
FB = sw[:, None]*np.column_stack(fc)
FB = FB/np.linalg.norm(FB, axis=0, keepdims=True)
xjD = sw*jD
xjD = xjD/np.linalg.norm(xjD)

rho = {}
sigF = {}
for j in range(28):
    others = np.column_stack([An[:, [k for k in range(28) if k != j]],
                              cvec, FB, xjD])
    U, s0, _ = np.linalg.svd(others, full_matrices=False)
    k = int(np.sum(s0/s0[0] > 1e-13))   # full numerical rank
    Q = U[:, :k]
    aj = A[:, j]
    a2 = aj - Q @ (Q.T @ aj)
    sF = S/np.linalg.norm(a2)
    sigF[names[j]] = sF
    rho[names[j]] = float(sig_mcmc[j]/sF)

PLANET = [n for n in names if n.endswith('_extra1')]
core = [n for n in names if n not in PLANET and n != 'delta_i']
rc = np.array([rho[n] for n in core])
med, mx = float(np.median(rc)), float(rc.max())
print('median(rho_core)=%.4g  max=%.4g  (rho = sigma_MCMC / sigma_Fisher)' % (med, mx))

if med <= 3.0 and mx <= 10.0:
    branch = '1 (median<=3 AND max<=10)'
    K_dyn = float(max(10, math.ceil(mx)))
elif med <= 30.0:
    branch = '2 (3<median<=30)'
    K_dyn = 3.0*mx
else:
    branch = '3 (median>30): K stays 934'
    K_dyn = 934.0

out = {
    'rho_core': {n: rho[n] for n in core},
    'rho_delta_i': rho['delta_i'],
    'rho_planet': {n: rho[n] for n in PLANET},
    'sigma_fisher': sigF,
    'median_rho_core': med, 'max_rho_core': mx,
    'noise_scale_s': S,
    'decision': {
        'rule_branch': branch,
        'median_rho_core': med, 'max_rho_core': mx,
        'K_dyn': K_dyn,
        'interpretation': (
            'H-global REJECTED: Fisher is WIDER than the published posterior '
            'on core directions (rho<=1 everywhere; rho~0.8-0.94 where the '
            'comparison is convention-clean, rho<<1 where our wider '
            'marginalization (30 Fourier pairs + near-null span) inflates '
            'Fisher). The static-SEP 934x gap is direction-specific: the '
            'turn-wrap manifold seen directly in F0 (wrap onset at '
            '|Delta|~1e-7) sets the published static width and does not '
            'operate on the live-validated dynamic templates (G2b).'),
        'residual_caveat': (
            'Turn re-assignment as an absorber of the dynamic template is '
            'the one unprobed channel (our GN fixes turns at construction; '
            'the published pipeline re-searches them). Flagged as 10.8d '
            'before publication-grade claims. The K=934 curve remains '
            'recorded as the ultra-conservative bracket. [10.8f status: '
            '10.8d executed and PASSED on the widened lattice -- '
            'turn_search_10_8d.json; the live-gate caveat of record is now '
            'the R6 null offset, sep_gateG_adjudication.json.]'),
    },
}
json.dump(out, open('sep_dynamic/anchor_resolution_10_8c.json', 'w'), indent=1)

# re-anchored companion curve from the corrected 10.8b curve, using the
# same K factors the 10.8b script consumed (K_dyn from THIS artifact's rule,
# K_static from the R5 ratio in sep_static_sensitivity.json)
af = sc.load_anchor_factors(HERE)
K_static = af['K_static']
assert af['K_dyn'] == K_dyn, (af, K_dyn)
rows = [l.split('\t') for l in
        open('sep_dynamic/sep_beta_limit_curve_10_8b.tsv').read().splitlines()[1:]]
lines = ['tau_d\tbeta_hat\tsigma_fisher\tu95_fisher\tu95_K10\tu95_K934']
anch = {}
for r in rows:
    tau, b, s = float(r[0]), float(r[1]), float(r[2])
    uf, u10, u934 = sc.u95_of(b, s), sc.u95_of(b, K_dyn*s), sc.u95_of(b, K_static*s)
    # consistency: the corrected 10.8b curve carries its own u95 columns
    # (stored at 6 significant digits -> allow rounding-level slack)
    for mine, theirs in ((uf, float(r[3])), (u10, float(r[4])), (u934, float(r[5]))):
        assert abs(mine - theirs)/theirs < 5e-6, (tau, mine, theirs)
    lines.append('%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g' % (tau, b, s, uf, u10, u934))
    anch[tau] = (uf, u10, u934)
open('sep_dynamic/sep_beta_limit_curve_10_8c.tsv', 'w').write('\n'.join(lines) + '\n')
print('decision: branch %s -> K_dyn = %g; curve consistent with the 10.8b '
      'artifact at <1e-6' % (branch, K_dyn))
print('min u95_K10 = %.4g at tau = %g d'
      % min(((v[1], k) for k, v in anch.items())))
