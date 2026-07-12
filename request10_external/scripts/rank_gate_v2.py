#!/usr/bin/env python3
"""Re-adjudicate the 10.7a Stage-2 gates with the corrected Jacobian (v2).

Recomputes, exactly as analysis_stage23.py did but on finite_jacobian_v2.npy:
  - principal cosines between the nuisance span and the 6D carrier space,
  - effective projection rank at the pre-registered cosine thresholds,
  - per-carrier-coordinate absorption,
  - dynamic-chi test-column span membership (survival) per tau.

Gate rules (from the 10.7a packet): collapse if effective rank = 6/6, or if
the dynamic-chi column lies in the finite nuisance span.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

base = np.load('baseline_planetGR.npz', allow_pickle=True)
errs = base['errs'].astype(np.float64)
t = base['toas'].astype(np.float64)
N = errs.size
J = np.load('finite_jacobian_v2.npy')
rank_v1 = json.load(open('carrier_projection_rank.json'))
OMS = [rank_v1['carriers_rad_per_day'][k] for k in ('Omega_in', 'Omega_out', 'Omega_dif')]

sw = 1.0 / errs
A = sw[:, None] * J
colnorm = np.linalg.norm(A, axis=0)
An = A / colnorm[None, :]
S = np.column_stack(sum([[np.cos(w*t), np.sin(w*t)] for w in OMS], []))
B = sw[:, None] * S

Uj, sj, _ = np.linalg.svd(An, full_matrices=False)
rel = sj / sj[0]
rank_J = int(np.sum(rel > 1e-10))
QJ = Uj[:, :rank_J]
rank_sens = {('1e-%d' % k): int(np.sum(rel > 10.0**(-k))) for k in (6, 8, 10, 12)}

QS, _ = np.linalg.qr(B)
cosines = np.linalg.svd(QJ.T @ QS, compute_uv=False)
absorb = []
for c in range(6):
    b = B[:, c] / np.linalg.norm(B[:, c])
    absorb.append(float(np.linalg.norm(QJ.T @ b)))
th_list = [0.9, 0.99, 0.999, 0.9999]
rank_eff = {str(th): int(np.sum(cosines >= th)) for th in th_list}

TAUS = [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 26.0, 52.0, 100.0, 200.0, 327.0]
span = {}
for tau in TAUS:
    s = np.zeros(N)
    for w in OMS:
        a = 1.0/(1.0 + 1j*w*tau)
        s += a.real*np.cos(w*t) - a.imag*np.sin(w*t)
    b = sw*s; b = b/np.linalg.norm(b)
    ab = float(np.linalg.norm(QJ.T @ b))
    span[str(tau)] = {'absorption': ab, 'survival': float(np.sqrt(max(0.0, 1-ab*ab)))}

out = {
    'jacobian': 'finite_jacobian_v2.npy',
    'principal_cosines_spanJ_vs_carrier6D': [float(x) for x in cosines],
    'principal_cosines_v1': rank_v1['principal_cosines_spanJ_vs_carrier6D'],
    'per_coordinate_absorption': dict(zip(
        ['cos_in', 'sin_in', 'cos_out', 'sin_out', 'cos_dif', 'sin_dif'], absorb)),
    'effective_rank_at_cosine_threshold': rank_eff,
    'jacobian_numerical_rank': rank_J,
    'rank_at_relative_sv_cutoff': rank_sens,
    'dynamic_chi_span_test': span,
    'gate_rule': 'collapse if effective rank = 6/6 or chi column inside nuisance span',
}
with open('carrier_projection_rank_v2.json', 'w') as fh:
    json.dump(out, fh, indent=1)

print('v2 principal cosines:', np.round(cosines, 6))
print('v1 principal cosines:', np.round(rank_v1['principal_cosines_spanJ_vs_carrier6D'], 6))
print('v2 rank_eff:', rank_eff, ' numerical rank %d/28' % rank_J)
print('v2 chi-column survival: ' + ', '.join(
    'tau=%s: %.3f' % (k, v['survival']) for k, v in span.items()))
print('RANK_GATE_V2_DONE')
