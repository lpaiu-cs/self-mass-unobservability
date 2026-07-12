#!/usr/bin/env python3
"""Request 10.8d turn-search robustness (10.8f-corrected re-run).

RULES OF RECORD (the amended criteria; the original pre-registered
wrapped-ramp-match and 25-chi2-margin rules were superseded by committed
amendments BEFORE the stages they governed -- see
notes/REQUEST10_8D_TURN_SEARCH_ROBUSTNESS.md):

Stage V: PASS iff the live alias-shift response proves turns are FIXED
(maxdev under a one-turn alias exceeds P/2). Recorded from the committed
turnV artifacts.

Stage L: turn-aliased candidate solutions (integer-turn sawtooths from
frequency/frequency-derivative aliases, wrapped with round-half-away
tie-breaking) are scanned over the lattice; cells within 25 chi2 of the
best solution are 'viable' (gap-hidden slips are EXPECTED to be chi2-free);
PASS iff the estimator's beta-hat is stable across ALL viable cells:
|beta_hat(cell) - beta_inj| <= max(0.1*beta_inj, 3*sigma_Fisher), and the
(0,0) recovery satisfies the same max() tolerance.

10.8f corrections: C1 causal templates (sep_common), C4 guard kept,
C8 tie-broken wrap + WIDENED lattice (n1 in [-12,12], n2 in [-6,6],
16 phase offsets, 3 quadratic references) + max() tolerance at both sites,
C10 amplitudes from the corrected 10.8b artifact.
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

inp = sc.load_inputs(HERE)
sw, res0, t, N = inp['sw'], inp['res0'], inp['t'], inp['N']
F_SPIN = 3.1618368394753872748e7
P_US = 86400e6/F_SPIN
tau0 = t.min()
T_SPAN = float(t.max() - t.min())

# ---- Stage V (amended; recorded from committed live measurements) ----
vc = np.load('sep_dynamic/turnV_columns.npz')
V = {}
for key, n in (('dres_n0.5', 0.5), ('dres_n1', 1.0)):
    meas = vc[key].astype(np.float64)
    V['n%g' % n] = {'maxdev_us': float(np.abs(meas).max()),
                    'exceeds_P_over_2': bool(np.abs(meas).max() > P_US/2)}
v_pass = V['n1']['exceeds_P_over_2']
out = {'stage_V': V, 'stage_V_pass': bool(v_pass),
       'rules_of_record': 'amended criteria; see module docstring and the 10.8d note'}
print('Stage V: turn-fixing %s (maxdev %.0f us vs P/2 = %.0f us)'
      % ('PROVEN' if v_pass else 'NOT ESTABLISHED', V['n1']['maxdev_us'], P_US/2))
if not v_pass:
    json.dump(out, open('sep_dynamic/turn_search_10_8d.json', 'w'), indent=1)
    raise SystemExit('STOP RULE: turn-fixing not established')

# ---- Stage L (causal templates, tie-broken wrap, widened lattice) ----
nu = sc.build_nuisance(inp)
Qb = nu['Q']
S_SCALE, _ = sc.noise_scale(Qb, sw, res0)
fit10_8b = json.load(open('sep_dynamic/sep_joint_fit_10_8b.json'))

N1 = range(-12, 13)
N2 = range(-6, 7)
PHI0 = np.linspace(0, 1, 16, endpoint=False)
TREF2 = (tau0, tau0 + 0.25*T_SPAN, tau0 + 0.5*T_SPAN)
CASES = [('fisher_2d', 2.0, fit10_8b['anchors']['tau_2']['u95_fisher']),
         ('Kdyn_2d', 2.0, fit10_8b['anchors']['tau_2']['u95_Kdyn']),
         ('Kstatic_18d', 18.0, fit10_8b['anchors']['tau_18']['u95_Kstatic'])]

out['stage_L'] = {}
l_pass = True
for label, tau, b_inj in CASES:
    Tc, Tb = sc.templates(inp, tau)
    SPAN = np.column_stack([Qb, sw*Tc/np.linalg.norm(sw*Tc),
                            sw*Tb/np.linalg.norm(sw*Tb)])
    Uq, sq, _ = np.linalg.svd(SPAN, full_matrices=False)
    Q = Uq[:, :int(np.sum(sq/sq[0] > 1e-12))]
    d = sw*(res0 + b_inj*Tb)
    d_perp = d - Q @ (Q.T @ d)
    chi2_0 = float(d_perp @ d_perp)
    Xs = np.column_stack([sw*Tc, sw*Tb])
    X2 = sc.proj_out(Qb, Xs)
    Minv2 = np.linalg.inv(X2.T @ X2)
    th0 = Minv2 @ (X2.T @ (d - Qb @ (Qb.T @ d)))
    b_rec = float(th0[1])
    sigF = S_SCALE*math.sqrt(Minv2[1, 1])
    tol = max(0.1*b_inj, 3.0*sigF)
    rec_ok = abs(b_rec - b_inj) <= tol
    min_margin = None
    worst = None
    n_viable = 0
    for n1 in N1:
        for n2 in N2:
            if n1 == 0 and n2 == 0:
                continue
            best_cell = None
            for tref2 in TREF2:
                for phi0 in PHI0:
                    ramp = (n1*P_US*(t - tau0)/T_SPAN
                            + n2*P_US*((t - tref2)/T_SPAN)**2 + phi0*P_US)
                    s = sc.wrap_half_away(ramp, P_US)
                    ds = d - sw*s
                    dsp = ds - Q @ (Q.T @ ds)
                    c2 = float(dsp @ dsp)
                    if best_cell is None or c2 < best_cell[0]:
                        best_cell = (c2, ds)
            margin = best_cell[0] - chi2_0
            if min_margin is None or margin < min_margin:
                min_margin, argmin = margin, (n1, n2)
            if margin < 25.0:
                n_viable += 1
                dsb = best_cell[1]
                thc = Minv2 @ (X2.T @ (dsb - Qb @ (Qb.T @ dsb)))
                dev = abs(float(thc[1]) - b_inj)
                if worst is None or dev > worst[0]:
                    worst = (dev, float(thc[1]), n1, n2, margin)
    ok = (worst is None) or (worst[0] <= tol)
    l_pass = l_pass and ok and rec_ok
    out['stage_L'][label] = {
        'tau': tau, 'beta_inj': b_inj, 'chi2_0': chi2_0,
        'beta_recovered_00': b_rec, 'recovery_tol': tol,
        'recovery_ok': bool(rec_ok),
        'lattice': {'n1': [-12, 12], 'n2': [-6, 6], 'n_phi0': len(PHI0),
                    'n_tref2': len(TREF2)},
        'min_lattice_margin': min_margin, 'argmin_lattice': list(argmin),
        'n_chi2_viable_cells': n_viable,
        'beta_stability_tol': tol,
        'worst_cell': (None if worst is None else
                       {'abs_dev': worst[0], 'beta_hat': worst[1],
                        'n1': worst[2], 'n2': worst[3], 'margin': worst[4]}),
        'pass': bool(ok)}
    print('Stage L %s: b_inj=%.4g rec dev=%.3g (tol %.3g) | viable %d | worst %s -> %s'
          % (label, b_inj, abs(b_rec - b_inj), tol, n_viable,
             'none' if worst is None else '%.3g at (%d,%d)' % (worst[0], worst[2], worst[3]),
             'PASS' if (ok and rec_ok) else 'FAIL'), flush=True)

out['stage_L_pass'] = bool(l_pass)
out['verdict'] = 'PASS' if (v_pass and l_pass) else 'FAIL'
out['scope_note'] = ('robustness is established for the TESTED lattice '
                     '(|n1|<=12, |n2|<=6, 48 phase/reference combos per cell), '
                     'not for arbitrary turn re-assignments')
json.dump(out, open('sep_dynamic/turn_search_10_8d.json', 'w'), indent=1)
print('10_8D_%s' % out['verdict'])
