#!/usr/bin/env python3
"""gateG2 (amendment G-2): the 10.8b live gates G1/G2a/G2b re-run with the
C4-COMPLETED harness -- the measurement span keeps the static-guard
residual direction, exactly matching the pipeline estimator (verified to
reproduce z_lin on res0 to machine precision) -- and with the headline
anchor tau = 2 d added. GN absorb semantics identical to R6/10.7d.
Criteria unchanged from the R6 registration."""
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
errs = base['errs'].astype(np.float64)
t = base['toas'].astype(np.float64)
scales = base['scales']
fmap = list(base['fmap'])
sc = np.array([float(scales[p]) for p in fmap])
w = 1.0/errs
N = res0.size

J = np.load('finite_jacobian_v2.npy')
gi = np.load('gateG2_inputs.npz')
gp = json.load(open('gateG2_params.json'))
SV_CUT = gp['sv_cut']
TEMPL = {'2': gi['T2'].astype(np.float64),
         '18': gi['T18'].astype(np.float64),
         '52': gi['T52'].astype(np.float64)}
TC = gi['Tc'].astype(np.float64)
# the committed static-guard column, staged inside gateG2_inputs.npz (the
# run of record used the scratch jac_sep/ copy, verified md5-identical to
# the committed sep_dynamic/col_SEP_D.npz)
jD = gi['jD'].astype(np.float64)

# ---- 10.8b nuisance span, C4-complete: guard residual kept ----
A = w[:, None]*J
An = A/np.linalg.norm(A, axis=0)
cvec = w/np.linalg.norm(w)
T_SPAN = float(t.max() - t.min())
fc = []
for j in range(1, 31):
    arg = 2.0*np.pi*j*(t - t.min())/T_SPAN
    fc += [np.cos(arg), np.sin(arg)]
FB = w[:, None]*np.column_stack(fc)
FB = FB/np.linalg.norm(FB, axis=0, keepdims=True)
xjD = w*jD
xn = xjD/np.linalg.norm(xjD)
B0 = np.column_stack([An, cvec, FB, xn])
U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
k0 = int(np.sum(s0/s0[0] >= SV_CUT))
Qk0 = U0[:, :k0]
g = xn - Qk0 @ (Qk0.T @ xn)
Qk = np.column_stack([Qk0, g/np.linalg.norm(g)])   # C4: rank k0+1
print('measurement span rank %d (guard residual %.4g kept)'
      % (Qk.shape[1], float(np.linalg.norm(g))), flush=True)
Ua, sa, Vat = np.linalg.svd(An, full_matrices=False)
ka = int(np.sum(sa/sa[0] >= SV_CUT))


def trunc_update(r_us):
    c = Vat[:ka].T @ ((Ua[:, :ka].T @ (w*r_us))/sa[:ka])
    return c/np.linalg.norm(A, axis=0)


def measure(r_us, key):
    Xs = np.column_stack([w*TC, w*TEMPL[key]])
    X2 = Xs - Qk @ (Qk.T @ Xs)
    Minv = np.linalg.inv(X2.T @ X2)
    yr = w*r_us
    th = Minv @ (X2.T @ (yr - Qk @ (Qk.T @ yr)))
    return float(th[1])


t0 = time.perf_counter()
fit = PyFittriple('parfile-planetGR-max-bestfit',
                  '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
fit.Compute_lnposterior()
print('constructed in %.0fs' % (time.perf_counter()-t0), flush=True)


def eval_at(dth):
    fit.Set_fitted_parameter_relativeshifts(dth/sc)
    fit.Compute_lnposterior()
    return fit.Get_time_residuals().copy()


STEPS = np.array(json.load(open('artifacts/finite_jacobian_meta.json'))['abs_steps'])
try:
    STEPS_V2 = json.load(open(os.path.expanduser(
        '~/work/nutimo_pilot/run_planetGR/finite_jacobian_v2_meta.json')))['abs_steps']
    STEPS = np.array(STEPS_V2)
except Exception:
    pass


def gn_absorb(beta_inj, tpl, label, cap=None):
    dth = np.zeros(28)
    r = res0 + beta_inj*tpl
    hist = [float(np.linalg.norm(w*r))]
    for it in range(1, 5):
        dth = dth - trunc_update(r)
        if cap is not None:
            lim = cap*STEPS
            dth = np.clip(dth, -lim, lim)
        r = eval_at(dth) + beta_inj*tpl
        hist.append(float(np.linalg.norm(w*r)))
        print('  %s GN %d: |w r| = %.3f  max|dth|/step=%.1f' %
              (label, it, hist[-1], float(np.max(np.abs(dth)/STEPS))), flush=True)
        if it >= 2 and abs(hist[-2]-hist[-1])/max(hist[-2], 1e-12) < 0.005:
            break
    return r, hist


ANCHORS = ('tau_2', 'tau_18', 'tau_52')
out = {'params': gp, 'harness': 'C4-complete (guard-kept measurement span)',
       'G1': {}, 'G2a': {}, 'G2b': {}}

# G1: live null (one absorb, measured at all three anchors)
r_null, h = gn_absorb(0.0, np.zeros(N), 'null')
g1_ok = True
for a in ANCHORS:
    ap = gp['anchors'][a]
    z_gn = measure(r_null, a.split('_')[1])/ap['sigma_fisher']
    ok = abs(z_gn - ap['z_lin']) <= 0.3
    g1_ok = g1_ok and ok
    out['G1'][a] = {'z_GN': z_gn, 'z_lin': ap['z_lin'],
                    'offset': z_gn - ap['z_lin'], 'pass': bool(ok)}
    print('G1 tau=%s: z_GN=%.4f z_lin=%.4f offset=%+.4f -> %s'
          % (a, z_gn, ap['z_lin'], z_gn - ap['z_lin'],
             'PASS' if ok else 'FAIL'), flush=True)
out['G1']['pass'] = bool(g1_ok)

# G2a: detection-regime calibration (4 sigma_F, undamped)
g2a_ok = True
for a in ANCHORS:
    ap = gp['anchors'][a]
    b_inj = 4.0*ap['sigma_fisher']
    r_fin, h = gn_absorb(b_inj, TEMPL[a.split('_')[1]], 'G2a %s' % a)
    bhat = measure(r_fin, a.split('_')[1])
    dev = (bhat - ap['beta_hat_lin'] - b_inj)/ap['sigma_fisher']
    ok = abs(dev) <= 0.5
    g2a_ok = g2a_ok and ok
    out['G2a'][a] = {'beta_inj': b_inj, 'beta_hat_GN': bhat,
                     'deviation_sigma_fisher': dev, 'W_hist': h, 'pass': bool(ok)}
    print('G2a tau=%s: b_inj=%.4g dev=%+.4f sigma_F -> %s'
          % (a, b_inj, dev, 'PASS' if ok else 'FAIL'), flush=True)
out['G2a']['pass'] = bool(g2a_ok)

# G2b: limit-amplitude regime (u95_Kstatic, damped GN, 2% relative)
g2b_ok = True
for a in ANCHORS:
    ap = gp['anchors'][a]
    b_inj = ap['u95_Kstatic']
    r_fin, h = gn_absorb(b_inj, TEMPL[a.split('_')[1]], 'G2b %s' % a, cap=30.0)
    bhat = measure(r_fin, a.split('_')[1])
    dev_rel = (bhat - ap['beta_hat_lin'] - b_inj)/b_inj
    ok = abs(dev_rel) <= 0.02
    g2b_ok = g2b_ok and ok
    out['G2b'][a] = {'beta_inj': b_inj, 'beta_hat_GN': bhat,
                     'deviation_relative': dev_rel, 'W_hist': h, 'pass': bool(ok)}
    print('G2b tau=%s: b_inj=%.4g dev_rel=%+.5f -> %s'
          % (a, b_inj, dev_rel, 'PASS' if ok else 'FAIL'), flush=True)
out['G2b']['pass'] = bool(g2b_ok)

with open('sep_gateG2.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('GATEG2_DONE G1=%s G2a=%s G2b=%s'
      % (out['G1']['pass'], out['G2a']['pass'], out['G2b']['pass']), flush=True)
