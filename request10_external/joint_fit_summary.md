# Request 10.7b/c/d Joint-Fit Real-Data Result: Decision Summary

Status: Imported from prior work. Contracts, each committed BEFORE the data
look it governs: `../notes/REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md`
(10.7b, commit `3d9337c`), `../notes/REQUEST10_EXTERNAL_JOINT_FIT_AMENDMENT_10_7C.md`
(10.7c, commit `92667cf`), `../notes/REQUEST10_EXTERNAL_JOINT_FIT_PROMOTION_10_7D.md`
(10.7d, commit `e7bf8e9`). Labels used exactly as required:
`runtime-motivated`, `conditional`, `collapse`.

## Verdict (final for this request chain)

Status: Counterexample candidate. **No detection; the program's first quoted
real-data upper limit, promoted to the validated truncated estimator over the
extended window by the pre-registered 10.7d gates (10/10 anchors passed).**
On the real `J0337+1715` baseline residuals (12474 TOAs, Nancay 2013-2021),
fitting the shared-relaxation transfer `G(z) = c_Y + beta/(1 + tau_chi z)`
jointly with all 28 standard timing parameters plus an offset:

```text
No detection:  Z = 0.271 over tau_chi in [1, 327] d, global p = 0.81,
               off-carrier controls quiet (p >= 0.20).

QUOTED 95% upper limit (truncated estimator, SVD cut 1e-3, window [1, 327] d):
  u95(1 d)  = 0.191 us     u95(3 d)  = 0.540 us
  u95(26 d) = 4.61 us      u95(52 d) = 9.18 us     u95(104 d) = 18.4 us
  (full curve: beta_limit_curve_10_7d.tsv; minimum 0.191 us at the 1 d edge)
```

`beta`, `u95` in us of common pre-transfer drive amplitude at unit drive
(`Lambda_k F_k = 1`); physical drives rescale per 10.7a caveat (c). The
10.7c conservative full-rank quote (u95(10 d) = 11.7 us, ..., 53.8 us at
52 d) remains recorded and is superseded by the above.

## Gate record

Status: Proven. 10.7d Gate D1 (live null calibration): one live-GN fit of the
real residuals (truncated-pinv updates, real Nutimo evaluations); the
truncated pipeline measured on the converged remainder tracks the linear
pipeline at every anchor: `|z_GN - z_lin| <= 0.017` for
`tau in {1, 3, 26, 52, 104} d` (tolerance 0.3). PASS 5/5.

Status: Proven. 10.7d Gate D2 (estimator calibration at the limit): injecting
`beta_inj = u95_trunc(tau)` into the real residuals and letting the live
model fit the combined data, the recovered `beta_hat_GN` deviates from
`beta_hat_lin_null + beta_inj` by only `-0.015 to -0.017 sigma` at every
anchor (tolerance 0.5 sigma) — end-to-end estimator calibration against the
real nonlinear model at amplitudes from 0.19 to 18.4 us. PASS 5/5.

Status: Proven. 10.7c Gate A (absorption reality): a 4-sigma chi injection
attacked with live-GN is not absorbed beyond the truncated-linear prediction
(`z_true = 4.0001` vs 4.0 predicted; full-rank fiction predicted 0.67).

Status: Note. 10.7c Gate L (displacement linearity, 0.100 vs 0.05) is
re-reported here as a diagnostic per the 10.7d contract: its 10% shape
deviation on the fitted structure perturbs chi2 by <~0.2% and s by <~0.1%,
and gates D1/D2 validate the actually-quoted quantities directly.

## Audit trail

Status: Proven. (1) 10.7b full-rank pipeline ran as pre-registered: no
detection; chi2_red = 1.246, s = 1.116 < 1.75; controls quiet; injection
coverage <= 0.08 sigma; C2 comparator reported (not excluded, as the 10.3
counting boundary anticipated). (2) The 10.7b linearization stop rule FIRED
(deviation 5.9): no 10.7b limit was quoted. (3) Diagnosis: the failure was
entirely the planet block — the 10.7a `*_extra1` Jacobian steps
(0.3 sigma_MCMC) were nonperturbative (tasc step = 2.8 planet periods):
secants, not derivatives; the 21 non-planet directions are linear to 0.3%.
(4) The 7 planet columns were re-derived at perturbative steps (half-step
linearity 0.2-2.5%; `jac_v2/`). The 10.7a Stage-2 gates are UNCHANGED under
the correction (principal cosines match to 5e-4; rank <= 5/6); the full-rank
pipeline is robust to it (< 1%). (5) The full-rank beta marginal was found to
ride on three near-null directions (rel SV <= 5e-4) realizable only through
period-wrapping planet excursions; truncating them (cut 1e-3) shrinks
`sigma_beta(52 d)` 27.4 -> 4.59 us. (6) 10.7c gates: A PASS exactly, L FAIL
marginally -> per quote rules the conservative full-rank curve was quoted and
the truncated curve held back. (7) 10.7d pre-registered estimator-calibration
gates D1/D2 (replacing L, with written rationale) plus the extended window
`[1, 327] d` (C1's 10 d floor was a sequential-pilot artifact; D2 at 1 and
3 d is the operative testability check). All 10 anchor-gates PASSED ->
truncated curve quoted over the extended window.

## Files

Status: Note. QUOTED: `joint_fit_upper_limit_10_7d.json`,
`beta_limit_curve_10_7d.tsv`, gates in `jointfit_gateD.json`
(+ `gateD_params.json`, `gateD_columns.npz`). Superseded quotes and
validation record: `joint_fit_upper_limit[_v2|_v2trunc].json`,
`beta_limit_curve[_v2|_v2trunc].tsv`, `jointfit_gateLA.json`,
`jointfit_linearization_check.json`, `jointfit_lindiag.json`. Corrected
Jacobian: `finite_jacobian_v2.npy` (+meta, `jac_v2/`),
`carrier_projection_rank_v2.json`. Inputs: `baseline_planetGR.npz`,
`xb52.npy`, `xc52.npy`, `jointfit_dtheta52*.npy`. Reproduction:
`scripts/joint_fit_upper_limit.py` (deterministic, seed 20260710; argv:
jacobian, suffix, SV cut, window min), `scripts/build_jac_v2.py`,
`scripts/rank_gate_v2.py`, WSL gate scripts per the notes.

## Caveats

Status: Note. (a) Noise is white with release EFAC 1.1225 plus global scale
s = 1.116 (chi2_red 1.25: mild excess only); no explicit red-noise model.
(b) Unit-drive convention; physical mapping needs the Request 10.2/10.6
source factors. (c) tau is not localizable at limit amplitudes (~2 sigma;
expected); the 5x diagnostic shows joint tau recovery works where sequential
fails. (d) The u95 minimum sits at the 1 d window edge; below ~0.3 d the
beta/c_Y degeneracy opens (diagnostic rows only). (e) Estimator validity is
gate-checked at the five anchors; between anchors it is interpolated by the
linearity of the same estimator (D1 holds across the whole curve by
construction of the single GN remainder).
