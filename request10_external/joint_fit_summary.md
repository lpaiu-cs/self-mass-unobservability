# Request 10.7b-e / 10.8b Joint-Fit Real-Data Results: Decision Summary

Status: Imported from prior work. Contracts, each committed BEFORE the data
look it governs: `../notes/REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md`
(10.7b, commit `3d9337c`), `../notes/REQUEST10_EXTERNAL_JOINT_FIT_AMENDMENT_10_7C.md`
(10.7c, commit `92667cf`), `../notes/REQUEST10_EXTERNAL_JOINT_FIT_PROMOTION_10_7D.md`
(10.7d, commit `e7bf8e9`), `../notes/REQUEST10_EXTERNAL_PHYSICAL_ANCHOR_REDNOISE_10_7E.md`
(10.7e, commit `e7f5f5a`). Labels used exactly as required:
`runtime-motivated`, `conditional`, `collapse`.

## Verdict (final for this request chain)

Status: Counterexample candidate. **No detection in any pipeline mode
(global p = 0.50-0.95, off-carrier controls quiet). The program's first
quoted real-data upper limit, promoted by the 10.7d gates (10/10 anchors) and
finalized by the pre-registered 10.7e red-noise rule.** On the real
`J0337+1715` baseline residuals (12474 TOAs, Nancay 2013-2021), fitting the
shared-relaxation transfer `G(z) = c_Y + beta/(1 + tau_chi z)` jointly with
all 28 standard timing parameters, an offset, and (headline mode) 30
low-frequency Fourier pairs:

```text
QUOTED 95% upper limit (truncated estimator, red-noise marginalized,
window [1, 327] d; unit drive):
  u95(1 d)  = 0.202 us     u95(3 d)  = 0.573 us
  u95(26 d) = 4.95 us      u95(52 d) = 9.98 us     u95(104 d) = 20.2 us
  (curve: beta_limit_curve_rn.tsv; minimum 0.202 us at the 1 d edge)

E1 — implied per-carrier pole-amplitude bounds (headline rn curve):
  inner and difference carriers:  A_pole,95 ~= 49-51 ns for ALL tau in [1, 104] d
  outer carrier:                  A_pole,95 = 0.20 us (1 d) ... 9.0 us (104 d)

E2 — dictionary-anchored physical bound (Conjectural, model-anchored,
Request-8 leading-order drive D_k = (0.96, 545, 0.25) us, O(1) geometry):
  dimensionless dynamic sensitivity-response amplitude
  beta_phys < 0.36-0.59 (95%) across tau_chi in [1, 104] d
  (red-marginalized; minimum 0.36 at tau = 26 d) — an order-unity lagged
  response of the effective sensitivity to the companion driving potential
  is excluded.
```

The white-noise 10.7d curve (min u95 = 0.191 us) is superseded per the
pre-registered 10% rule — the maximum red/white ratio on the quoted window
reached 1.1001 (at tau = 297 d), and the projected-residual periodogram
confirms real low-frequency excess (median amplitude ratio 12.6 between
f < 10/T and the mid-band). The 10.7c conservative full-rank quote
(u95(52 d) = 53.8 us) remains recorded; both are superseded by the above.

## Request 10.8b: dynamic free-fall SEP channel (HEADLINE)

Status: Counterexample candidate. **No detection; first upper limit on a
dynamically-responding SEP violation.** Contract:
`../notes/REQUEST10_8B_REALDATA_DYNAMIC_SEP_FIT.md` (pre-registered commit
`2282924`; G2 amendment `6e55569` committed before the gate rerun). Signal:
the shared one-state pole transfer applied to the SEP parameter `Delta` in
the free-fall sector, templates from the T2 MEASURED integrator response
columns (patched build, A=0 gate bit-exact; time convention resolved by the
quasi-static ramp with the integral-kernel 1/2 signature).

```text
No detection: Z = 0.318 over tau_chi in [2, 500] d, global p = 0.95;
              anti-causal (advanced-response) control quiet (p = 0.48).

QUOTED 95% upper limit (anchored estimator, K = 934; window [2, 500] d):
  beta_ff < 1.09e-7 (tau = 2 d, curve minimum)
  anchors: 1.09 / 1.20 / 1.31 / 1.73 / 4.73 x 1e-7 at tau = 2/5/18/52/200 d
  (Fisher statistical floor: 1.2e-10, reported alongside; curve in
   sep_dynamic/sep_beta_limit_curve_10_8b.tsv)
```

`beta_ff` = dimensionless amplitude of the `Delta`-oscillation pole residue
at unit drive per carrier. Compared with the same data's clock-sector bound
(`beta_phys < 0.36-0.59`, 10.7e), the free-fall channel is **~3 x 10^6
deeper**, and it probes a coupling the published static SEP analyses cannot
see by construction (an oscillating `Delta` is 99.94% -> 11-28% LESS
absorbable than a static one; F2).

Status: Proven (10.8b gates, 6/6). G1 live null: `|z_GN - z_lin| <= 0.019`
at both anchors (tol 0.3). G2a detection-regime (4 sigma_F injections,
undamped live GN): recovery within 0.06 sigma_F (tol 0.5). G2b
limit-amplitude (u95_anchored injections ~ 48,000 whitened units, damped GN
capped at 30x validated steps): the live model absorbs essentially nothing
(|w r| unchanged) and recovery holds to 0.16-0.23% relative (tol 2%). The
original G2 (undamped at 1830 sigma_F) segfaulted the integrator and was
amended BEFORE rerun (commit `6e55569`) — recorded, not hidden.

Status: Counterexample candidate (10.8c re-anchor, executed per the
pre-registered rule of `../notes/REQUEST10_8C_ANCHOR_RESOLUTION.md`, commit
`26abfa1`). The ratio spectrum `rho_j = sigma_MCMC/sigma_Fisher` over the 20
core parameters REJECTS the global-optimism hypothesis: `rho <= 1`
everywhere (median ~0; 0.8-0.94 on convention-clean directions, i.e. Fisher
calibrated-to-conservative), so the static 934x gap is direction-specific —
the turn-wrap manifold seen directly in F0 (onset at `|Delta| ~ 1e-7`) sets
the published static width and was proven inoperative on the dynamic
templates (G2b). Rule branch 1 fires: `K_dyn = max(10, ceil(max rho)) = 10`:

```text
RE-ANCHORED 95% upper limit (K_dyn = 10, window [2, 500] d):
  beta_ff < 1.17e-9 (tau = 2 d);  1.29 / 1.40 / 1.86 / 5.07 x 1e-9
  at tau = 5/18/52/200 d   (curve: sep_beta_limit_curve_10_8c.tsv;
  Fisher floor and the K = 934 ultra-conservative bracket recorded alongside)
```

Status: Note. Caveats: (a) ONE absorber remains unprobed — turn
re-assignment (our GN fixes turns at construction; the published pipeline
re-searches them): flagged as 10.8d, REQUIRED before publication-grade
claims; (b) drive normalization is unit-drive per carrier (physical source
factors as in 10.7e E2 pending); (c) zero-drive-phase convention (a
phase-marginalized two-quadrature variant is a mechanical extension); (d)
tau below 2 d diagnostic only; (e) white + 30-pair Fourier noise model as in
the 10.7e headline.

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

Status: Proven. (8) 10.7e ran the pre-registered physical-anchor and
red-noise modes (all deterministic, same seed; no detection in any mode):
the red/white u95 ratio reached 1.1001 > 1.10 at tau = 297 d, so per the
pre-registered rule the red-marginalized curve became the headline quote
(anchor-level inflation is only 6-10%); the projected-residual periodogram
(`residual_periodogram.json`) shows a 12.6x low-frequency amplitude excess,
independently justifying the marginalization. Gate carry-over per the 10.7e
note: the added Fourier columns and the dictionary re-weighting are exactly
linear; the D1/D2-validated 28-parameter model content is unchanged.
E1 arithmetic in `pole_amplitude_bounds.json` (both the letter-of-contract
10.7d version and the headline-consistent rn version). E2 template:
`D_k = (0.9565, 544.9, 0.2526) us`, phases locked to `tasc_p/tasc_b` with the
10.5 phase-lock for the combination carrier.

## Files

Status: Note. QUOTED (headline): `joint_fit_upper_limit_rn.json`,
`beta_limit_curve_rn.tsv`; physical anchor `joint_fit_upper_limit_phys_rn.json`
(+ `_phys` white cross-check); E1 `pole_amplitude_bounds.json`; periodogram
`residual_periodogram.json`. Promotion record: `joint_fit_upper_limit_10_7d.json`,
`beta_limit_curve_10_7d.tsv` (white-noise version, superseded by the 10% rule),
gates in `jointfit_gateD.json` (+ `gateD_params.json`, `gateD_columns.npz`).
Superseded quotes and validation record:
`joint_fit_upper_limit[_v2|_v2trunc].json`, `beta_limit_curve[_v2|_v2trunc].tsv`,
`jointfit_gateLA.json`, `jointfit_linearization_check.json`,
`jointfit_lindiag.json`. Corrected Jacobian: `finite_jacobian_v2.npy`
(+meta, `jac_v2/`), `carrier_projection_rank_v2.json`. Inputs:
`baseline_planetGR.npz`, `xb52.npy`, `xc52.npy`, `jointfit_dtheta52*.npy`.
Reproduction: `scripts/joint_fit_upper_limit.py` (deterministic, seed
20260710; argv: jacobian, suffix, SV cut, window min, rn-pairs, template),
`scripts/build_jac_v2.py`, `scripts/rank_gate_v2.py`, WSL gate scripts per
the notes.

## Caveats

Status: Note. (a) The headline noise model is white (release EFAC 1.1225 +
global scale) plus 30 unconstrained low-frequency Fourier pairs — a
prior-free red-noise marginalization, not a spectral-model fit. (b) The E2
physical anchor sets all O(1) geometric projection factors to 1 and adopts
the Request-8 leading-order clock drive (`Conjectural` per the ledger); its
numbers scale inversely with the true D_k. (c) tau is not localizable at
limit amplitudes (~2 sigma; expected); the 5x diagnostic shows joint tau
recovery works where sequential fails. (d) The unit-drive u95 minimum sits at
the 1 d window edge; below ~0.3 d the beta/c_Y degeneracy opens (diagnostic
rows only). (e) Estimator validity is gate-checked at the five anchors under
the white-noise span; the 10.7e extensions add exactly-linear columns only
(carry-over argument in the 10.7e note).
