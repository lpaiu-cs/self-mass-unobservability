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

[10.8f-corrected record: the first run fit the ADVANCED (anticausal)
quadrature as "causal"; an independent review caught the swap and every
data-touching stage was re-run on the corrected causal templates per the
pre-committed `../notes/REQUEST10_8F_REVIEW_RESPONSE.md` (commit
`52153b3`). Numbers below are the corrected zero-phase record.]

```text
No detection: Z = 1.328 over tau_chi in [2, 500] d, global p = 0.39;
              anti-causal (advanced-response) control quiet
              (Z = 0.598, p = 0.87).

Zero-phase curve (window [2, 500] d; 95% statistical CL x anchor tiers;
curve in sep_dynamic/sep_beta_limit_curve_10_8b.tsv):
  K_dyn = 10:  beta_ff < 1.36e-9 (tau = 2 d, curve minimum)
     anchors: 1.36 / 1.45 / 1.61 / 2.44 / 7.00 x 1e-9 at tau = 2/5/18/52/200 d
  K = 934 bracket: 1.26e-7 (2 d); Fisher-only floor: 1.93e-10 (min).
  Nuisance rank 71/90 with the static-guard residual kept explicitly
  (guard norm 1.6e-7); s = 1.1149.
```

`beta_ff` = dimensionless amplitude of the `Delta`-oscillation pole residue
at unit drive per carrier. Compared with the same data's clock-sector bound
(`beta_phys < 0.36-0.59`, 10.7e), the free-fall channel is far deeper in
raw `Delta` amplitude (~3 x 10^8 at the K10 tier), and it probes a coupling
the published static SEP analyses cannot see by construction (an
oscillating `Delta` is 99.94% -> 8-27% LESS absorbable than a static one;
F2 causal shared-tau peak 21.1% at tau = 5 d).

Status: Proven (live gates of record: gateG2, amendment G-2 of the 10.8f
note, C4-complete harness, anchors tau = 2 / 18 / 52 d INCLUDING the
headline anchor; sep_dynamic/sep_gateG2.json). **9/9 PASS under the
unchanged registered criteria**: G1 offsets -0.002 / -0.015 / -0.035
sigma_F (tol 0.3); G2a +0.009 / +0.056 / +0.057 sigma_F (tol 0.5); G2b
-0.109% / -0.024% / +0.128% (tol 2%). The live-refit systematic carried
on Fisher-tier numbers is <= 0.035 sigma_F as measured.

Audit trail (recorded, superseded): the first causal-branch gate run (R6,
sep_gateG.json) recorded G1 FAIL (+0.342/+0.599 sigma_F) and G2a FAIL at
52 d, and a first adjudication mis-attributed this to the live refit
walking to the unregularized LS point. The pre-committed amendment G-2
located the defect in the gate HARNESS: its measurement span omitted the
C4 guard-residual direction — measuring the untouched baseline residuals
(no live refit) in that span reproduces the entire offset
(+0.344/+0.607), and the guard-kept span reproduces the pipeline z_lin to
machine precision (sep_gateG_adjudication.json carries the two-layer
decomposition). The pre-C1 anticausal gate had read the same span
artifact at only 0.01-0.02 sigma_F, which is how the quadrature swap
originally passed. No criterion was rescored at any point. The original
G2 (undamped at 1830 sigma_F) segfaulted the integrator and was amended
BEFORE rerun (commit `6e55569`) — recorded, not hidden.

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
RE-ANCHORED zero-phase limit (95% statistical CL x K_dyn = 10,
window [2, 500] d; 10.8f-corrected):
  beta_ff < 1.36e-9 (tau = 2 d);  1.45 / 1.61 / 2.44 / 7.00 x 1e-9
  at tau = 5/18/52/200 d   (curve: sep_beta_limit_curve_10_8c.tsv,
  regenerated by scripts/anchor_resolution_10_8c.py; Fisher floor and the
  K = 934 ultra-conservative bracket recorded alongside)
```

Status: Proven on the tested lattice (10.8d turn-search robustness,
executed per the pre-registered-and-amended
`../notes/REQUEST10_8D_TURN_SEARCH_ROBUSTNESS.md`, commits
`349c70c/5df4870/94c04c3`; 10.8f re-run: causal templates, tie-broken
round-half-away wrap, lattice WIDENED to |n1| <= 12, |n2| <= 6 with 48
phase/reference combos per cell, tolerance max(0.1 b_inj, 3 sigma_F)).
Stage V measured that the runtime fixes turns (alias-shift maxdev
1522 us > P/2 = 1366 us), so wrap arithmetic is definitional. Stage L found
the REAL structure: 7 chi2-viable turn-alias solutions exist at every test
amplitude (integer-turn slips hidden inside observing gaps — the classic
phase-connection ambiguity, quantified here against the full span), but the
estimator's beta-hat is STABLE across all of them: worst deviation 1.01e-10
vs tolerance 2.07e-10 at the Fisher and K10 amplitudes, 8.5e-11 vs 1.5e-8
at K_static. Pulse-number re-assignment on the tested lattice cannot mimic
or hide the dynamic template (arbitrary re-assignments beyond the lattice
are not claimed). The former "unprobed absorber" caveat is RESOLVED at this
scope; the K_dyn = 10 quote stands with this record.

Status: Counterexample candidate (10.8e, pre-registered commit `3827649`,
executed). BOTH remaining review items resolved:

```text
FINAL HEADLINE (10.8f-corrected; phase-marginalized over the FULL
time-origin domain [0, P_out = 327.26 d) in 4,821 steps, worst phase per
tau and per anchor tier; window [2, 500] d; no detection under the
registered rule: E1 Z = 2.285, p = 0.26, with the anti-causal control
equally elevated Z = 2.340, p = 0.28 — elevation common to both
conventions, not causal-response-like; E2 Z = 1.815, p = 0.47):

  amplitude of any lag-responding Delta-oscillation
  |delta Delta| < 1.68e-9   (95% statistical CL x K_dyn = 10; tau_chi = 2 d)
  1.85 / 1.98 / 2.60 / 7.16 x 1e-9  at tau = 5/18/52/200 d
  (zero-phase curve 1.36e-9 recorded; Fisher-only floor 2.79e-10;
   full-rank (no-truncation) bracket 3.53e-9; K934 bracket 1.57e-7;
   window edge 1.74e-8 at 500 d; curve: sep_limit_curve_10_8e.tsv)
```

Note that `Delta` is itself the dimensionless SEP parameter, so this
amplitude bound is a physical statement independent of any drive model. The
dictionary-anchored COUPLING version (Conjectural, Request-8 `U/c^2` drive,
f_w ~ 1e-10..1e-11, tasc-locked phases): `beta_phys < 18-23` for
`tau in [2, 52] d` (59.7 at 200 d) — the free-fall channel's drive lacks
the clock channel's 1/Omega integration amplification, so in coupling units
the two channels bound different composites (clock: 0.4-0.6 on the
rate-coupling; free-fall: ~20 on the potential-coupling but ~3e8 deeper in
raw Delta amplitude).

Status: Note. Remaining caveats: (a) tau below 2 d diagnostic only;
(b) white + 30-pair Fourier noise model as in the 10.7e headline; (c) E2
coupling numbers carry O(1) geometric drive factors; (d) the live-refit
displacement measured by gate G1 with the corrected harness is
<= 0.035 sigma_F at all anchors including the headline lag (gateG2,
amendment G-2) and is carried at that size on Fisher-tier numbers;
(e) turn robustness is established for the tested alias lattice, not for
arbitrary re-assignments.

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
