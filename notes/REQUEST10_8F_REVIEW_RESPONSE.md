# Request 10.8f: Review-Response Corrections And Re-Runs

Status: Note. Two independent reviews of PR #1 (a 15-finding multi-agent
code review and an external reviewer's 5-finding report) identified
implementation defects and over-claims in the 10.8 chain. This note fixes
the correction set and the re-run plan BEFORE the re-runs are executed.
Detection rules are NOT changed (the pre-registered E1-with-control form is
restored where the implementation drifted from it).

## Accepted findings and fixes

C1 (quadrature convention swap; review #1). `col_s` is the `+pi/2` (`-sin`)
drive response, so the causal pole template is `g1 col_c - g2 col_s`.
10.8e already used this; 10.8b, the F2/F3 gates, the 10.8d lattice, and the
G-gate inputs used `+g2` (the advanced template) mislabeled as causal.
FIX: one shared implementation in `scripts/sep_common.py` with the causal
sign; re-run F2/F3, 10.8b, 10.8d, and the live G-gates on the causal branch.
The generation-phase convention is now recorded in code.

C2 (merged test suite; review #2/#3). Restore the 18 frozen-theorem tests
dropped by the union merge, correct the stale rank-4 assertions
(44/19 -> the corrected 50/25 family, matching the shipped modules and
`docs/family-class-table.md`), and fix the main() banner.

C3 (worst-phase definition; review #4, confirmed numerically). The
pre-registered `u95_pm(tau) = max over t_off of u95` is now evaluated
per-K (the Fisher-argmax shortcut understated 4/5 anchors by up to 4.5%).

C4 (static-leakage guard; review #5). The `j_Delta` guard direction fell
below the SVD cut and was inactive. FIX: the guard is kept explicitly —
the truncated basis is augmented with the guard's residual direction after
truncation — and its retained fraction is reported.

C5 (phase-scan completeness; external P2-4). `t_off in [0, P_in)` wraps the
fast carriers but moves the outer phase by only 1.8 deg. FIX: scan the full
common-origin domain `t_off in [0, P_out)` on a grid that keeps the inner
carriers' phase step below 15 deg (implemented via a 6-projected-column
factorization so the scan stays cheap). Expected effect (external reviewer's
independent computation): tau = 2 d weakens ~1.3e-9 -> ~1.33e-9 and
tau = 200 d -> ~5.3e-9.

C6 (span-truncation robustness; external P1-1). Report the limit curve at
the full-rank span alongside the truncated one as an explicit bracket, and
state the physical justification of the truncation (the live Gate A
measurement: the real integrator absorbed the injected template exactly as
the TRUNCATED model predicted, z = 4.0001 vs the full-rank prediction 0.67).
Expected full-rank bracket ~2.6x looser at tau = 2 d.

C7 (statistical framing; external P1-2). The quoted quantities are
re-labeled: the Fisher-based u95 (with live injection calibration) is the
STATISTICAL 95% CL under the stated noise model and span; K = 10 and
K = 934 are SYSTEMATIC safety factors, not probability statements. Headline
becomes "conservative bound = statistical 95% CL x safety factor", and
"validated 95% CL" phrasing is removed from K-scaled numbers.

C8 (turn-search scope and spec; review #8/#9, external P1-3). Tie-broken
wrap (round-half-away) replaces banker's rounding; the recovery tolerance
uses the registered max() form; the module docstring documents the amended
(current) rules; the lattice is WIDENED to `n1 in [-12,12]`,
`n2 in [-6,6]`, 16 phase offsets, 3 quadratic references; paper language is
softened from "any turn-aliased solution"/"proof" to the tested lattice.

C9 (detection-rule drift; review #7). The 10.8e detection logic is restored
to the pre-registered form: E1 requires global p < 0.003 AND a quiet
anti-causal control; E2 is a limit-only quote with no standalone trigger.

C10 (anchor/constant unification; review #10/#12). K is read from the
committed anchor artifact (not retyped); the noise scale is computed
per-span (never hardcoded); OMS/TASC constants are loaded from
`carrier_projection_rank.json`; the F1 white-mode sigma is recomputed with
the white-mode scale.

C11 (provenance; review #11/#13, external P2-5). The WSL-side runtime
scripts (gates, T2 columns, calibrations, patch applier) are committed under
`scripts/wsl/`; a `scripts/anchor_resolution_10_8c.py` reproduces the 10.8c
diagnosis; `sep_convention_calibration.json` gains a SUPERSEDED marker;
`prep_steps.py`/`jacobian_runner.py` docstrings are corrected to the actual
`frac = 1.0` steps convention.

C12 (paper corrections; review #6/#13/#14/#15). Window-edge claims fixed
(the curve reaches ~1.2e-8 at 500 d; exclusion depth ~2.2 orders at the
edge); the inner-carrier response figure corrected (8.5e10 us per unit, not
2.6e12); grid count 65; zero-phase digit 1.16e-9; table-width
renormalization and fence-aware abstract parsing in the builder; Makefile
target for the Paper B build.

## Re-run matrix (in order)

R1 sep_feasibility_gates (causal) -> F2/F3 re-adjudicated.
R2 sep_joint_fit_10_8b (causal, guard-kept, K-from-artifact) -> new
   gateG inputs.
R3 sep_phase_marg_10_8e (guard-kept, full t_off domain, per-K max,
   registered detection rule, cut bracket) -> NEW HEADLINE.
R4 turn_search_10_8d (causal, tie-broken wrap, widened lattice).
R5 sep_static_sensitivity (mode-correct scales) -> K artifact refreshed.
R6 live G-gates re-run in WSL on the corrected causal templates.
R7 Paper B numbers/claims regenerated; tex rebuilt.

Deferred to 10.9 (future work, stated in the paper): nonlinear likelihood /
physical priors on the truncated directions, full pulse-number
marginalization beyond the widened lattice, an independent holdout data set,
environment container/CI.

## Outcomes (recorded post-run, 2026-07-12; the matrix above was the
## pre-committed plan)

R1 PASS: F2 best survival 0.2113 at tau = 5 d; F3 min anchored 6.867e-8.
R2 done: rank 71/90 with the guard direction kept (guard residual 1.65e-7);
   s = 1.1149; Z = 1.328 (p = 0.39), anticausal control Z = 0.598 (p = 0.87);
   min u95_Kdyn = 1.365e-9 at tau = 2 d.
R3 done -> HEADLINE: phase-marginalized curve over the full t_off domain
   (4821 offsets), per-K worst-phase tracking; min u95pm_K10 = 1.68e-9 at
   tau = 2 d (Fisher-only 2.794e-10; full-rank bracket 3.534e-9; K934
   1.569e-7). E1 Z = 2.285 p = 0.26 with control equally loud; E2 Z = 1.815
   p = 0.47: no detection under the registered rules.
R4 PASS: widened lattice (|n1| <= 12, |n2| <= 6, 48 combos/cell),
   tie-broken wrap; 7 viable cells, worst beta deviation 1.01e-10 against
   tolerance 2.07e-10.
R5 done: K_static ratio refreshed = 934.44 (R5-ext sigma 1.9263e-9).
R6 recorded: **G1 FAIL (both anchors), G2a FAIL (tau = 52 d), G2b PASS.**
   The registered criteria stand; no post-hoc rescoring. Diagnosis
   (sep_gateG_adjudication.json, scripts/gateG_adjudicate.py): one
   common-mode quantity explains every excess -- the live truncated-GN
   refit of the NULL data displaces beta_hat by +0.342 / +0.599 sigma_F
   (tau = 18 / 52 d) relative to the frozen linear pipeline; subtracting
   the GN null point, the injected-signal response is accurate to
   +0.004 / +0.021 sigma_F (G2a) and +0.04% / +0.08% (G2b). The registered
   G2a/G2b deviations referenced the LINEAR null (beta_hat_lin), so the
   offset entered them whole (0.342 + 0.004 = 0.346; 0.599 + 0.021 =
   0.619 -- exact). The pre-C1 anticausal gate run (git 52153b3) measured
   the SAME GN null residual with the swapped quadrature and read offsets
   of only +0.012 / +0.019 sigma_F, which is how the swap masked the
   displacement. Consequence: the live-refit null offset (<= 0.6 sigma_F)
   is carried as a systematic on all Fisher-only numbers (worst-case
   <= +30% on a Fisher u95) and is covered by the K_dyn = 10 headline
   scaling with an order of magnitude to spare; the non-detection verdicts
   are unaffected. Paper B states the gate outcome as FAIL-with-diagnosis;
   template-shape validation at limit amplitude rests on G2b (PASS,
   <= 0.1% differential).
R7 in progress: Paper B numbers/claims regenerated from the corrected
   artifacts; README and joint-fit summary updated to the same numbers.

## Amendment G-2 (2026-07-12; committed BEFORE the gateG2 run it governs)

Post-R6 decomposition (span_mismatch_check, local, no live evaluations):
measuring the SAME baseline residuals res0 -- no GN refit at all -- in the
R6 gate harness's measurement span reproduces the entire "G1 offset":
span-only delta = +0.3439 / +0.6068 sigma_F at tau = 18 / 52 d, against
the recorded G1 offsets +0.3423 / +0.5987. The R6 harness's measurement
span was truncation-only (rank 70): it dropped the static-guard residual
direction that correction C4 added to every pipeline script -- C4 was
never propagated into the WSL gate harness. The genuine live-refit
displacement is the remainder: z_GN - z(gate span, res0) =
-0.0015 / -0.0081 sigma_F. The R6 adjudication's attribution (an
LS-vs-MAP refit walk projecting onto the causal template) is therefore
WRONG as stated; the walk exists (max|dth| ~ 14 steps) but its projection
onto the estimator is ~0.01 sigma_F, not 0.6.

Verified locally before this amendment: appending the guard-residual
direction to the gate span (rank 71) reproduces z_lin on res0 to machine
precision (delta ~ 1e-16) at both R6 anchors.

Registered plan (gateG2):
- Harness: identical to R6 except the measurement span is the guard-kept
  span (C4 completed in the harness); the GN absorb is unchanged.
- Anchors: tau = 2 d (the headline anchor, newly added), 18 d, 52 d.
  Inputs from the committed 10.8b artifact via
  scripts/make_gateG2_inputs.py -> sep_dynamic/gateG2_inputs.npz,
  gateG2_params.json.
- Criteria: unchanged from the R6 registration -- G1 |z_GN - z_lin| <=
  0.3; G2a (4 sigma_F injection, undamped GN) |dev| <= 0.5 sigma_F;
  G2b (u95_Kstatic injection, GN capped at 30x validated steps)
  |dev_rel| <= 2%. The dev formulas keep the beta_hat_lin reference for
  continuity with R6.
- Record rule: gateG2 supersedes R6 as the live-gate record of record;
  the R6 run remains recorded with the corrected diagnosis above. If any
  gateG2 criterion still fails, the failure is carried at its measured
  size as a genuine live-refit systematic -- no rescoring.
- Registered expectation (stated for honesty, not a criterion): offsets
  of order 0.01 sigma_F at 18/52 d per the decomposition; tau = 2 d is
  the genuinely new measurement (different quadrature mix: the outer
  carrier is near-instantaneous, g1 = 0.999).

### G-2 addendum: worst-phase headline template (committed BEFORE its run)

Review point (PR #2): the gateG2 tau = 2 d anchor exercised the
ZERO-PHASE template/params pair (t_off = 0, 10.8b artifact), while the
quoted headline is the 10.8e worst-phase row (t_off = 104.3494 d,
u95pm_K10 = 1.68e-9). Registered addendum: gate the exact worst-phase
pair. Inputs: T2wp/Tc2wp = causal templates at tau = 2 d,
t_off = 104.3494 d; params (beta_hat_lin = +0.0643 sigma_F, sigma_fisher,
u95 tiers) recomputed on the pipeline span by make_gateG2_inputs.py and
ASSERTED to reproduce the committed u95pm_K10 to < 1e-6. Criteria
unchanged (G1 0.3 sigma_F on |z_GN - z_lin|; G2a 0.5 sigma_F at
4 sigma_F injection; G2b 2% at the K_static-scaled amplitude). Runner:
scripts/wsl/wsl_gateG2wp.py -> sep_dynamic/sep_gateG2wp.json; the
adjudication derives its verdict from the recorded output.

### G-2 outcome (recorded post-run, 2026-07-12)

gateG2 (sep_gateG2.json): **9/9 PASS** under the unchanged registered
criteria, with the headline anchor directly gated.
G1 offsets: -0.0018 / -0.0152 / -0.0350 sigma_F at tau = 2 / 18 / 52 d
(tol 0.3). G2a: +0.0085 / +0.0562 / +0.0570 sigma_F (tol 0.5). G2b:
-0.109% / -0.024% / +0.128% (tol 2%). The live-refit systematic carried
on Fisher-tier numbers is therefore <= 0.035 sigma_F AS MEASURED
(replacing the earlier 0.6 sigma_F carry, which the span decomposition
proved to be the R6 harness artifact). The R6 artifact and its FAIL
verdicts remain recorded; sep_gateG_adjudication.json now carries the
full two-layer diagnosis (differential + span decomposition) and the
supersession. The live gate battery record of record: PASS at all
anchors including tau = 2 d.
