# Request 10.8g: Second-Round Review Corrections (Text/Framing; No Data Re-Runs)

Status: Note. A second independent adversarial review (2026-07-12, of the
post-10.8f manuscript) produced eight findings. Every quoted number was
checked against the committed artifacts before any edit; **no quoted number
changes** in this round. The findings are label, framing, and scope
corrections, plus one new design-only measurement. Detection rules, limits,
and the gate record are untouched.

## Adjudicated findings and fixes

D1 (anchor-tier label; review #1). The manuscript's phrase "95% statistical
CL x K_dyn = 10 anchor" misdescribes the computation: the pipeline
(`sep_common.u95_of`, called as `u95_of(beta_hat, K * sigma_fisher)` in
`sep_phase_marg_10_8e.py`) evaluates each tier as the **95% interval of a
Gaussian centered on the fitted amplitude with the width inflated by K**,
i.e. `u95(beta_hat, K sigma_F)` (~ 1.96 K sigma_F when |beta_hat| is small)
— NOT K times the statistical 95% limit. Because the worst phase is tracked
per tier and the beta_hat term matters more at small widths, tier ratios
differ from K: at tau = 2 d, u95pm_K10 / u95pm_fisher = 6.0, and both
K-tiers share the common base u95pm_K10 / 10 = u95pm_K934 / 934.4
= 1.6795e-10 (verified from `sep_phase_marg_10_8e.json`). FIX: state the
tier definition once in Section 4.1, propagate the corrected label to the
abstract, Sections 4.2, 6.2, 7, and Appendix C. Numbers unchanged.

D2 (ramp calibration "exactly"; review #4-part). Section 3.2 called
`a_hat(days) = 0.0219` "exactly" the predicted `1/(2 x 24.077) = 0.0208`;
the residual is +5.3% (and +0.4% on the 0.502 vs 0.5 branch). FIX: quote
the residuals and note they are small against the factor ~24 separating the
two unit hypotheses; drop "exactly". Also make explicit that the
integral-kernel (secular) reasoning is used only for the convention/
normalization calibration, never in the production templates.

D3 (manuscript vs C12 numeric conflicts; review #4-part). The review found
the manuscript (window edge 1.74e-8; zero-phase 1.36e-9) conflicting with
the 10.8f C12 note (~1.2e-8; 1.16e-9). Adjudication against artifacts:
`sep_limit_curve_10_8e.tsv` row tau = 500 gives u95pm_K10 = 1.736e-8, and
`sep_joint_fit_10_8b.json` gives min_u95_Kdyn = 1.3647e-9 — the
**manuscript is correct**; the C12 values were pre-re-run expectations
superseded by the R2/R3 outcomes recorded in the same note. FIX: none in
the manuscript; this note records the supersession (the 10.8f note itself
is a committed record and is not edited).

D4 (causal/advanced discriminability; review #2). The phase-marginalized
anticausal control is as loud as the causal statistic (Z = 2.340 vs 2.285),
and the manuscript's "common to both conventions" reading was asserted, not
quantified. NEW MEASUREMENT (design-only; consumes template columns, TOA
times and weights, and the nuisance span — never the residuals):
`scripts/sep_quadrature_overlap_10_8g.py` ->
`sep_dynamic/sep_quadrature_overlap_10_8g.json`. In the estimator's
nuisance-projected metric, after marginalizing the co-fitted instantaneous
direction, the causal and advanced beta templates overlap at
|cos| = 0.47 (tau = 2 d, worst headline phase), 0.86-0.91 (tau = 5-18 d),
0.29 (52 d), and 0.93-0.99 with opposite sign (200-500 d, where the two
conventions span nearly the same line). FIX: report these numbers in
Section 6.1 and state the reading plainly — the channel has limited
leverage on the lag *sign*; the quoted limit bounds the amplitude of an
oscillating-Delta response in the causal single-pole convention, and the
advanced template functions as a common-mode control, not a lag-sign
discriminator.

D5 (headline-tier depth framing; review #3). "Three orders below the static
bound" was quoted only at the most favorable tier. FIX: Section 6.2 states
the depth per tier (3.0 dex at K_dyn = 10 truncated; 2.7 dex full-rank
bracket; 1.1 dex at K = 934) at tau = 2 d.

D6 (turn-lattice scope; review #6). The lattice indices were never defined
in the manuscript. From `turn_search_10_8d.py`: a candidate reassignment
shifts pulse numbers by the wrapped smooth ramp
`n1 P (t - tau0)/T + n2 P ((t - t_ref2)/T)^2 + phi0 P` over the full span
(16 phase offsets, 3 quadratic references per cell), so the scanned family
is smooth linear-plus-quadratic reassignments. An isolated slip introduced
in one observing gap and undone in another is not in this family. FIX:
define the indices and the scope in Section 5; the existing "beyond the
tested lattice is not claimed" caveat stands.

D7 (companion cross-reference). Paper A's corrected electric basis is
five-dimensional (gradient-kinematics correction, same date); Paper B's
references to the companion are dimension-agnostic ("finitely many static
sensitivity coordinates") and need no change. Recorded for completeness.

D8 (span-equality framing; review #8). The gate harness/pipeline
span-equality is evidenced by the recorded guard-kept reproduction of
`z_lin` to ~1e-16 at the gated anchors; the manuscript already states this.
No stronger standing invariant is claimed in this round; a config-level
span-equality assertion across the full (tau, phase) grid is listed as
future hardening. FIX: none (scope note only).

## Ordering

This note is committed together with the D4 measurement script/artifact and
before the manuscript edits it governs (same commit series, note first in
the diff); no data-touching stage is re-run, so no pre-registration gate is
triggered.
