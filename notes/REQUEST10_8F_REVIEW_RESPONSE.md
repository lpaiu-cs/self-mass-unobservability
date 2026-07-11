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
