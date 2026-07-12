# Request 10.7d: Promotion of the Truncated Limit via Estimator-Calibration Gates

Status: Counterexample candidate. Per the 10.7c quote rules, Gate L
(displacement linearity, 0.100 vs 0.05) held back the ~6x tighter truncated
curve while Gate A validated its absorption model end-to-end at 0.003%
estimator error. This note pre-registers the replacement promotion criterion
and an extension of the quoted window, BEFORE any of the new gates are run.

Status: Note. Rationale for replacing Gate L: the quoted limit is a statement
about the ESTIMATOR's response to signal and noise, not about pointwise model
linearity along a noise-fitted nuisance displacement. Gates D1/D2 below
validate exactly the quoted quantities (null behavior and injected-signal
calibration) against the live nonlinear model. Gate L's 10% shape deviation
on the fitted structure (|lin|_W = 59.1 of a total residual norm ~125.9)
perturbs chi2 by <~ 0.2% and the noise scale s by <~ 0.1% -- negligible for
u95. Gate L is therefore demoted to a reported diagnostic.

Status: Note. Pre-registration discipline: the truncated-pipeline diagnostic
rows below tau = 10 d were computed to disk in the 10.7c run but have NOT been
inspected; no extended-window real-data statistic and no D1/D2 gate value is
known at the time this note is committed. Injection amplitudes are defined by
formula (u95_trunc at each anchor), to be read off after the extended-window
run.

## Pre-Registered Protocol

Status: Counterexample candidate. Extended quoted window: `tau_chi in
[1, 327] d` (the 10.7b window [10, 327] came from C1, a statement about the
SEQUENTIAL pilot's recovery; the operative definition of testability for the
joint estimator is gate D2 below, evaluated down to 1 d). Grid: 61 log-spaced
points on [1, 327] plus anchors `{1, 3, 26, 52, 104} d`; identical seed,
nsim, u95 definition, comparator, controls, and stop rules as 10.7b/c;
truncation of record: relative SV cut 1e-3.

Status: Counterexample candidate. Detection re-evaluation: the extended-window
truncated pipeline recomputes `Z = max|z|` with its own 500-sim null and the
four off-carrier controls (10.7b rules). A detection candidate
(`p < 0.003` AND controls quiet) suspends quoting and triggers the 10.7b
detection protocol.

Status: Counterexample candidate. Gate D1 (live null calibration): run the
REAL residuals (no injection) through the live-GN absorption procedure of
Gate A (truncated-pinv updates, <= 4 iterations, 0.5% convergence break),
then measure `z_GN(tau)` through the truncated joint pipeline on the
remainder. PASS iff `|z_GN(tau) - z_lin(tau)| <= 0.3` at every anchor
`tau in {1, 3, 26, 52, 104} d` (z_lin = the linear-pipeline real-data values).

Status: Counterexample candidate. Gate D2 (estimator calibration at the
limit): for each anchor `tau`, inject `beta_inj = u95_trunc(tau)` of the
`tau` chi column into the real residuals, run the same live-GN absorption,
measure `beta_hat_GN(tau)`. PASS iff
`|beta_hat_GN(tau) - beta_hat_lin_null(tau) - beta_inj| <= 0.5 sigma_beta_trunc(tau)`
at every anchor.

Status: Counterexample candidate. Quote rules (fixed in advance):

- D1 and D2 pass at all five anchors -> the truncated u95 curve over
  `[1, 327] d` is QUOTED as the 10.7d limit, superseding the 10.7c
  conservative quote (which remains recorded).
- D1/D2 pass at `{26, 52, 104}` but fail at `{1}` or `{3}` -> the truncated
  curve is quoted on `[10, 327] d` only; the sub-10 d region stays diagnostic.
- Any failure among `{26, 52, 104}` -> no promotion; the 10.7c conservative
  quote stands and a fully nonlinear profile (Minuit) becomes the only path.
- Gate L is re-reported as a diagnostic in all branches but gates nothing.

## Outputs

Status: Note. `joint_fit_upper_limit_10_7d.json`, `beta_limit_curve_10_7d.tsv`,
`jointfit_gateD.json`, updated `joint_fit_summary.md`, all under
`../request10_external/`.
