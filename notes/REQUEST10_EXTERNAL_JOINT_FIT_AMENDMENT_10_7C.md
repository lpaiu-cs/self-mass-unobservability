# Request 10.7c: Amendment — Validated-Span Limit After the 10.7b Stop Rule

Status: Counterexample candidate. The 10.7b pre-registered linearization check
FAILED (relative deviation 5.93 vs < 0.05) at the full-rank fitted nuisance
displacement, so per the 10.7b stop rules no 10.7b limit is quoted. Diagnosis
(`../request10_external/jointfit_lindiag.json`, `jac_v2/rederive_report.json`):

- The 21 non-planet parameter directions are linear to 0.3% at the fitted
  displacement; `delta_i` alone is linear to 5e-5.
- The 10.7a planet-block (`*_extra1`) Jacobian columns used steps of
  0.3 sigma_MCMC, which for these weakly-constrained elements are
  nonperturbative (tasc_extra1 step = 2.8 planet periods; oman step = 1 rad):
  those columns were secants, not derivatives. They were re-derived with
  perturbative steps (per-column half-step linearity 0.2-2.5%), giving
  `finite_jacobian_v2.npy`.
- With v2 columns the 10.7a Stage-2 gates are UNCHANGED (principal cosines
  match v1 to 5e-4; effective rank still <= 5/6), and the 10.7b full-rank
  pipeline numbers are unchanged to < 1%.
- However the full-rank beta marginal is dominated by three near-null nuisance
  directions (rel sv <= 5e-4) whose realization requires planet-element
  excursions of order the element values (period-wrapping tasc shifts) — the
  regime where the response is demonstrably nonlinear/oscillatory. Removing
  them (SVD cut `1e-3`, rank 26/29) preserves the fitted residual structure to
  1.7% but shrinks `sigma_beta(52 d)` from 27.4 to 4.59 us.

Status: Note. Pre-registration discipline: at the time this note is written,
the truncated-span REAL-DATA statistics (`beta_hat_trunc`, `Z_trunc`) have NOT
been computed, and neither WSL gate below has been run. Only truncated
`sigma_beta` values (data-independent) and the full-rank results are known.

## Pre-Registered 10.7c Protocol

Status: Counterexample candidate. Truncation of record: relative singular
value cut `1e-3` on the weighted, column-normalized nuisance block (29
columns; rank 26). Sensitivity at cuts `{3e-2, 1e-2, 3e-3}` is reported
alongside. The truncated pipeline is otherwise IDENTICAL to 10.7b
(same grids, seed, nsim, u95 definition, comparator, controls, stop rules).

Status: Counterexample candidate. Gate L (linearity of the truncated
displacement): full-Nutimo evaluation at the truncated `tau = 52 d` fitted
displacement `dtheta_trunc`; PASS iff
`||res(p0+dtheta_trunc) - res0 - J_v2 dtheta_trunc||_W / ||J_v2 dtheta_trunc||_W < 0.05`.

Status: Counterexample candidate. Gate A (reality of absorption): inject
`beta_test = 4 x sigma_beta_trunc(52 d)` of the `tau = 52 d` chi column into
the real residuals; attempt to absorb it with the REAL timing model by
Gauss-Newton iteration (J_v2-truncated update steps, live Nutimo residuals,
<= 4 iterations). Measure the surviving carrier signal
`z_true = beta_hat / sigma_beta_trunc` through the truncated joint pipeline on
the post-GN remainder. Predictions: truncated-linear `z ~= 4.0`; full-rank
linear `z ~= 0.66`. PASS iff `z_true >= 0.7 x 4.0 = 2.8` (the real model
cannot absorb the signal much beyond the truncated-linear prediction).

Status: Counterexample candidate. Quote rules (fixed in advance):

- Gate L PASS and Gate A PASS -> quote the truncated-span u95 curve as the
  10.7c upper limit (tighter; validated span).
- Gate L PASS and Gate A FAIL -> quote the full-rank u95 curve as a
  conservative upper limit (the full linear span over-absorbs by construction,
  so its limit is valid-conservative; documented as such).
- Gate L FAIL -> quote the full-rank u95 curve as a conservative upper limit
  ONLY with the over-absorption argument stated, and flag a 10.7d for a fully
  nonlinear profile; the truncated curve is NOT quoted.
- In every branch the detection statement is re-evaluated on the truncated
  pipeline with the same global-null calibration as 10.7b before any limit is
  quoted; a truncated-pipeline detection candidate suspends quoting and
  triggers the 10.7b detection protocol instead.

## Outputs

Status: Note. `joint_fit_upper_limit_v2trunc.json`, `beta_limit_curve_v2trunc.tsv`,
`jointfit_gateLA.json`, updated `joint_fit_summary.md`, all under
`../request10_external/`.
