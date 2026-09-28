# Request 11.1 — nuisance audit, registered design

Date: 2026-09-09. User authorization: audit omitted nuisance directions, then coverage, then physical matching. Baseline manuscript commit: 4897038. Before-task checkpoint: 4d5b34b. OSK's original project/theory/timing nodes have been corrected before this work.

This is a new analysis of previously inspected, public-derived stored arrays, not a blinded experiment. Prior information includes the 71/90 ranks and the tabulated full/truncated interval ratios. No new Nutimo build or timing integration is required. REQUEST10 data and verdicts remain immutable. Code is inspected before outcome computation; the design and implementation will be committed before running this new audit.

## Questions and fixed comparisons

1. Reconstruct the existing normalized nuisance matrix from its actual timing Jacobian metadata, offset, Fourier pairs and static-SEP guard. Reuse `sep_common` for input, projection and template conventions.
2. Report singular values, numerical rank, original weak-mode loadings and the dimension of the full-space component removed by the reference truncation. Distinguish mixed singular directions from individual physical parameters. Loading shares use normalized column coordinates and are not physical priors.
3. Compare relative cuts `1e-2, 1e-3, 1e-4, 1e-6, 0` with the same explicit guard. Independently check the full projector by QR and SVD; report orthogonality, guard retention and projector agreement. Do not identify a small nonzero singular value with a physically forbidden direction.
4. At lags `2,5,18,52,200,500` days, evaluate all original 4821 origins and report Gaussian K=1/K=10 envelopes, worst origins and conditional beta information. Match stored anchor values before interpreting new comparisons. Per-grid changes are sensitivity analyses, not newly calibrated empirical limits.
5. At the stored two-day reference origin and each new full/truncated envelope origin, measure leakage of a unit-norm omitted nuisance residual into beta after co-fitting the instantaneous column. Give the exact worst-case standardized bias for residual norm budgets `R=1,3,10,30` in whitened noise units.
6. Show ridge penalties on the omitted orthonormal residual coordinates for prior widths `a/s=0,1,3,10,infinity`. These are explicitly artificial prior sensitivity scenarios, not astrophysical priors. The end points must recover hard truncation and full nuisance marginalization respectively.

## Interpretation and stop rules

Status: Proven. For an unbiased linear amplitude estimator l under unit white noise, sigma=norm(l). An omitted orthonormal nuisance basis D gives worst bias `R*norm(D.T@l)/norm(l)` in sigma units for a residual of norm at most R. Finiteness and small singular values alone provide no bound on its physical amplitude.

Status: Conjectural. Hard truncation is defensible as a physical result only if independently supported parameter restrictions bound the omitted directions. No such prior will be invented from their small singular values or from the desired beta limit.

Abort interpretation if input digests differ from the unified revision manifest; columns/norms are invalid; projector residual/orthogonality exceeds 1e-7; the original K=10 anchors fail a relative 1e-5 reproduction tolerance. If singular values reach numerical precision, report numerical-rank uncertainty rather than silently promoting the 90-column result.

Output: a deterministic script, JSON/CSV evidence, and a result note classifying the supported nuisance treatment and what stage 2 must test. New coverage design is registered only after this stage is interpreted; physical matching follows coverage. Negative results and exact failure boundaries count as research progress.

Execution amendment: the first run stopped at import because the active Python has no SciPy; no arrays were loaded and no outcomes computed. Use the standard-library erf for the same normal CDF instead. Design, grid, thresholds and interval definition are unchanged.
