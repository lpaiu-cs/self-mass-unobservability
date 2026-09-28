# Request 11.2 — registered conditional coverage experiment

Date: 2026-09-09. Stage 1 result commit: 3506522. Before-task checkpoint: 871f039. This experiment uses the frozen linear response arrays and does not certify nonlinear timing or the actual astrophysical noise distribution. Its scenarios are stress tests, not fitted nuisance priors. No REQUEST10 artifacts are overwritten.

## Fixed design

- Seed: 2026090902. 8192 independent Gaussian realizations per fixed scenario, with common random numbers across estimators/amplitude comparisons. Report dependence of comparisons; do not pool repeated scenarios as independent observations.
- Lags: 2,5,18,52,200,500 days. True origins: 0 and the stage-1 truncated/full K=10 envelope origins (deduplicated).
- Injected beta: 0, ±2, ±5, ±20, ±50 times the **full-space unit-noise sigma at that true cell**. Report dimensional amplitudes. Large signals may exceed the prior nonlinear validation range; conclusions then concern the frozen linear model only. Instantaneous coefficient is zero in generation and freely co-fitted.
- Scenarios: white noise; deterministic omitted-nuisance mean of whitened norm R=3,10,30, with orientation opposing beta in the truncated estimator; additional Gaussian Fourier noise at frequencies 31–60, spectral variances proportional to j^-4 and unprojected weighted RMS 0.25 or 1 relative to white noise. These extra frequencies lie outside the fitted first 30 pairs. For beta=0 choose one fixed positive orientation of the omitted mean.
- Estimators: rank-71 diagonal weighting, full rank-90 diagonal weighting, and full rank-90 GLS using the **known generating covariance**. GLS is an oracle positive control, not a covariance learned from the observations.
- For each replicate, re-estimate the residual noise scale exactly as in the stored procedure (after nuisance projection, before fitting cY/beta). Calculate the two-column estimator and report coverage of the original symmetric Gaussian-mass interval at K=1 and K=10, and signed central intervals as a diagnostic. Also report known-scale coverage to test the analytic Gaussian result independently of noise-scale estimation.
- Use a finite orthonormal basis of the union of the omitted subspace, six response columns and projected extra Fourier columns, plus an independent chi-square residual norm with the remaining degrees of freedom. This is an exact sufficient-coordinate simulation for the stated linear Gaussian experiment, not a reduced sample of TOAs. Verify orthogonality, full/truncated estimators, covariance and direct-vs-compressed mean/noise-scale algebra on generated vectors.

## Registered-grid envelope check

At tau=2 and 200 days, the full-space worst origin, beta=+50 full-space sigma, and the white/R=30/extra-Fourier-RMS=1 scenarios, use 512 realizations and all original 4821 fitted origins. Report envelope coverage for rank-71 and full rank-90 at K=1,10. Use the same reference noise-scale estimator per realization. This tests the actual max-over-grid operation; it does not claim coverage for origins outside the registered domain. Pointwise failures need not imply envelope failures.

## Interpretation and checks

Status: Proven. For fixed known Gaussian covariance, an unbiased Gaussian beta estimator gives at least 95 percent repeated-sampling coverage for the interval [-U,U] defined by 95 percent Gaussian mass about beta_hat. Coverage is 100 percent for sufficiently small absolute true beta and tends to 95 percent at large signal. Taking a maximum over a grid containing the true origin cannot reduce this pointwise coverage. K>=1 only increases U under these fixed assumptions.

Status: Proven. An unbounded omitted mean coupled to beta defeats any fixed K. Stage 1 establishes that this coupling is nonzero for the reference truncation. The stage-2 finite stresses measure examples, not a universal amplitude bound on those means.

Post-result wording clarification (Request 11.2 result; no change to the registered runs): the preceding fixed-K impossibility is at fixed nominal Gaussian width. It is not a theorem for every data-dependent residual-scale rule. The registered simulations separately test the actual re-estimated-scale procedure.

Report Wilson 95 percent intervals for every coverage fraction. For a nominal-model regression screen, flag coverage below 0.95 - 5*sqrt(0.95*0.05/8192); inspect its cause without changing the seed, model or threshold. This screen is not a familywise certification or an astrophysical noise validation. A mismatch stress may fail; it must remain in the report.

Stop for altered frozen input hashes, compression residuals above 1e-7, nonpositive Gram/covariance, or failure of the known-Gaussian algebra controls. Preserve any failure and register implementation amendments before rerunning. Do not tune K to pass the simulated cases and then claim independent calibration. The next stage is physical matching regardless of whether the old interval passes or fails these conditional tests.
