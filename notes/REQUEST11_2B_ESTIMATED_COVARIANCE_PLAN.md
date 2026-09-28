# Request 11.2b — estimated covariance follow-through

This extension is registered after inspecting Request 11.2's planned results. It does not replace those results. Request 11.2 found diagonal full-space pointwise U coverage as low as 58.25% under its RMS=1 extra-Fourier stress, while known-covariance GLS retained its expected coverage. This motivates testing estimated rather than oracle covariance. No new family will be selected after this follow-through's outcomes.

Fixed model: the same rank-90 nuisance space and Fourier 31–60 covariance with spectral variance j^-4. Covariance is sigma^2*(I+a^2 L L.T), where L uses Request 11.2's normalization. Estimate a on the predeclared grid `{0} union geomspace(0.01,4,100)` and profile sigma by REML after fitting both instantaneous and beta columns. Include logdet(C), logdet(X.T C^-1 X), and residual degrees of freedom N-90-2. End-point hits are reported, never silently expanded.

Seed 2026090903, 8192 independent realizations per condition, common random numbers across comparisons. Same 18 lag/origin cells, amplitudes 0, ±2, ±5, ±20, ±50 in full diagonal unit-noise sigma, and generating extra-Fourier RMS a=0,0.25,1. Full nuisance absorbs deterministic omitted means; their exact invariance is checked separately, not counted as fresh independent simulations. K=1 only; do not tune K from these outcomes.

Report original Gaussian-mass U coverage, signed central coverage and Wilson intervals. The same nominal-model regression screen is used. A data-fitted covariance is not an independently known covariance, so its success is Monte Carlo evidence in this family, not a theorem of exact coverage. Keep all failures.

Fit this specified REML model to stored residuals at those same preselected cells and report a, sigma, beta, and conditional U. These are local linear-model fits, not a newly searched global timing solution. The selected cells and covariance family do not validate all lag/phase domains or astrophysical noise possibilities. No new detection statistic or universal SEP upper bound will be declared.

Controls: compare a=0 REML to ordinary full-space residual fitting; demonstrate signal-mean invariance of the REML objective; compare sufficient-coordinate and direct residual quadratic forms. Hash input arrays and retain Request 11.2 output unchanged. Stage 3 begins after this extension's result is recorded.
