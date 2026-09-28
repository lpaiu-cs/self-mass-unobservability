# Request 11.2 — coverage and estimated-covariance result

Design commits: 9412142 (main simulation), 12d9f9a (estimated-covariance follow-through, registered after the main outcomes). Inputs and stage-1 geometry are frozen. Commands: `python verification/coverage_audit.py` and `python verification/estimated_covariance_audit.py`, four BLAS threads. Evidence: the corresponding JSON files under `outputs/research-completion/`.

Status: Proven. For an unbiased estimator b~Normal(beta,sigma^2) with known sigma, the symmetric Gaussian-mass interval [-U,U] defined by the manuscript's Equation (interval), K=1, has frequentist coverage at least 0.95. U is increasing in |b| and U(0)=1.959964 sigma. If beta>U(0), let h satisfy Phi((beta-h)/sigma)+Phi((beta+h)/sigma)=1.95. Failure occurs when |b|<h, with probability Phi((beta+h)/sigma)-Phi((beta-h)/sigma)<=0.05, since the first CDF is at most one. Symmetry handles negative beta; smaller |beta| is always covered. K>=1 enlarges U. A maximum over intervals containing the true grid cell cannot reduce coverage. This is not a theorem for estimated or misspecified covariance.

Status: Imported from prior work. The new registered simulation uses 8192 independent realizations per condition, 18 preselected lag/origin cells and nine signal amplitudes. Common random numbers correlate comparisons; rows must not be pooled as independent experiments. The generating model is the stored linear response, not a new timing integration. Compression/direct-projection checks pass. The oracle-covariance positive control's minimum U coverage is 0.94812, above the registered regression-screen threshold 0.93796.

| Status | Generating stress | Minimum truncated diagonal U coverage, K=1 | Minimum full diagonal U coverage, K=1 | Minimum oracle GLS U coverage, K=1 |
| --- | --- | ---: | ---: | ---: |
| Imported from prior work | White noise | 0.94678 | 0.94812 | 0.94812 |
| Imported from prior work | Omitted mean, norm 3 | 0.09558 | 0.94812 | 0.94812 |
| Imported from prior work | Omitted mean, norm 10 | 0 | 0.94812 | 0.94812 |
| Imported from prior work | Omitted mean, norm 30 | 0 | 0.94812 | 0.94812 |
| Imported from prior work | Extra Fourier RMS 0.25 | 0.93103 | 0.74695 | 0.94849 |
| Imported from prior work | Extra Fourier RMS 1 | 0.81323 | 0.58252 | 0.94897 |

Status: Imported from prior work. These minima are descriptive across tested cases, not simultaneous confidence guarantees. The worst full diagonal K=1 cell has 4772/8192 coverage, Wilson 95% interval [0.57180,0.59316], at lag 18 days, origin 192.20119 days, extra-Fourier RMS 1 and signal 20 full-space unit-noise sigma. Its actual estimator noise is 9.21 times the nominal unit-noise sigma, whereas its median residual scale is only 1.19. Scalar width rescaling does not describe this directional covariance.

Status: Imported from prior work. Truncated pointwise K=10 coverage can also fail: at lag 2 days, origin zero, omitted norm 30 and signal -20 full-space sigma, coverage is 0/8192 (Wilson upper endpoint 0.000469). This does NOT establish failure of the original K=10 grid envelope. The separately registered 512-realization complete-grid envelope tests at lags 2 and 200 days all give K=10 coverage at least 511/512; the smallest full K=1 envelope coverage is 428/512 under extra-Fourier RMS 1 at 200 days. A finite set of successful inflated envelopes does not certify uniform robustness.

Status: Proven. The stage-1 impossibility statement concerns an unbounded omitted mean at fixed nominal Gaussian width and nonzero estimator overlap. It does not establish a universal impossibility for every data-dependent residual-scale inflation rule. The actual re-estimated-scale failures above are numerical counterexamples to the particular tested procedure. This clarification narrows wording, not the registered computation or frozen outcomes.

Status: Imported from prior work. In the follow-through, the full 90-direction model estimates covariance sigma^2(I+a^2 LL^T) by REML, co-fitting both signal coefficients, with the predeclared Fourier 31–60 spectral variance j^-4 and fixed a grid. Across 486 conditions, the minimum K=1 U coverage is 7749/8192=0.945923, Wilson [0.940813,0.950615]. Minimum coverage by true a=0,0.25,1 is 0.947266, 0.947144, 0.945923. No fit reaches the upper a-grid endpoint. Signal-translation and zero-a ordinary-fit controls pass. This supports approximately nominal coverage within this covariance family; it is not exact coverage or independent validation of astrophysical noise.

Status: Imported from prior work. Fits to stored residuals at the same preselected cells select a=0.07368 or 0.08316. At lag 2 days and the full-space reference origin 174.27780 days, the local conditional U is 5.25094e-10 at K=1. These are local fits, not a newly searched envelope or a detection analysis. Lag/phase look-elsewhere calibration, chromatic/instrumental effects, nonlinear derivative errors and untested noise families remain outside this experiment.

Decision: full nuisance plus explicitly fitted covariance is the tested baseline. Preserve the original K=10 table as historical sensitivity data, and present the simulation as conditional method validation. Do not turn a safety multiplier into a universal SEP limit. Proceed to stage 3's force-level matching.

Classification: theorem progress (known-Gaussian coverage statement) and loophole progress (tested nuisance/covariance failure modes and a conditional remedy).
