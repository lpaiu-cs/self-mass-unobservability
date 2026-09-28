# Request 11.1 — result

Design/implementation commits: e790643, b8d931b. Input manuscript: 4897038. Command: `python verification/nuisance_audit.py` with four BLAS threads. Evidence: `outputs/research-completion/nuisance-audit.json` and `nuisance-intervals.csv`.

Status: Imported from prior work. The frozen normalized matrix has numerical rank 90, minimum relative singular value 4.18436e-7, full SVD/QR subspace residual 1.29e-10, and reference truncated rank 71. The removed subspace has dimension 19; its directions are not floating-point null directions of these arrays. This is numerical validation of the stored matrix, not independent certification of the finite-difference derivative accuracy.

Status: Imported from prior work. The weak singular modes mix actual timing parameters with Fourier columns. Examples include delta_i/masspar_p/period_o; proper motion; sky position; planet parameters mixed with low Fourier terms; spinfreq/spinfreq1; and period_i/static-SEP/period_o. Their normalized-coordinate loadings are listed in the JSON. Twenty low singular vectors lose a total of nineteen dimensions because the explicit guard recovers one linear combination. There are not nineteen individually excluded physical parameters.

| Status | Relative cutoff | Retained rank | K=10 U at 2 d | K=10 U at 200 d |
| --- | ---: | ---: | ---: | ---: |
| Imported from prior work | 1e-2 | 67 | 1.66816e-9 | 7.08946e-9 |
| Imported from prior work | 1e-3 | 71 | 1.67953e-9 | 7.15521e-9 |
| Imported from prior work | 1e-4 | 80 | 3.43587e-9 | 5.03656e-8 |
| Imported from prior work | 1e-6 plus guard | 90 | 3.53361e-9 | 1.24127e-7 |
| Imported from prior work | Full | 90 | 3.53361e-9 | 1.24127e-7 |

Status: Imported from prior work. All original full/truncated K=10 anchors reproduce within the registered 1e-5 relative tolerance. The full-space finite-grid envelope at 500 d is 3.02245e-7. The residual noise scale changes by only about 0.02 percent between ranks 71 and 90; the interval changes arise mainly from signal/nuisance geometry rather than scalar noise rescaling. The worst origins also move substantially.

Status: Proven. If l is the truncated beta estimator and D an orthonormal basis for the omitted residual directions, the worst normalized bias for norm budget R is `R*norm(D.T l)/norm(l)`. At fixed data model, removal of D by full nuisance fitting eliminates this particular bias exactly. At fixed nominal Gaussian width, no finite K protects the truncated estimator uniformly against an unbounded omitted mean. Request 11.2 separately tests the actual data-dependent residual-scale procedure.

Status: Imported from prior work. At the old two-day reference origin, unit omitted-residual norm gives maximum bias 0.5556 sigma. At the full-space envelope origins, the corresponding coefficients are 0.9323 at 2 d and about 0.99–0.999 at the longer lags. Thus an omitted residual with norm 30 in whitened units can shift the truncated estimate by roughly 17–30 sigma. These are adversarial linear-model bounds, not measurements of such a physical residual. The observed omitted-component norm is 4.108; it does not provide an externally justified bound on future nuisance means.

Status: Imported from prior work. Artificial Gaussian prior widths on the omitted residual coordinates smoothly interpolate between hard removal and full marginalization. At the old two-day origin, the unit-noise conditional sigma increases from 7.686e-11 to 9.266e-11 while beta moves from 5.51e-12 to -9.15e-11. These coordinate priors have no established astrophysical calibration.

Decision: use the full stored 90-direction matrix as the primary baseline for stage 2; retain rank-71 only as a registered sensitivity/undercoverage comparison. Do not convert numerical smallness into a physical zero or select a prior for a favorable limit. No externally justified prior excluding the nineteen-dimensional space has been established.

Limits: actual derivative error and nonlinear nuisance curvature are not certified by an SVD/QR agreement. Existing planet-column finite-step metadata already reports nonzero differences. Coverage on frozen arrays will therefore be explicitly a linear-model validation; it cannot certify the nonlinear instrument/astrophysical likelihood. No runtime or frozen REQUEST10 result was changed.

Classification: theorem progress (explicit bias boundary) and loophole progress (conditional nuisance treatment).
