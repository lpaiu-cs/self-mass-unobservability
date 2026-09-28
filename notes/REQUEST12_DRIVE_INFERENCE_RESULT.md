# Corrected drive and simultaneous confidence regions

## Physical input and truncation

Status: Imported from prior work. The live parameter_set=6 convention gives masses (1.4378144085, 0.1975363853, 0.4101027069) solar masses and relative semimajor axes (4.7761915344e9, 1.7648750789e11) metres at Newtonian order. The source computes the inner inclination from the outer inclination plus delta_i in degrees, rather than using the unused apcosi_i parameter. The reconstructed outer mass agrees with the stored derived parameter to roundoff. Runtime source and input hashes are retained separately.

Status: Proven. Normalize the leading potential drive by Ustar=sum_k A_k, using delta U/c^2. Its positive amplitudes then sum to one. This choice defines beta; changing Ustar rescales beta and is not a physical improvement in sensitivity. The instantaneous response is co-fitted.

Status: Imported from prior work. Ustar=1.7499915637e-10 and normalized amplitudes are (0.2439492541, 0.6919568193, 0.0640939266). The corrected phases retain the earlier nonzero closure. Exact coplanar Kepler-potential grids at 128^2 and 256^2 give omitted-input RMS 3.374269% of the leading-drive RMS. This is not a bound on omitted timing response or additional scalar-force terms. The corrected leading-drive pointwise U values range from 4.49e-10 near 18 days to 1.20e-8 at 500 days for the specified covariance fit; these are conditional phenomenological coefficients, not a restored EOS bound.

## Joint-region inference

Status: Proven. Write the full-nuisance contrast data as y=C theta+noise, where theta has six real carrier coefficients, and use an orthonormal factorization C=Q R. Estimate the prespecified noise covariance after fitting all six coefficients. Translating y by any signal in col(C) leaves the fitted residual, covariance choice and white scale unchanged. Consequently the distribution of the coefficient estimation error is independent of the true carrier amplitudes, phases and lag, within this fixed linear mean/covariance model.

Status: Proven. A confidence region E for theta can be inverted through any declared family theta=W(phase,tau) b. If theta_true is in E, the true parameter point is in that inverse image. Thus projecting the entire region controls phase/lag selection without treating separate grid cells as independent tests. Evaluating only six lag sections is not a numerical computation of the entire projection. Intersections may be empty and must not be silently replaced by a zero upper bound.

Status: Proven. For known white-noise scale and covariance Sigma(a)=I+a^2 L L^T, a declared a<=a_max implies Sigma(a)<=Sigma(a_max) in positive-semidefinite order. GLS with Sigma(a_max) produces an estimation-error covariance bounded by its nominal covariance. Its six-dimensional quadratic error is stochastically bounded by chi-square(6), so the 0.95 quantile gives a conservative region for every a in the declared interval. This does not cover unbounded/misspecified covariance or an estimated white scale automatically.

Status: Imported from prior work. Calibration seed 2026090912 uses 8192 draws at each a=0,.25,1,4. The frozen region threshold is 12.8241766227, the maximum of the nominal chi-square threshold and the four registered 95% order statistics. Commit 6d7f3f3 freezes it before validation. Independent seed 2026090913 gives inclusion rates 0.951660, 0.954346, 0.952759 and 0.954102. A separate shallower Fourier-spectrum stress gives 0.959229. These finite experiments support the specified estimator in these cases; they do not prove coverage for every intermediate covariance or all astrophysical noise. The covariance endpoint is genuinely part of the experiment, not evidence that the true covariance is bounded by it.

Status: Imported from prior work. The recorded data's six-coefficient omnibus statistic is 16.3525, above the calibrated threshold. This null is that every residual carrier coefficient vanishes. It is not the null beta=0 with the physical instantaneous term freely fitted, and it is not evidence uniquely identifying a relaxation pole. Every reported physical lag section includes beta=0. At 2 days the joint-region |beta| endpoint is 4.10217e-10, at 500 days 8.31672e-9. Such section endpoints can be smaller than pointwise intervals because a joint-region intersection also spends distance on lack of fit in other coefficient directions. They must not be advertised as universally stronger constraints.

Status: Conjectural. Omitted harmonics/forces, numerical derivative changes and a free transient can move the mean outside the six-column space. Region coverage above does not automatically cover those extensions. Their dedicated audits are separate gates.

Programs: verification/physical_drive_completion.py and verification/simultaneous_inference.py. Classification: theorem progress and conditional inference progress.
