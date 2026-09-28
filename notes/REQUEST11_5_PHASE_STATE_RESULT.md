# Request 11.5 — phase, initial-state and numerical validation boundaries

Design fb269b9; coarse results and subsequent refinement registration recorded at 5b9dba2. Commands: `python verification/phase_state_audit.py`, then `python verification/phase_refinement.py`, with four BLAS threads. Evidence: `phase-state-audit.json` and `phase-refinement.json` under `outputs/research-completion/`. Coarse outcomes remain unchanged.

## Known-phase envelope versus unknown drive

Status: Proven. For the unit-drive benchmark with co-fitted instantaneous column, put d_tau=min_c ||w_beta-c*w_0|| in six-dimensional carrier-coefficient space. Orthogonal phase rotations preserve d_tau. With A the full-nuisance, fixed-covariance whitened map, I_beta(phi)>=s_min(A)^2*d_tau^2 for every phase. The beta estimator direction is in col(A), so

    U(phi) <= |beta_hat(phi)| + z_0.975*sigma_beta(phi)
           <= [||Proj_col(A) y|| + z_0.975*s] / [s_min(A)*d_tau].

Here s is the same pre-signal-fit residual scale used in the stored construction, recomputed in the specified covariance metric. This is an analytic upper envelope for all three continuous phases, evaluated numerically; it is not a tight supremum or an interval-arithmetic rounding certificate. Its statistical interpretation still depends on the covariance/estimator premises, separately from the deterministic inequality.

Status: Imported from prior work. The audit compares the 4821 historical origins, 72x72 orbital longitudes with zero closure, and a 24^3 independent-phase grid at six lags and three fixed covariance amplitudes. The diagonal legacy K=1 results reproduce Request 11.1. The grids are not nested: the coarse three-phase grid misses narrow long-lag maxima already present in the historical grid. Reporting its smaller maximum as a stronger bound would be incorrect.

Status: Imported from prior work. The registered follow-through starts from each grid's maximum and refines the three independent phases with 17 angular scales. It preserves all seed maxima; no move cap is hit and the largest final-scale relative improvement is 2.58e-9. Across the tested lag/covariance cells the local-refined envelope is 1.00073–1.03721 times the historical envelope. This small increase is observed local optimization, not proof that every global maximum was found. The analytic continuous upper envelope remains the guarantee on the fixed construction.

| Status | Lag (days) | Local-refined U, a=0.0831611 | Analytic all-phase upper envelope |
| --- | ---: | ---: | ---: |
| Imported from prior work | 2 | 6.30659e-10 | 1.42264e-8 |
| Imported from prior work | 200 | 3.39331e-8 | 4.69431e-8 |
| Imported from prior work | 500 | 8.26579e-8 | 1.13172e-7 |

Status: Imported from prior work. The fixed middle covariance amplitude is inherited from earlier local REML fits, not re-estimated as a global optimum here. Neither these unit-drive U values nor the phase grids restore the withdrawn physical-drive beta interpretation. Expanded-domain detection significance and coverage with re-estimated covariance are not newly calibrated in this experiment.

## Initial state: an exact missing-response boundary

Status: Proven. The independent homogeneous state is chi_h exp[-(t-t0)/tau]. Settling to a fraction epsilon of its initial amplitude requires elapsed time tau*log(1/epsilon). Without a bound on chi_h, decay alone gives no uniform absolute signal bound. For tau=500 days, one-percent and one-per-mille decay require 2302.59 and 3453.88 days; the stored 2987.86-day span ends at 0.00253967 of an initial amplitude placed at its start. The tau=2 terminal ratio underflows to zero in the floating-point JSON; mathematically it is exp(-2987.86/2), not exactly zero.

Status: Proven. Knowing a linear response on the constant drive and six sine/cosine inputs does not determine it on an exponential transient. The causal local differential operator

    A(D)=D product_k(D^2+omega_k^2)

annihilates all seven known inputs but acts on exp(-t/tau) with nonzero factor (-1/tau) product_k(1/tau^2+omega_k^2). L and L+eta*A(D) therefore agree on the stored response columns and can differ arbitrarily on that transient. This is insufficiency of the finite response record even within causal linear maps, not a claim that all such maps are realizable by the particular timing theory. A fixed, validated forward operator can supply the missing response; the six columns alone cannot.

Status: Imported from prior work. On the actual sampling, the fraction of an exponential drive outside the constant-plus-six-periodic input span is 0.7984–0.9973 over the six lags. This checks input-space independence; it is not a calculation of the corresponding TOA transient. No fake exponential residual column is fitted as if it were a physical timing response.

## Single-cycle slips across observing gaps

Status: Imported from prior work. All 565 adjacent-TOA gaps longer than one day were tested with an idealized permanent one-cycle step, using the frozen spin frequency and error weights. Each fit includes all 90 nuisance directions and all six freely adjustable harmonic coefficients, deliberately allowing more signal absorption than a two-coefficient pole model. Both signs are evaluated against the stored residuals. Direct and sufficient-coordinate residual quadratic forms agree.

Status: Imported from prior work. Under diagonal covariance, the weakest step leaves a unit-noise residual norm 6.43368 and minimum delta-chi-square 27.9406. Under the a=1 extra-Fourier stress, the same gap leaves norm 1.11486 and minimum delta-chi-square 0.152887. It is the 223.37084-day gap between relative days 2066.37293 and 2289.74377. This stress covariance is not the measured best noise model. The result neither establishes an actual slip nor uniformly rules one out. Arbitrary multiple slips, isolated session offsets and nonlinear pulse-number reconnection were not searched.

## Numerical precision is not derivative accuracy

Status: Imported from prior work. Full-space QR/SVD agreement remains 1.29e-10. The superseded coarse-step v1 and corrected v2 nuisance spaces have maximum principal sine 0.999776, and changed planet columns differ by up to 1.045 in weighted relative norm. This illustrates why a known invalid coarse-step derivative set is not an alternative valid prior. It is not an estimate of the remaining v2 error. The stored v2 half-step summaries range from 0.0016004 to 0.0248096; they do not supply all projected error vectors or an all-parameter/all-domain error bound.

Status: Proven. For a full-rank normalized nuisance matrix B and perturbation ||E||<=epsilon<s_min(B), equal-rank projector rotation is bounded by epsilon/[s_min(B)-epsilon] (capped at one). Indeed a unit vector in col(B+E) has a coefficient vector of norm at most 1/[s_min(B)-epsilon], and its residual outside col(B) is bounded by epsilon times that norm. A useful certificate thus requires an actual matrix-error budget relative to the weak singular scale. Numerical agreement for a fixed B cannot supply it.

Decision: retain the known-phase analytic envelope, report the local refinement, explicitly leave the initial-state response uncalibrated, and preserve the weak gap under the stress covariance. Do not call successful periodic gates a validation of transients, all pulse assignments or nonlinear derivatives.

Classification: theorem progress (continuous phase envelope and causal transient obstruction) and loophole progress (expanded-phase sensitivity, single-cycle gap and derivative-scope audit).
