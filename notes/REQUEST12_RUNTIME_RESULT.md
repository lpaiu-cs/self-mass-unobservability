# Live transient, derivative and pulse-count follow-through

Status: Imported from prior work. Archived and isolated rebuilt engines each reproduce the frozen zero-drive residual array exactly. Source, shared-library, interface, parfile and TOA SHA-256 values are recorded. The isolated exponential-capable engine is used for new production. REQUEST10 bytes are unchanged. The copied runtime retains the original timing theory and prescribed SEP pair-coupling patch; it is not a newly matched scalar-tensor timing theory.

## Dedicated transient response

Status: Imported from prior work. At tau=2,52,500 days, central differences at coupling amplitudes 1e-8 and 5e-9 supply six full residual-pair artifacts. After all 90 nuisance directions, derivative changes are 0.389791%, 0.005752% and 0.000732%, passing the registered 5% convergence gate. Maximum absolute perturbations are 0.231, 6.055 and 57.896 microseconds at the larger amplitude, below the 1000-microsecond wrap guard. A final zero-amplitude recovery passes. The exponential is defined at the integrator epoch and uses its calibrated time unit, not inserted directly into residuals.

Status: Imported from prior work. Co-fitting this measured transient with the corrected leading instantaneous and relaxed drive increases beta's fitted standard error by factors 1.00096--1.00523 across the nine registered lag/covariance combinations. Thus this particular added state does not cause a large variance penalty. This is a local linear result at three lags and the stated amplitudes; the six-coefficient region calibration does not automatically apply after adding a transient coefficient.

## Derivative convergence

Status: Imported from prior work. All 28 timing derivatives were recomputed at half their archived corrected steps; seven planet derivatives also have quarter-step values. The maximum half-step weighted relative change is 0.665177% (oman_extra1). Some quarter-step differences grow rather than shrink, so a uniform asymptotic finite-difference convergence claim is unsupported. Full difference vectors are retained, not only scalar summaries.

Status: Imported from prior work. Both full nuisance matrices have rank 90, but their maximum principal sine is 0.999913. In the archived normalization, the change norm 0.006857 exceeds the smallest singular value 1.12207e-6 by thousands of times. The fixed-matrix QR/SVD agreement does not settle the physical derivative accuracy. For the corrected physical template under diagonal covariance, replacing all timing derivatives by half-step values changes standard errors by factors 0.8680--0.8785 and shifts fitted coefficients. This is a sensitivity result; the half-step version is not promoted as a more accurate scientific baseline or used to replace previous verdicts.

Status: Proven. The sufficient projector bound requiring an error norm below the smallest singular value cannot be instantiated from these differences. A measured difference is not itself an upper bound on true error. The empirical derivative certificate remains incomplete.

## Pulse assignments and nonlinear domain

Status: Imported from prior work. All 201 assignments with at most two nonzero +/-1 steps on the ten longest gaps were evaluated (including the zero assignment). The weakest nonzero assignment remains the 223.37-day single step. Minimum delta-chi-square is 27.9406 for diagonal noise, 11.6183 at a=0.0831611 and 0.152887 at a=1. These are fixed-covariance finite-lattice comparisons, not complete pulse reconnection.

Status: Imported from prior work. The unconstrained full-nuisance compensation at a=1 requires a planet eccentricity parameter change delta eta=26.0400. At fractions .25,.5,1, the resulting eccentricities are 6.5805,13.0884,26.1073: outside the bound-orbit domain. The runtime would clamp them, so the first such computation was stopped and no clamped result is accepted. This rejects that particular linear displacement as a physical nonlinear fit; it does not rule out another nonlinear solution for the pulse assignment.

Status: Imported from prior work. The registered post-failure local fractions .001,.003,.01 yield weighted nonlinear discrepancies 0.4797%,1.5042%,4.0799%, passing the 5% local gate. They cover at most one percent of the full required compensation and therefore do not rescue the original assignment. Remaining errors after nuisance projection are 0.00885,0.00923,0.01487 unit-noise norm. Arbitrary multiple slips and a constrained nonlinear likelihood search remain incomplete.

Programs: verification/prepare_runtime12.py, runtime_completion.py, gap_pair_audit.py, runtime12_analysis.py. Classification: measured loophole progress with failed empirical-promotion gates retained.
