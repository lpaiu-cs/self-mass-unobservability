# Request 10.7a External Nutimo Pilot: Decision Summary

Status: Imported from prior work. Contract: `notes/REQUEST10_EXTERNAL_NUTIMO_HANDOFF_PACKET.md`
(branch `lpaiu/minimal-nonlinear-sideband`). Labels used exactly as required:
`runtime-motivated`, `conditional`, `collapse`.

## Verdict

Status: Counterexample candidate. The dynamic-chi branch for the named
`J0337+1715` / `Nutimo` implementation class is

```text
runtime-motivated AND conditional  (NOT collapse)
```

No hard stop rule triggered. Promotion to a real external runtime experiment
is granted **conditionally** (conditions C1–C2 below).

## Stage results

Status: Proven (preflight). A pre-provisioned Nutimo runtime was declared and
exercised: Zenodo 13899771 release (code `nutimo.tar.bz2`, data + `Analysis`),
built in WSL Ubuntu-22.04 with g++-9 / tempo2 @ f9fd985 (2021-03) / Minuit2
5.34.14 / boost 1.55.0. Residuals and finite gradients are computable through
`python_Fittriple_interface` (~3 min per full residual evaluation, 12474 TOAs).
Baseline `parfile-planetGR-max-bestfit` reproduces the released residuals up to
a benign linear drift (0.67 us rms, corr(diff,t) = -0.988, carrier content
<= 0.11 us) attributable to clock/EOP runtime files absent from the public
tarball; a linear drift lies inside the fitted span (spinfreq, spinfreq1).
Details: `baseline_fit_summary.json`.

Status: Proven (configuration closure — PASS). `specialcase` is empty in the
baseline parfile: the `RN_PL` per-harmonic amplitude/phase special case is NOT
active on any target carrier. The circum-ternary planet is modelled as an
integrated 4th body (`integrator 3`, `*_extra1` finite physical parameters),
i.e. shared state/geometry projection, not carrier-local Fourier nuisance.
The `PL*GR` parfiles of the same release enable `RN_PL` and remain the
Request 10.6 collapse-class comparator; they were not used as baseline.
Carriers `Omega_in = 3.856137`, `Omega_out = 0.019200`,
`|Omega_in - Omega_out| = 3.836937` rad/day are positive, distinct and
nonresonant (ratio 200.84; 1:1, 2:1, 1:2 all excluded).
Details: `configuration_manifest.json`.

Status: Proven (named Jacobian rank gate — PASS, no collapse). The central
finite-difference residual Jacobian over the 28 named fitted parameters
(steps = 1 sigma from the released `MCMC-planetGR.npz` chain) was projected
onto the six real coordinates of the three carrier complex amplitudes.
Principal cosines between the 28-parameter nuisance span and the 6D carrier
space:

```text
[1.000000, 0.999998, 0.999979, 0.996707, 0.977524, 0.275334]
```

Effective projection rank is <= 5/6 at every threshold (5 at cos >= 0.9,
4 at cos >= 0.99, 3 at cos >= 0.999); rank 6/6 is never reached. The
unabsorbable direction is concentrated on the outer-dipole combination
carrier: per-coordinate absorption cos_dif = 0.856, sin_dif = 0.600 versus
0.9999 for the inner carrier. This is precisely the finite shared
phase/geometric manifold behaviour predicted by the Request 10.6 source
audit, now confirmed at runtime. Robustness: column-normalized rank 28/28
stable across SV cutoffs 1e-6..1e-12; recomputing the most nonlinear column
(`delta_i`) with a 10x smaller step moves the 6th cosine by only 2e-4.
Details: `carrier_projection_rank.json`, `finite_jacobian.npy`.

Status: Proven (dynamic-chi test column — PASS). Test columns
`G(z) = c_Y + beta/(1 + tau_chi z)` (c_Y = 0, beta = 1, unit drive, zero
phase) for tau_chi in [0.05, 327] d have W-metric absorption <= 0.9946 by the
finite nuisance span; the surviving fraction is 10-22% for every tau_chi.
The dynamic-chi column does NOT lie in the named finite nuisance span.
Details: `dynamic_chi_span_test.json`, `dynamic_chi_test_column.tsv`.

Status: Proven (minimal synthetic injection — CONDITIONAL). After a
linearized refit of the 28 standard nuisance parameters (absorbing 97.8-99.2%
of the injected column in W-norm): for injected tau_chi >= 10 d the surviving
carrier relation is still best described by one shared tau_chi against the
admitted real shared-coefficient derivative comparator (chi2 better by
1.9x-17x than real degree-2 at equal parameter count). For tau_chi <= 5 d the
relation is not recovered. Two honest limitations: (i) the recovered tau_hat
saturates near ~112 d (~2/Omega_out) because the post-refit remainder is
dominated by the single surviving carrier direction — sequential
project-then-fit loses tau identifiability; (ii) the deliberately generous
complex degree-1 comparator (4 real parameters vs 3) fits the post-refit
amplitudes marginally better (chi2 ratio 0.85-0.93) — the Request 10.3
counting boundary `K = N + 2 = 3` is empirically marginal after refit.
Details: `synthetic_injection_recovery.json`.

## Stop rules

Status: Proven. None triggered:
- dependency repair did not replace the science gate (environment provisioned
  once, time-boxed, from the named public release);
- no RN_PL-like soak-up on target carriers in the baseline configuration;
- carrier inventory complete (3/3 samples, nonresonant);
- projection rank never 6/6;
- dynamic-chi column outside the finite nuisance span for all tested tau_chi.

## Promotion conditions

Status: Counterexample candidate.

- C1 (sensitivity window): a real runtime experiment should target
  `tau_chi >~ 10 d`, strongest near `tau_chi ~ 1/Omega_out ~= 52 d`. The
  branch is not runtime-testable in this dataset for `tau_chi <~ 5 d` after
  standard-nuisance refit.
- C2 (joint fitting + decisive comparator): the runtime experiment must fit
  `(c_Y, beta, tau_chi)` jointly with the finite nuisance parameters (not
  sequentially, which biases tau_hat toward ~112 d), and must pre-register the
  complex degree-1 per-carrier comparator as the alternative to exclude,
  since it is the only tested comparator class not decisively beaten by the
  shared-tau relation in this pilot.

## Caveats

Status: Proven. (a) The refit is linearized (Gauss-Newton around the
baseline); a full nonlinear refit could shift absorption fractions at the
percent level. (b) Nutimo prints `Bad combination of quadrupole and
integrator type` for the inactive `quadrupole1 = 0` (unfitted) parameter —
no effect on this analysis. (c) The pilot uses unit drive amplitudes
(`Lambda_k F_k = 1`); physical drive amplitudes rescale carrier weights by
known factors and do not change span-membership conclusions. (d) The 2020
Zenodo release (3778978) was not needed; all inputs came from the 2025
release of record.
