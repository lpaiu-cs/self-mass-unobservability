# Frozen settings before Request 13 production

## Timing accuracy pilot

Compare columns 21,22,27 at steps h and h/2 using four independent settings: control (integration 1e-16, mesh 250, inverse-time termination 1e-14 days), integration (1e-18,250,1e-14), mesh (1e-16,500,1e-14), inversion (1e-16,250,1e-17). Reproduce the original control and zero-shift recovery. The isolated source corrects bounds-check order in both forward and backward integration. No accuracy certificate is inferred from these comparisons.

## Stellar matching

Status: Counterexample candidate. Compute the unscalarized GR branch of massless DEF gravity, A(phi)=exp(beta phi^2/2), beta=-4, phi_infinity=0, with the published LAL SLy table. At zero background scalar the scalar perturbation decouples at linear order; the gravitational mass target is 1.4378144085 solar masses. This is a specified theory/EOS response calculation, not a claim that J0337 has these parameters or nonzero scalar charge.

Status: Imported from prior work. LAL's table reader specifies pressure then energy density in inverse square metres. The downloaded original table and reader source have SHA-256 provenance in sources/manifest.json. Use monotone log-pressure/log-energy interpolation with explicit domain checks. Solve TOV with two tolerance levels and center-radius checks, then the regular scalar radial equation. Match the static exterior exactly and the frequency-domain outgoing exterior at successively larger radii. First reproduce the beta=0 constant static solution and the published Khalil prototype coefficients. Do not use gravitational mass as baryon mass.

Status: Conjectural. The outgoing-wave scattering response can provide a dynamical matching check beyond equilibrium information. A low-frequency polynomial fit is accepted only over a demonstrated convergence window; fit-derived damping/inertia are not exact full-star mode data. Companion charges and physical timing mapping remain separate open gates.

## Preliminary nonlinear timing/noise production

All 28 timing parameters are evaluated through the live fixed-turn residual engine, with eccentricity disks e^2<1-1e-6, positive periods, positive mass ratio and projected orbital scales. Initial runs: zero, +1 and -1 pulse step at the longest gap, plus a displaced physical start for zero. The covariance is sigma^2 diag(TOA_error^2)+a^2 F diag(k^-gamma) F^T with 30 Fourier pairs, gamma in [0,7] and a/sigma in [exp(-10),exp(8)]. A constant offset is profiled. This is an explicit alternative to the earlier free deterministic Fourier nuisance envelope, not a retroactive replacement of its verdict.

Use the Gaussian determinant in the likelihood and jointly profile the white scale, red amplitude and slope on every live timing evaluation. Propose constrained Gauss-Newton/Broyden steps with a live likelihood line search; the initial 16-iteration production budget is a pilot, not a completion or convergence rule. A fresh full Jacobian, additional starts/assignments and injection checks are required for scientific promotion. Numerical precision and physical signal gates remain open.
