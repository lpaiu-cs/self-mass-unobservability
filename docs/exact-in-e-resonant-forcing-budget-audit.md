# Exact-In-E Resonant Forcing Budget Audit

Status: Counterexample candidate. This note asks whether exact-in-e orbital
forcing can provide the usable frequency samples needed by the second-order
resonant comparator theorem.

## Harmonic Dictionary

Status: Conjectural. Treat exact-in-e eccentric forcing
`F(Y(t))=(a/r(t))^p` as providing positive harmonics

```math
\Omega_k=k n,\qquad k=1,2,\ldots,H,
```

up to a usable harmonic cutoff `H`. This is a counting model: it assumes the
relevant Fourier coefficients are nonzero and measurable after projection
calibration.

Status: Counterexample candidate. Define the dimensionless phase-wrap
location

```math
\rho=\frac{\Omega_{\rm wrap}}{n}
=
\frac{1}{n}\sqrt{\omega_\chi^2/\mu_\chi}.
```

## Budget And Bracketing Rule

Status: Proven. For comparator degree `N` and finite shared projection
nuisance count `K_Lambda`, the second-order line-shape theorem requires

```math
H \ge N+2+K_\Lambda .
```

Status: Proven. A single eccentric system brackets the phase-wrap denominator
only if there are usable harmonics on both sides of `rho`.

For non-integer `rho>1`, this requires

```math
H\ge \lfloor\rho\rfloor+1.
```

For integer `rho=r`, the sample at `r n` lies exactly on the denominator zero,
so a strict bracket requires

```math
H\ge r+1
```

and a lower harmonic exists only for `r>1`.

Status: Proven. If `rho<=1`, a single positive-harmonic dictionary cannot
bracket the phase-wrap denominator; one needs another system, another forcing
scale, or a lower-frequency drive.

## Minimum Harmonic Cutoff

Status: Proven. The exact-in-e resonance-aware minimum is

```math
H_{\min}
=
\max\left(N+2+K_\Lambda,\ H_{\rm bracket}(\rho)\right),
```

provided `rho>1`. If `rho<=1`, `H_min` does not exist for a single eccentric
system.

## Example Classifications

| Case | `rho` | Comparator/projection | `H_min` | Status | Verdict |
| --- | ---: | --- | ---: | --- | --- |
| near half-integer wrap | `3/2` | `N=1`, acceleration `K=1` | `4` | Proven | budget-breaking if first four harmonics are usable |
| near half-integer wrap | `3/2` | `N=1`, range `K=2` | `5` | Proven | budget-breaking if first five harmonics are usable |
| near half-integer wrap | `3/2` | `N=2`, range `K=2` | `6` | Proven | budget-breaking if first six harmonics are usable |
| integer wrap | `4` | `N=1`, acceleration `K=1` | `5` | Proven | need one harmonic above the exact hit |
| integer wrap | `4` | `N=1`, range `K=2` | `5` | Proven | budget and bracket both require five harmonics |
| sub-fundamental wrap | `3/4` | any single-system positive-harmonic design | none | Proven | no single-system bracket |

## Verdict

Status: Counterexample candidate. Exact-in-e forcing can in principle supply
enough samples for the second-order line-shape theorem, unlike the current
small-e `n,2n` dictionary.

Status: Counterexample candidate. This is only a conditional positive result:
the required harmonics must be usable with adequate amplitude and calibrated
projection, and the internal phase-wrap location must satisfy `rho>1`.

Status: Proven. If the resonance lies below the fundamental orbital frequency
or projection nuisance is arbitrary per harmonic, exact-in-e forcing does not
produce a finite-sample novelty claim.

Status: Counterexample candidate. The next empirical-theory boundary is
therefore not the comparator theorem. It is whether a concrete system supplies
the required `H_min` usable harmonics around the internal resonance.

Status: Counterexample candidate. The first amplitude-weighted version of that
usable-harmonic boundary is recorded in
[`amplitude-weighted-resonant-design-theorem.md`](amplitude-weighted-resonant-design-theorem.md).
