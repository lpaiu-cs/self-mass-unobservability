# Amplitude-Weighted Resonant Design Theorem

Status: Counterexample candidate. This note upgrades the exact-in-e resonant
forcing audit from a harmonic-counting theorem to a usable-sample theorem.

## Exact-In-E Harmonic Coefficients

Status: Proven. For

```math
F(t)=(a/r(t))^p=(1-e\cos E)^{-p},
\qquad
M=E-e\sin E,
```

the positive cosine harmonic coefficient is

```math
A_k(p,e)=
\frac{1}{\pi}\int_0^{2\pi}
(1-e\cos E)^{1-p}
\cos[k(E-e\sin E)]\,dE .
```

Status: Proven. This is an exact-in-e definition of the forcing coefficient.
The symbolic audit evaluates it deterministically by periodic quadrature; the
theorem does not rely on a small-e truncation.

## Observable Weight

Status: Counterexample candidate. For a second-order internal mode, the
dimensionless harmonic weight is

```math
W_k
=
|A_k(p,e)|
\left|
\frac{\rho^2}
{\rho^2-k^2+i\delta k}
\right|
|\Lambda_k|,
```

where

```math
\rho=\Omega_{\rm wrap}/n,
\qquad
\delta=\gamma_\chi/(\mu_\chi n).
```

Status: Conjectural. A harmonic is usable if either an absolute calibrated
noise floor is beaten or, in the current dimensionless audit,

```math
W_k/W_{\max}\ge \eta .
```

Here `eta` is a relative amplitude cutoff.

## Amplitude-Breaking Criterion

Status: Proven. A harmonic-counting design is not sufficient. The usable set

```math
U_\eta=\{k: W_k/W_{\max}\ge \eta\}
```

must satisfy both

```math
|U_\eta|\ge N+2+K_\Lambda
```

and

```math
\exists k_- , k_+ \in U_\eta
\quad
k_-<\rho<k_+ .
```

Status: Proven. If either condition fails, the result is an
amplitude-underbudget no-go for the stated cutoff and comparator class.

Status: Proven. If `rho<=1`, the positive-harmonic single-system design still
cannot bracket the resonance, regardless of amplitude.

## Example Audit

Status: Counterexample candidate. The default audit uses `p=2`,
`rho=3/2`, damping `delta=0.2`, relative cutoff `eta=10^-3`, and harmonics
through `H=6`.

Status: Counterexample candidate. Under those dimensionless assumptions,
moderate eccentricities can supply enough usable harmonics for the
acceleration channel and, depending on projection, for the range channel. This
is not an empirical detection claim; it is a proof-of-design count under an
explicit amplitude cutoff.

## Collapse And No-Go Boundaries

Status: Proven. If high-k coefficients fall below the cutoff before the sample
budget is exceeded, the design is observationally underbudget even though the
formal exact-in-e dictionary contains infinitely many nonzero harmonics.

Status: Proven. If projection nuisance is arbitrary per harmonic, no finite
amplitude-weighted design establishes the shared second-order line shape.

Status: Counterexample candidate. The next concrete boundary is to replace
the dimensionless cutoff `eta` by a channel-specific calibrated noise model.

