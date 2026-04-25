# Non-Adiabatic Regime

Status: Counterexample candidate. The non-adiabatic regime is the finite `Omega tau_chi` regime in which the internal state cannot be replaced by a static sensitivity coefficient without losing phase information.

## Monochromatic Solution

Status: Proven. For

```text
tau_chi dot chi + chi = alpha F0 cos(Omega t)
```

the steady periodic solution is

```text
chi(t) = A cos(Omega t) + B sin(Omega t)
```

with

```text
A = alpha F0 / (1 + Omega^2 tau_chi^2)
B = alpha F0 Omega tau_chi / (1 + Omega^2 tau_chi^2).
```

Status: Proven. The amplitude and phase-lag form is

```text
chi(t) = R cos(Omega t - phi)
R = alpha F0 / sqrt(1 + Omega^2 tau_chi^2)
phi = arctan(Omega tau_chi).
```

Status: Proven. With the complex convention `F = Re[F0 exp(i Omega t)]`, the transfer function is

```text
H(Omega) = alpha / (1 + i Omega tau_chi),
arg H = -arctan(Omega tau_chi).
```

## First New-Observable Candidate

Status: Counterexample candidate. The mass readout contains the quadrature coefficient

```text
delta(m/m0)|_quadrature
  = c_chi alpha F0 Omega tau_chi
    / (1 + Omega^2 tau_chi^2).
```

Status: Counterexample candidate. This is a new observable relative to an instantaneous static sensitivity basis because the static basis has no delayed component.

Status: Proven. This is not a standalone proof of novelty relative to a broader derivative EFT, because a single `dot F` coefficient can mimic the leading quadrature at one frequency.

Status: Counterexample candidate. The stronger observable is the frequency-dependent pair

```text
A(Omega), B(Omega)
```

with a common `tau_chi`, because no finite instantaneous static coefficient vector reproduces the full rational response over varying `Omega`.

## Basis Audit Result

Status: Proven. The full readout transfer is

```text
G(Omega) = c_Y + c_chi alpha/(1 + i Omega tau_chi).
```

Status: Proven. A single known drive frequency can be fit by

```text
a0 F + a1 dot F
```

with frequency-local coefficients, so one quadrature point alone is not enough
against a derivative-EFT comparator.

Status: Counterexample candidate. The nonadiabatic observable target is the
shared rational transfer relation across changing `Omega`, not the mere
existence of a sine quadrature in one monochromatic experiment.

## Frequency-Sweep Distinguishability

Status: Proven. A finite derivative comparator of degree `N` has the form
`P_N(i Omega) = sum_{n=0}^{N} d_n (i Omega)^n`.

Status: Proven. With freely chosen shared coefficients, this comparator can
match the relaxation transfer at `N + 1` distinct sampled frequencies if the
coefficients are allowed to be complex.

Status: Counterexample candidate. With real derivative coefficients and
positive drive frequencies, matching `floor((N+1)/2)+1` distinct frequencies is
already an exact non-adiabatic distinguishability test.

Status: Counterexample candidate. With the stronger complex-coefficient
comparator, matching `N + 2` distinct frequencies is the conservative exact
test, because the residual is proportional to the pole strength
`alpha c_chi tau_chi^(N+1)`.

## Forcing Dictionary Result

Status: Counterexample candidate. The first concrete frequency source is a
small-eccentricity orbital harmonic expansion of a normalized invariant
`(a/r)^p`.

Status: Proven. The `n` and `2n` harmonics give the first real-coefficient
degree-`1` and degree-`2` finite-derivative obstruction in one system, provided
the observable projection is not arbitrary frequency by frequency.

## Projection Channel Audit

Status: Counterexample candidate. The acceleration-like channel has constant
projection, so it cannot create or hide a frequency-dependent pole except by
setting the channel scale to zero.

Status: Proven. The range-like channel transfer is

```text
T_R(z) =
  Gamma [c_Y(1 + tau_chi z) + beta]
  /
  [(1 + tau_chi z)(kappa^2 + z^2)].
```

Status: Proven. At the relaxation pole `z=-1/tau_chi`, the numerator is
`Gamma beta`; therefore the pole survives unless the channel is blind or the
projection has a zero at the same point.

## Triple Carrier Bridge

Status: Counterexample candidate. Existing triple inner/outer carriers sample
the same transfer at `G(Omega_in)` and `G(Omega_out)`.

Status: Counterexample candidate. This two-carrier inventory breaks the
physical real shared-coefficient degree-`1` and degree-`2` derivative
comparators when `Omega_in != Omega_out`, `F_in F_out != 0`, and
`beta tau_chi != 0`.

Status: Proven. It does not by itself break higher real derivative order or
complex degree-`1` comparators. The bridge is therefore a low-order positive
route, not a final universal theorem.

## Triple Three-Carrier Inventory

Status: Counterexample candidate. Adding the existing GR outer-dipole
combination carrier gives the candidate set

```text
Omega_in, Omega_out, |Omega_in - Omega_out|.
```

Status: Counterexample candidate. For nonresonant positive frequencies and
nonzero finite-dimensional projection, these three samples distinguish real
degree `N <= 4` derivative comparators and complex degree `N <= 1` polynomial
comparators.

Status: Proven. Real degree `N >= 5`, complex degree `N >= 2`, and arbitrary
per-carrier projection nuisance still require more carrier information or an
order/prior assumption.

## Triple Projection-Nuisance Realism

Status: Counterexample candidate. The three-carrier bridge observes

```text
O_k = Lambda_k(theta) G(z_k) F_k,
z_k = i Omega_k.
```

Status: Counterexample candidate. If `Lambda_k(theta)` is calibrated, a common
real scale, or a finite shared projection model, the dynamic-chi transfer law
is not automatically reduced to a frequency-local nuisance. The result remains
conditional for unknown finite projection models because the projection
parameters consume rank.

Status: Proven. If `Lambda_k` is an arbitrary complex nuisance independently
assigned to each carrier, then

```text
Lambda_k = O_k/(G(z_k)F_k)
```

fits every carrier pointwise and removes the shared-`tau_chi` target.

## Physical Projection-Manifold Gate

Status: Counterexample candidate. A physical outer-dipole projection supplies
a phase-locked combination carrier:

```text
Lambda_in  = A_in  exp(i phi_in)
Lambda_out = A_out exp(i phi_out)
Lambda_c   = A_c   exp(i(phi_in - phi_out)).
```

Status: Proven. The corresponding real Jacobian rank is `5`, while arbitrary
complex projection over the same three carriers has rank `6`.

Status: Counterexample candidate. This means the dynamic-chi three-carrier
bridge is runtime-worthy only for timing models that keep the combination
phase and geometry tied to the inner and outer carrier phases.

Status: Proven. A model that fits the combination carrier as an independent
complex amplitude destroys the bridge before any external runtime is useful.

## Large-Frequency Behavior

Status: Proven. For `x = Omega tau_chi >> 1`,

```text
A = alpha F0 / x^2 + O(x^-4)
B = alpha F0 / x + O(x^-3)
R = alpha F0 / x + O(x^-3)
phi = pi/2 - 1/x + O(x^-3).
```

Status: Counterexample candidate. The high-frequency limit is dominated by quadrature and amplitude suppression, which is qualitatively different from a constant static sensitivity coefficient.

## Two-Frequency Result

Status: Proven. For the linear model and linear readout, a drive

```text
F(t) = F1 cos(Omega1 t) + F2 cos(Omega2 t)
```

produces only `Omega1` and `Omega2` in `chi(t)`.

Status: Proven. The minimal one-state linear model does not generate sidebands.

Status: Proven. The minimal nonlinear drive extension

```text
tau_chi dot chi + chi = alpha F + beta_F2 F^2
```

generates sidebands at `2Omega1`, `2Omega2`, `Omega1+Omega2`, and
`|Omega1-Omega2|`.

Status: Proven. The minimal nonlinear readout extension

```text
lambda_Fchi F chi + lambda_chi2 chi^2
```

generates sum and difference sidebands. These sideband phases inherit
`tau_chi` through the linear `chi` coefficients.

Status: Counterexample candidate. Sideband generation is a stronger observable
class than one-frequency quadrature because an LTI projection cannot create
frequencies absent from the input.
