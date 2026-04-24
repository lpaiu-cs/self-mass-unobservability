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

Status: Conjectural. Sidebands become a next escalation only after adding nonlinear readout, nonlinear forcing, or a non-linear internal equation.
