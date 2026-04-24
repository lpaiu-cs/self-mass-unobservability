# Minimal Nonlinear Sideband Test

Status: Counterexample candidate. M5 tests the smallest nonlinear extensions
after the linear one-state model failed to generate sidebands.

## Baseline Failure

Status: Proven. The linear model

```text
tau_chi dot chi + chi = alpha F,
m/m0 = 1 + c_Y F + c_chi chi
```

is linear time-invariant. For

```text
F(t) = F1 cos(Omega1 t) + F2 cos(Omega2 t),
```

it produces only `Omega1` and `Omega2`.

## Nonlinear Drive Extension

Status: Counterexample candidate. The first minimal nonlinear internal equation
is

```text
tau_chi dot chi + chi = alpha F + beta_F2 F^2.
```

Status: Proven. For two-tone forcing, the nonlinear drive contains

```text
F^2 =
  (F1^2/2) [1 + cos(2 Omega1 t)]
  + (F2^2/2) [1 + cos(2 Omega2 t)]
  + F1 F2 cos((Omega1 + Omega2)t)
  + F1 F2 cos((Omega1 - Omega2)t).
```

Status: Proven. Therefore `beta_F2 F^2` drives new frequencies
`2 Omega1`, `2 Omega2`, `Omega1 + Omega2`, and `|Omega1 - Omega2|` whenever
the corresponding amplitudes are nonzero.

Status: Proven. A harmonic drive `S cos(nu t)` in the relaxation equation
produces

```text
chi_nu(t) =
  [S/(1 + nu^2 tau_chi^2)] cos(nu t)
  + [S nu tau_chi/(1 + nu^2 tau_chi^2)] sin(nu t).
```

Thus every nonzero nonlinear-drive sideband inherits its own relaxation phase
lag.

## Nonlinear Readout Extension

Status: Counterexample candidate. The second minimal nonlinear extension keeps
the internal equation linear but reads out

```text
m/m0 = 1 + c_Y F + c_chi chi
       + lambda_Fchi F chi + lambda_chi2 chi^2.
```

Status: Proven. If

```text
chi_j(t) = A_j cos(Omega_j t) + B_j sin(Omega_j t),
A_j = alpha F_j/(1 + Omega_j^2 tau_chi^2),
B_j = alpha F_j Omega_j tau_chi/(1 + Omega_j^2 tau_chi^2),
```

then `F chi` and `chi^2` both generate components at
`Omega1 + Omega2` and `Omega1 - Omega2`.

Status: Proven. For the sum frequency, the nonlinear readout coefficients are

```text
C_sum =
  lambda_Fchi (F1 A2 + F2 A1)/2
  + lambda_chi2 (A1 A2 - B1 B2),

S_sum =
  lambda_Fchi (F1 B2 + F2 B1)/2
  + lambda_chi2 (A1 B2 + B1 A2).
```

Status: Counterexample candidate. A nonzero pair `(C_sum, S_sum)` is a genuine
sideband candidate because `Omega1 + Omega2` is absent from the linear input
and absent from any linear static/derivative response to that input.

## Orbital-Harmonic Specialization

Status: Counterexample candidate. For the M3 drive

```text
F(Y(t)) = (a/r(t))^p
```

the `O(e^2)` input contains `n` and `2n`. Nonlinear mixing of those two
harmonics generates `3n`, which is absent from the linear `O(e^2)` input.

Status: Proven. The nonlinear-drive mixed `3n` forcing amplitude is

```text
beta_F2 * [p e] * [p(p+3)e^2/4]
  = beta_F2 p^2 (p+3) e^3 / 4.
```

Status: Proven. This `3n` sideband collapses when

```text
beta_F2 = 0,  e = 0,  p = 0,  p = -3.
```

## Projection Robustness

Status: Proven. A linear time-invariant projection cannot create a frequency
absent from its input. It can only multiply each existing harmonic by
`Lambda(i nu)`.

Status: Counterexample candidate. Therefore a detected sideband at
`Omega1 + Omega2`, `|Omega1 - Omega2|`, or `3n` is stronger than a one-frequency
quadrature anomaly: it cannot be produced by an arbitrary linear projection of
the original linear response unless the projection itself is time-varying or
nonlinear.

Status: Proven. Projection can still hide a true sideband if
`Lambda(i nu_sideband)=0`, or if the sideband frequency is outside the measured
band.

## M5 Boundary

Status: Counterexample candidate. The minimal nonlinear extension succeeds as
a sideband loophole if at least one nonlinear coefficient is nonzero and the
selected forcing has at least two nonzero harmonics.

Status: Proven. It collapses back to the linear no-sideband result when

```text
beta_F2 = lambda_Fchi = lambda_chi2 = 0
```

or when the drive has only one usable harmonic and no self-harmonic is retained
in the observable band.
