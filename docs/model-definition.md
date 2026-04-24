# Model Definition

Status: Imported from prior work. The model attacks A4, the assumption that the free-fall sector contains no orbital-timescale internal state variable.

## Drive Choice

Status: Conjectural. The MVP chooses one scalar drive `F(Y)` from the external-field invariant set `Y`.

Status: Conjectural. For the first derivation, `F(Y)` is represented by a deterministic monochromatic drive

```text
F(t) = F0 cos(Omega t)
```

where `F0` is small enough that the linear response model is meaningful.

Status: Conjectural. The simplest dimensional convention is to normalize `F(Y)` so that `F` is dimensionless. If the selected invariant has dimensions, absorb the reference scale into the definition of `F`.

Status: Counterexample candidate. The first concrete forcing template is
`F(Y(t)) = (a/r(t))^p`, a normalized scalar invariant sampled along a
small-eccentricity orbit.

Status: Proven. This template supplies harmonics at `n` and `2n` through
`O(e^2)` when `e != 0`, `p != 0`, and `p != -3`.

## One-State Relaxation Model

Status: Counterexample candidate. The minimal dynamic loophole model is

```text
tau_chi * d chi_A / dt + chi_A = alpha * F(Y)
m_A(Y, chi_A) = m_A^(0) * [1 + c_Y F(Y) + c_chi chi_A]
```

Status: Conjectural. `chi_A` is an internal visibility state, not an extra static primitive family.

Status: Proven. If `F` and `chi_A` are dimensionless, then `alpha`, `c_Y`, and `c_chi` are dimensionless, while `tau_chi` has dimensions of time.

Status: Proven. If `F` carries dimension `[F]`, then `alpha` has dimension `[chi]/[F]`, `c_Y` has dimension `[F]^{-1}`, and `c_chi` has dimension `[chi]^{-1}`.

## Parameters

| Symbol | Status | Meaning |
| --- | --- | --- |
| `m_A^(0)` | Conjectural | Baseline body mass entering the free-fall-style worldline response. |
| `Y` | Conjectural | External invariant or small invariant set used to build the scalar drive. |
| `F(Y)` | Conjectural | Chosen scalar drive, normalized to dimensionless form for the MVP. |
| `chi_A` | Counterexample candidate | One internal state with orbital-timescale relaxation. |
| `tau_chi` | Counterexample candidate | Internal relaxation time; novelty requires finite `Omega tau_chi`. |
| `alpha` | Conjectural | Drive-to-state susceptibility. |
| `c_Y` | Conjectural | Instantaneous static response coefficient. |
| `c_chi` | Counterexample candidate | State readout coefficient in the mass response. |

## Comparator EFT

Status: Conjectural. The strict static comparator contains only instantaneous sensitivity coordinates multiplying local powers of `F(Y)`.

Status: Conjectural. The broader local comparator may include derivative Wilson coefficients such as `dot F`, `ddot F`, and higher derivatives.

Status: Proven. The one-state relaxation model is exactly equivalent to an infinite derivative expansion only inside the formal inverse `(1 + tau_chi d/dt)^(-1)`.

Status: Proven. A finite derivative comparator of order `N` is represented in
frequency space by a polynomial `P_N(i Omega)`.

Status: Proven. For real Wilson coefficients and positive frequencies, this
polynomial splits into a real even channel and an imaginary odd channel in
`Omega^2`.

Status: Counterexample candidate. The frequency-sweep comparison requires the
same polynomial coefficients to apply across all sampled drive frequencies.

Status: Counterexample candidate. The observable projection is modeled as
`O_hat(Omega_k) = Lambda(Omega_k) G(i Omega_k) F_hat(Omega_k)`. The pole test
requires `Lambda` to be known, calibrated, or finite-dimensional rather than an
arbitrary complex nuisance at each sampled frequency.

Status: Counterexample candidate. The first concrete projection channels are
`delta a_hat = Gamma q_hat` and `ddot R + kappa^2 R = Gamma q_A`, representing
acceleration-like and range-like linearized free-fall readouts.

Status: Counterexample candidate. For the triple carrier bridge, the measured
carrier samples are written as

```text
O_k = Lambda_k(theta) G(i Omega_k) F_k,
Omega_k in {Omega_in, Omega_out, |Omega_in-Omega_out|}.
```

Status: Conjectural. The projection-nuisance realism gate is the distinction
between finite shared `theta` and arbitrary independent complex `Lambda_k`.
Only the former can leave the shared-`tau_chi` transfer relation as a
runtime-worthy target.

## Minimal Nonlinear Extensions

Status: Counterexample candidate. Before adding a second-order mode, M5 tests
minimal nonlinear first-order extensions:

```text
tau_chi dot chi_A + chi_A = alpha F(Y) + beta_F2 F(Y)^2
```

and

```text
m_A/m_A^(0) =
  1 + c_Y F + c_chi chi_A
  + lambda_Fchi F chi_A
  + lambda_chi2 chi_A^2.
```

Status: Proven. These nonlinear terms are sufficient to generate sidebands
under two-frequency forcing.

## Deferred Secondary Model

Status: Counterexample candidate. If first-order relaxation is insufficient, the next model is

```text
mu_chi * d^2 chi_A / dt^2
  + gamma_chi * d chi_A / dt
  + omega_chi^2 (chi_A - alpha F(Y)) = 0
```

Status: Conjectural. This second-order mode is deferred until the first-order model's collapse boundary is exhausted.
