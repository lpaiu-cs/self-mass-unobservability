# Projection Channel Audit

Status: Counterexample candidate. M4 fixes a concrete free-fall-style
observable channel and asks whether the projection factor is calibrated,
finite-dimensional, or arbitrary.

## Channel Definition

Status: Counterexample candidate. The minimal body-response variable is

```text
q_A(t) = delta m_A(t) / m_A^(0),
q_hat(z) = G(z) F_hat(z),
G(z) = c_Y + beta/(1 + tau_chi z),
beta = alpha c_chi.
```

Status: Counterexample candidate. The first projection channel is an immediate
linear acceleration-like readout,

```text
delta a_hat(z) = Gamma q_hat(z).
```

Status: Counterexample candidate. The second projection channel is a
linearized range-like response,

```text
ddot R + kappa^2 R = Gamma q_A(t).
```

Status: Proven. In frequency space with `z=i Omega`,

```text
R_hat(z) = [Gamma/(kappa^2 + z^2)] q_hat(z),
Lambda_R(z) = Gamma/(kappa^2 + z^2).
```

## Projection Verdicts

Status: Proven. If `Gamma != 0` is known or calibrated, the acceleration-like
readout preserves the same transfer up to an overall scale.

Status: Proven. If `Gamma` and `kappa` are known or calibrated, the range-like
readout can be deprojected by multiplying by `(kappa^2 + z^2)/Gamma`, so it
recovers `G(z)` exactly away from projection poles.

Status: Counterexample candidate. If `Gamma` and `kappa` are unknown but shared
across all sampled frequencies, the range channel introduces a finite nuisance
basis rather than an arbitrary frequency-local nuisance. The pole test remains
structural, but the exact sample-count theorem must include these nuisance
parameters or use independent calibration.

Status: Proven. The range projection does not cancel the relaxation pole. The
observed range transfer is

```text
T_R(z) =
  Gamma [c_Y(1 + tau_chi z) + beta]
  /
  [(1 + tau_chi z)(kappa^2 + z^2)].
```

At `z=-1/tau_chi`, the numerator is `Gamma beta`, so the relaxation pole
survives when `Gamma beta != 0` and `kappa^2 + tau_chi^-2 != 0`.

## Collapse Conditions

Status: Proven. Projection collapse occurs when `Gamma = 0`; the channel is
blind to `q_A`.

Status: Proven. Deprojection fails at a sampled projection pole,
`kappa^2 + z_k^2 = 0`, unless that resonance is separately modeled.

Status: Proven. A general projection can erase the relaxation pole if it has a
zero at `z=-1/tau_chi`.

Status: Proven. If `Lambda(z_k)` is an unconstrained complex nuisance at every
sampled frequency, then the observed transfer

```text
O_hat(z_k) = Lambda(z_k) G(z_k) F_hat(z_k)
```

can be fit point by point, and the pole relation is not tested.

## M4 Boundary

Status: Counterexample candidate. The range-like channel is a viable
projection dictionary only under one of these assumptions:

1. `Gamma` and `kappa` are known from the baseline channel.
2. `Gamma` and `kappa` are fitted as shared finite nuisance parameters.
3. The analysis uses ratios or deprojection constraints that remove `Gamma`
   and keep `kappa` fixed.

Status: Counterexample candidate. If none of these assumptions is justified,
M4 yields a projection no-go for this channel rather than an observable
dynamic-loophole claim.
