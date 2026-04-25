# Observable Targets

Status: Imported from prior work. A static parameter redefinition does not count as theorem-breaking novelty.

## What Does Not Count

Status: Proven. An absorbed secular shift does not count. Example: if `Omega = 0`, then `chi = alpha F` and the readout is equivalent to `c_Y -> c_Y + c_chi alpha`.

Status: Proven. A purely in-phase single-frequency amplitude shift does not count by itself, because it can be fit by a static sensitivity coefficient at that frequency.

Status: Proven. A one-frequency quadrature term is static-basis-degenerate if the comparator EFT allows an unconstrained local `dot F` Wilson coefficient.

Status: Conjectural. A mere enlargement of the static primitive family does not count unless it produces a response relation that cannot be represented by finite static coordinates.

## What Counts

Status: Counterexample candidate. A phase-lagged quadrature response counts against an instantaneous static sensitivity basis when `Omega tau_chi` is finite and the quadrature is tied to the same parameters that determine the in-phase response.

Status: Counterexample candidate. A frequency-dependent transfer function counts when observations at more than one drive frequency require the same relaxation time `tau_chi` and cannot be matched by a finite static coefficient vector.

Status: Counterexample candidate. Sidebands count if the response contains frequencies absent from the linear input drive, such as `Omega_1 + Omega_2` or `|Omega_1 - Omega_2|`.

Status: Counterexample candidate. Hysteresis or delayed threshold response counts if it cannot be represented by analytic local coefficients around a single equilibrium.

## MVP Classification Test

Status: Proven. For the first-order model under `F(t)=F0 cos(Omega t)`, write

```text
chi(t) = A cos(Omega t) + B sin(Omega t)
```

Then

```text
A = alpha F0 / (1 + Omega^2 tau_chi^2)
B = alpha F0 Omega tau_chi / (1 + Omega^2 tau_chi^2)
```

Status: Counterexample candidate. `B != 0` is the first quadrature candidate.

Status: Proven. `B != 0` is not enough by itself against a comparator that already contains a free `dot F` Wilson coefficient.

Status: Counterexample candidate. The stronger target is the relation

```text
B / A = Omega tau_chi
```

across drive frequencies, together with the amplitude suppression

```text
sqrt(A^2 + B^2) / F0 = |alpha| / sqrt(1 + Omega^2 tau_chi^2).
```

## Basis Audit Update

Status: Proven. Against an instantaneous static `F` coefficient, the in-phase
component is degenerate and the finite-`Omega tau_chi` quadrature remains
outside the static basis.

Status: Proven. Against a one-frequency basis `{F, dot F}` with unconstrained
coefficients, the same quadrature is degenerate because `dot F` supplies the
missing sine quadrature at that one frequency.

Status: Counterexample candidate. The surviving target is therefore not a
single quadrature point by itself. It is the shared frequency-dependent
transfer relation `alpha c_chi/(1 + i Omega tau_chi)`, which cannot be
represented exactly by a finite shared-coefficient polynomial derivative basis
over varying drive frequency.

## Frequency-Sweep Target

Status: Proven. A degree-`N` complex derivative polynomial can absorb `N + 1`
frequency samples by interpolation if its coefficients are freely fit to that
sample set.

Status: Counterexample candidate. The first exact frequency-sweep observable
requires at least `floor((N+1)/2)+1` distinct positive frequencies with one
shared real coefficient vector.

Status: Counterexample candidate. If the comparator is granted complex
coefficients, the conservative exact target is `N + 2` distinct frequencies.

Status: Proven. The exact obstruction is the nonzero interpolation residual
derived in [`frequency-sweep-distinguishability.md`](frequency-sweep-distinguishability.md).

## Concrete Projection Target

Status: Counterexample candidate. A normalized orbital invariant
`F(Y(t))=(a/r(t))^p` supplies a single-system two-frequency target at `n` and
`2n` through `O(e^2)`.

Status: Proven. With a calibrated linear projection
`O_hat_k = Lambda G(i Omega_k) F_hat_k`, the observable recovers the same
transfer pair after dividing by `Lambda F_hat_k`.

Status: Proven. If `Lambda(Omega_k)` is an arbitrary complex nuisance at every
frequency, the observable target collapses because the nuisance can absorb the
pole relation.

## Projection Channel Target

Status: Counterexample candidate. An acceleration-like channel
`delta a_hat = Gamma q_hat` preserves the pole relation up to one shared scale.

Status: Counterexample candidate. A range-like channel
`ddot R + kappa^2 R = Gamma q_A` has
`Lambda_R(z)=Gamma/(kappa^2+z^2)`, a finite rational projection rather than an
arbitrary frequency-local nuisance.

Status: Proven. If `Gamma` and `kappa` are known or calibrated, the range-like
channel deprojects exactly to `G(z)` away from projection poles.

Status: Counterexample candidate. If `Gamma` and `kappa` are unknown but shared,
the projection remains finite-dimensional; exact distinguishability then needs
independent calibration or additional frequency samples.

## Triple Shared-Tau Carrier Bridge

Status: Counterexample candidate. A hierarchical triple can provide two
existing GR carrier frequencies, `Omega_in` and `Omega_out`, without requiring
a new tensor harmonic.

Status: Counterexample candidate. If the two frequencies are distinct and both
carrier amplitudes are nonzero, they are enough to test a real
shared-coefficient degree-`1` or degree-`2` derivative EFT.

Status: Proven. The same two carriers are not enough against real degree
`N >= 3` comparators or against an overpowered complex degree-`1` polynomial
comparator. Those require additional carrier frequencies or a justified
derivative-order prior.

## Triple GR Carrier Inventory

Status: Counterexample candidate. The existing GR outer-dipole combination
carrier `|Omega_in-Omega_out|` can provide the third sample needed to strengthen
the shared-`tau_chi` test.

Status: Counterexample candidate. With the three-carrier set
`Omega_in`, `Omega_out`, and `|Omega_in-Omega_out|`, a nonresonant triple can
pressure real shared-coefficient derivative comparators through degree
`N <= 4` and the deliberately generous complex comparator through degree
`N <= 1`.

Status: Proven. This three-carrier target collapses if the combination carrier
is not distinct, has zero projection, or is absorbed by an arbitrary complex
projection nuisance assigned independently to each carrier.

## Triple Projection-Nuisance Gate

Status: Counterexample candidate. The three-carrier bridge remains a positive
route if the triple projection factors have calibrated or finite shared
structure:

```text
O_k = Lambda_k(theta) G(i Omega_k) F_k.
```

Status: Counterexample candidate. A calibrated or common real projection scale
preserves the shared-`tau_chi` transfer relation because it can be deprojected
or removed by carrier ratios.

Status: Counterexample candidate. Finite real or complex shared projection
nuisance is a conditional gate: it does not erase the bridge point by point,
but it consumes rank and requires projection priors, additional carriers, or a
Jacobian audit.

Status: Proven. Arbitrary per-carrier complex projection nuisance collapses
the target because each carrier can choose
`Lambda_k=O_k/(G(i Omega_k)F_k)`.

## Physical Projection-Manifold Gate

Status: Counterexample candidate. The runtime-worthy dynamic-chi target is the
phase-locked triple projection manifold

```text
Lambda_c = A_c exp(i(phi_in - phi_out)).
```

Status: Proven. This manifold has generic rank `5` in the six-real-dimensional
three-carrier complex-amplitude vector, so it is finite and shared rather than
pointwise arbitrary.

Status: Proven. If the combination carrier phase is floated independently as
`A_c exp(i phi_c)`, the projection rank becomes `6` and the bridge collapses
to arbitrary per-carrier complex nuisance.

Status: Conjectural. External timing runtime is justified only if the named
timing model preserves the phase/geometric links instead of fitting
independent complex carrier amplitudes.

## Named Timing-Model Projection Audit

Status: Counterexample candidate. The public `Nutimo` triple timing core is
source-inspected as a finite fitted-parameter, state-vector, and geometry
projection model. That class is compatible with the runtime-worthy side of
the Request 10.5 gate.

Status: Proven. An explicit sinusoid special case with per-harmonic amplitudes
and phases is also present in the source tree. If that nuisance is enabled on
the target carrier frequencies, the dynamic-chi bridge collapses before a
runtime fit can be scientifically informative.

Status: Counterexample candidate. The named-code observable target is
therefore not a generic `J0337` detectability claim. It is a conditional
hand-off target: test a dynamic-chi column against the standard finite
three-body timing Jacobian while excluding arbitrary harmonic soak-up terms.

## External Runtime-Worthiness Pilot

Status: Counterexample candidate. Runtime is now justified only as a narrow
gate sequence:

```text
configuration closure
-> named Jacobian rank gate
-> minimal synthetic shared-tau_chi injection
```

Status: Proven. A runtime configuration with `RN_PL`-like per-harmonic carrier
soak-up enabled is not an observable test of dynamic-chi; it is the collapse
comparator already identified in Request 10.6.

Status: Counterexample candidate. The first runtime observable target is a
finite Jacobian question, not a posterior question: whether the named
three-carrier projection stays below rank `6/6` and whether a dynamic-chi
test column remains outside the standard fitted-parameter span.

## Sideband Test

Status: Proven. A linear time-invariant one-state model driven by `F1 cos(Omega1 t) + F2 cos(Omega2 t)` cannot create sidebands in the linear readout.

Status: Proven. A minimal nonlinear internal drive
`tau_chi dot chi + chi = alpha F + beta_F2 F^2` generates
`2Omega1`, `2Omega2`, `Omega1+Omega2`, and `|Omega1-Omega2|`.

Status: Proven. A minimal nonlinear readout
`lambda_Fchi F chi + lambda_chi2 chi^2` generates sum and difference
sidebands even when the internal equation remains linear.

Status: Counterexample candidate. For the M3 orbital forcing with `n` and `2n`,
the nonlinear drive creates a `3n` sideband with amplitude
`beta_F2 p^2(p+3)e^3/4`, which is absent from the linear `O(e^2)` input.

Status: Proven. A linear time-invariant projection cannot synthesize a
frequency absent from its input; it only multiplies an already generated
sideband by `Lambda(i nu)`.
