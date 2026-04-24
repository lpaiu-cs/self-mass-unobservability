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

## Sideband Test

Status: Proven. A linear time-invariant one-state model driven by `F1 cos(Omega1 t) + F2 cos(Omega2 t)` cannot create sidebands in the linear readout.

Status: Conjectural. Sideband novelty requires a nonlinear readout, nonlinear drive construction, or a non-linear internal equation.
