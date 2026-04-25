# Second-Order Projection Audit

Status: Counterexample candidate. This note lifts the one-state projection
audit to the second-order internal mode.

## Internal Transfer

Status: Imported from prior work. The second-order internal transfer is

```math
H_2(z)=
\frac{\alpha\omega_\chi^2}
{\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2}.
```

Define

```math
D_\chi(z)=\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2.
```

## Acceleration-Like Readout

Status: Proven. For an acceleration-like channel,

```math
\hat O_a(z)=\Gamma H_2(z)\hat F(z).
```

The finite-shared-nuisance diagnostic is

```math
\frac{\hat O_a(z)}{\hat F(z)}D_\chi(z)
-\Gamma\alpha\omega_\chi^2=0.
```

Status: Proven. If `Gamma` is known or fitted as one shared nuisance
parameter, the internal quadratic denominator survives the projection.

## Range-Like Readout

Status: Proven. For a range-like channel,

```math
(\partial_t^2+\kappa^2)R=\Gamma q,
\qquad
\hat O_R(z)=
\frac{\Gamma H_2(z)}{\kappa^2+z^2}\hat F(z).
```

The finite-shared-nuisance diagnostic is

```math
\frac{\hat O_R(z)}{\hat F(z)}
(\kappa^2+z^2)D_\chi(z)
-\Gamma\alpha\omega_\chi^2=0.
```

Status: Proven. If `Gamma` and `kappa^2` are known or fitted as finite shared
nuisance parameters, the internal resonance is not erased by the range
projection.

## Collapse Boundaries

Status: Proven. The projection audit collapses if `Gamma=0`.

Status: Proven. It also collapses as an identifiability claim if each
frequency is assigned an independent complex projection nuisance.

Status: Proven. The range channel has an additional projection singularity at
`\kappa^2+z^2=0`; samples exactly on this projection pole are not valid
deprojection points.

Status: Counterexample candidate. A deliberately pole-cancelling projection
could hide the internal denominator, but that is not produced by the
acceleration-like or range-like channels above.

## Verdict

Status: Counterexample candidate. Under finite shared projection nuisance, the
second-order resonance/phase-wrapping target survives projection. The next
burden is not projection survival; it is whether the resonant line shape beats
a static shared comparator budget.

