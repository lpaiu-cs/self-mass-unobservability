# Shared-Tau Observable Targets

Status: Counterexample candidate. This note defines the observable targets for
the shared-tau theorem/no-go track.

## Primary Targets

Status: Counterexample candidate. The two-tone target is

```math
R_{\Omega_1,\Omega_2}
=
\frac{\hat O(\Omega_1+\Omega_2)}
{\hat O(\Omega_1)\hat O(\Omega_2)}
```

after deprojection and forcing normalization.

Status: Counterexample candidate. The orbital target is

```math
R_{3n}
=
\frac{\hat O(3n)}{\hat O(n)\hat O(2n)}
```

when eccentric forcing supplies the input harmonics `n` and `2n`.

Status: Proven. These ratios are not claimed to be dynamic merely because the
numerator exists. Static nonlinear local terms can create the numerator.

## What Counts

Status: Counterexample candidate. A positive M7 observable requires all of the
following:

- the linear lines fit one `G(z)=c_Y+\beta/(1+\tau_\chi z)`;
- the generated line fits `C_side/(1+\tau_\chi z_side)`;
- the same `tau_chi` appears in both fits;
- the number of distinct samples exceeds the finite static comparator
  interpolation budget.

Status: Proven. A single ratio point does not count if an unconstrained complex
nonlinear coefficient is allowed for that point.

Status: Proven. A finite set of ratio points does not count if the static
comparator has enough independent polynomial coefficients to interpolate all
of them.

## What Counts As Failure

Status: Proven. Failure mode 1 is static mimicry: the observed line and
sideband data can be fit by finite `P_N(z)` and `Q_M(u,v)` within their
coefficient budgets.

Status: Proven. Failure mode 2 is nuisance collapse: projection or nonlinear
readout nuisance is allowed to vary independently at every frequency or pair.

Status: Proven. Failure mode 3 is dynamical collapse: `tau_chi=0`, `C_side=0`,
or the sampled band satisfies only an adiabatic derivative expansion to the
claimed accuracy.

## Output Rule

Status: Counterexample candidate. A successful M7 result should report the
minimum sample count, comparator class, surviving ratio law, and collapse
condition together. Otherwise it is only another sideband existence statement.

