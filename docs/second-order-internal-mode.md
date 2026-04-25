# Second-Order Internal Mode

Status: Counterexample candidate. This note starts the next dynamic loophole
after the one-state relaxation program reached an observational sample-budget
ceiling for the current forcing dictionaries.

## Model

Status: Counterexample candidate. The minimal second-order internal state is

```math
\mu_\chi \ddot\chi
+\gamma_\chi \dot\chi
+\omega_\chi^2(\chi-\alpha F)=0 .
```

Equivalently, for derivative variable `z`,

```math
\chi(z)=H_2(z)F(z),
\qquad
H_2(z)=
\frac{\alpha\omega_\chi^2}
{\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2}.
```

Status: Counterexample candidate. Compared with the one-state relaxation
transfer, `H_2` has a quadratic denominator and can carry resonance and phase
wrapping rather than a single relaxation rolloff.

## Monochromatic Response

Status: Proven. For `F(t)=F_0 cos(Omega t)`,

```math
H_2(i\Omega)
=
\frac{\alpha\omega_\chi^2}
{\omega_\chi^2-\mu_\chi\Omega^2+i\gamma_\chi\Omega}.
```

Writing `\chi=A cos(Omega t)+B sin(Omega t)`, the normalized components are

```math
A/F_0=
\frac{\alpha\omega_\chi^2(\omega_\chi^2-\mu_\chi\Omega^2)}
{(\omega_\chi^2-\mu_\chi\Omega^2)^2+\gamma_\chi^2\Omega^2},
```

```math
B/F_0=
\frac{\alpha\omega_\chi^2\gamma_\chi\Omega}
{(\omega_\chi^2-\mu_\chi\Omega^2)^2+\gamma_\chi^2\Omega^2}.
```

Status: Proven. The phase tangent is

```math
\frac{B}{A}
=
\frac{\gamma_\chi\Omega}
{\omega_\chi^2-\mu_\chi\Omega^2}.
```

This expression changes sign through the resonance denominator and is the
first place where phase wrapping can appear.

## Resonance

Status: Proven. The amplitude denominator is minimized at

```math
\Omega_{\rm peak}^2
=
\frac{\omega_\chi^2}{\mu_\chi}
-\frac{\gamma_\chi^2}{2\mu_\chi^2},
```

provided

```math
\gamma_\chi^2 < 2\mu_\chi\omega_\chi^2 .
```

Status: Counterexample candidate. This resonance/line-shape structure is the
main reason the second-order mode is more promising than the first-order
relaxation state for budget-breaking observables.

## Collapse Limits

Status: Proven. If `mu_chi=0`, then

```math
H_2(z)
=
\frac{\alpha}{1+(\gamma_\chi/\omega_\chi^2)z},
```

so the model collapses to first-order relaxation with
`\tau_\chi=\gamma_\chi/\omega_\chi^2`.

Status: Proven. In the low-frequency regime,

```math
H_2(z)
=
\alpha
-\alpha\frac{\gamma_\chi}{\omega_\chi^2}z
+\alpha\left(
\frac{\gamma_\chi^2}{\omega_\chi^4}
-\frac{\mu_\chi}{\omega_\chi^2}
\right)z^2
+O(z^3).
```

Thus, over a sufficiently small band, the second-order mode again collapses
order-by-order into local derivative coefficients.

Status: Counterexample candidate. Non-collapse requires sampling near enough
to the internal mode that the quadratic denominator, resonance, or phase
wrapping cannot be replaced by the finite local derivative comparator within
the chosen sample budget.

## Next Observable Target

Status: Counterexample candidate. The next theorem target is a
sample-budgeted line-shape test: determine the minimum number of frequency
samples needed to distinguish `H_2(z)` from a finite static polynomial
comparator with projection nuisance included.

