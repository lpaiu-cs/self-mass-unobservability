# Nonanalytic Counterexample Program

## Target

- Status: Proven. The target of this boundary-stress program is no longer family-envelope completeness or fixed-order operator finiteness.
- Status: Proven. The target is the smallest local nonanalytic monopole model that keeps the finite family and operator envelope intact but breaks analytic monopole jet collapse.
- Status: Proven. The program is confined to the parity-even, nonspinning, local MVP free-fall theorem domain, with `A5` the only stress point.

## Layer Separation

- Status: Proven. Finite primitive-family envelope closure concerns which irreducible scalar, vector, and STF primitive families are admitted.
- Status: Proven. Finite operator-space closure concerns how many parity-even local scalar operators exist at fixed order `\Delta \le 4`.
- Status: Proven. Analytic monopole jet collapse concerns whether the local monopole response `m_A(Y)` is captured by finitely many Taylor coefficients in the corrected scalar normal-form coordinates `Y^I`.
- Status: Proven. Wilson-data bookkeeping concerns the residual higher-multipole coefficient set and remains separate from monopole-response data.

## Canonical Boundary Test

- Status: Proven. The sharpest local `A5` stress test uses a single corrected scalar normal-form coordinate `Y` and keeps all primitive-family and reduction data fixed.
- Status: Proven. The canonical smallest explicit candidate is the smooth flat one-coordinate monopole model

```math
m_A(Y) = m_0 + \alpha \,\phi_{\mathrm{flat}}(Y),
\qquad
\phi_{\mathrm{flat}}(Y) =
\begin{cases}
0, & Y \le 0,\\
e^{-1/Y^2}, & Y > 0.
\end{cases}
```

- Status: Proven. This model is local, parity-even, and one-coordinate, but it is nonanalytic at the reference point `Y = 0`.
- Status: Proven. The threshold comparison model

```math
m_A(Y) = m_0 + \alpha\,\Theta(Y-Y_c)\sqrt{Y-Y_c}
```

is also local and finite-family, but it is rougher because it adds an explicit activation point and a branch singularity.

## Current Verdict

- Status: Proven. The smallest explicit local nonanalytic counterexample to the analytic jet step is the one-coordinate smooth flat monopole model.
- Status: Proven. It keeps the irreducible family envelope, fixed admitted family content, fixed-order operator finiteness, and finite normal-form quotient intact.
- Status: Proven. It breaks only the analytic monopole jet step of [`../lemmas/55-monopole-jet-collapse.md`](../lemmas/55-monopole-jet-collapse.md).
- Status: Proven. The replacement bookkeeping is non-Taylor monopole germ data, not Wilson coefficients.
