# chi-state loophole

## Minimal Action

- Status: Counterexample candidate. Introduce one explicit orbital-timescale internal coordinate `\chi_A(\tau)` with inertial scale `\mu_\chi`:

```math
S_A
=
\int d\tau \left[
-m_0 \left(1 + s_1 Y(z_A) + \frac12 s_2 Y(z_A)^2\right)
+ \frac{\mu_\chi}{2} \dot{\chi}_A^2
- \frac{\mu_\chi \omega_\chi^2}{2} \chi_A^2
+ g_\chi \chi_A Y(z_A)
\right].
```

- Status: Counterexample candidate. The action remains local and analytic in the fields kept explicitly.
- Status: Counterexample candidate. The only directly violated theorem assumption is A4: no orbital-timescale internal state variable.

## Explicit Dynamics

- Status: Counterexample candidate. The Euler-Lagrange equation is

```math
\ddot{\chi}_A + \omega_\chi^2 \chi_A = \frac{g_\chi}{\mu_\chi} Y(z_A).
```

- Status: Counterexample candidate. The formal solution is

```math
\chi_A
=
\chi_{A,\mathrm{hom}}
+ \frac{g_\chi}{\mu_\chi}\left(\omega_\chi^2 + D_\tau^2\right)^{-1} Y.
```

- Status: Counterexample candidate. The homogeneous piece `\chi_{A,\mathrm{hom}}` carries independent initial data.
- Status: Counterexample candidate. The free-fall response therefore depends on more than the instantaneous local invariant `Y`; it depends on the extra state and its initial conditions.

## Why This Escapes Instantaneous Sensitivity Collapse

- Status: Counterexample candidate. When `\omega_\chi \sim \Omega_{\rm orb}`, the resolvent `(\omega_\chi^2 + D_\tau^2)^{-1}` is not reducible to a short local derivative expansion over the orbital timescale.
- Status: Counterexample candidate. In that regime the effective monopole response is not an algebraic function of finitely many instantaneous sensitivity coordinates.
- Status: Counterexample candidate. This is a genuine loophole to the theorem candidate if A4 is false, even though locality and analyticity remain intact.

## Adiabatic-Collapse Limit

- Status: Counterexample candidate. If `\omega_\chi \gg \Omega_{\rm orb}` and the homogeneous mode is absent or decays, then

```math
\chi_A
=
\frac{g_\chi}{\mu_\chi \omega_\chi^2}
\left[
1 - \frac{D_\tau^2}{\omega_\chi^2}
+ O\!\left(\frac{D_\tau^4}{\omega_\chi^4}\right)
\right] Y.
```

- Status: Counterexample candidate. Substituting back gives a local expansion of the form

```math
\Delta L_{\rm eff}
=
\frac{g_\chi^2}{2\mu_\chi \omega_\chi^2} Y^2
+ \frac{g_\chi^2}{2\mu_\chi \omega_\chi^4} (\dot Y)^2
+ O\!\left(\omega_\chi^{-6}\right),
```

up to total derivatives.

- Status: Counterexample candidate. In this heavy-state limit the loophole collapses back into ordinary local sensitivity and Wilson-coefficient data.
- Status: Counterexample candidate. The `chi` loophole is therefore minimal in a precise sense: it fails only when the extra state survives on the orbital timescale.

## Smallest-Model Reason

- Status: Counterexample candidate. Only one new degree of freedom is added, and the rest of the worldline EFT remains local and analytic.
