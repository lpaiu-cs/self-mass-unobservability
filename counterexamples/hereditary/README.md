# hereditary loophole

## Target

- Status: Proven. This counterexample class now targets the sharp `A3` boundary after the positive finite-family collapse theorem and the local nonanalytic `A5` boundary are both already explicit.
- Status: Proven. The correct target is not any retarded kernel whatsoever, but the smallest genuinely hereditary kernel that cannot be collapsed into finitely many local auxiliary states.

## Kernel Split

- Status: Proven. The local analytic control keeps the theorem domain untouched:

```math
m_A(Y(\tau)) = m_0 + \alpha Y(\tau) + \beta Y(\tau)^2.
```

- Status: Proven. The exponential-memory control

```math
m_A[Y](\tau)
=
m_0 + \alpha Y(\tau)
+ \lambda_{\exp}\int_{-\infty}^{\tau} d\tau' \,
\frac{e^{-(\tau-\tau')/T_h}}{T_h}\,Y(\tau')
```

is nonlocal in kernel form but finitely Markovianizable.

- Status: Proven. The smallest explicit genuine hereditary branch is the one-coordinate causal power-law kernel

```math
m_A[Y](\tau)
=
m_0 + \alpha Y(\tau)
+ \lambda_{\gamma}\int_{-\infty}^{\tau} d\tau' \,
K_{\gamma}(\tau-\tau')\,Y(\tau'),
```

with

```math
K_{\gamma}(s)
=
\frac{\Theta(s)}{\Gamma(1-\gamma)\,T_h}
\left(\frac{s}{T_h}\right)^{-\gamma},
\qquad 0<\gamma<1.
```

## Why The Power-Law Model Is Sharp

- Status: Proven. The exponential control is not a pure `A3` success because it can be rewritten with one auxiliary local state `\chi_A`.
- Status: Proven. The power-law kernel keeps the primitive-family envelope and fixed-order local operator closure intact, because it changes only the monopole-response functional on the already admitted coordinate `Y`.
- Status: Proven. It breaks the reduction to a local monopole function `m_A(Y(\tau))` of finitely many instantaneous normal-form coordinates.
- Status: Proven. Therefore the exact theorem layer broken is the local monopole reduction prerequisite to Lemmas 55 and 56, not fixed-order operator finiteness.
