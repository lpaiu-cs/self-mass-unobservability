# Lemma 60: Hereditary Kernel Construction

## Statement

- Status: Proven. Keep the irreducible scalar/vector/STF primitive-family envelope fixed and keep the reduced scalar normal-form basis `Y^I` already established in the local theorem domain.
- Status: Proven. Restrict attention to a single corrected scalar coordinate `Y(\tau)` and modify only the monopole response law.
- Status: Proven. Then the following three kernel models are explicit and sufficient for the first `A3` boundary-stress program:
  a local analytic control,
  an exponential-memory control,
  and a power-law hereditary candidate.

## Explicit Models

- Status: Proven. The local analytic control is

```math
m_A(Y(\tau)) = m_0 + \alpha Y(\tau) + \beta Y(\tau)^2.
```

- Status: Proven. The exponential-memory control is

```math
m_A[Y](\tau)
=
m_0 + \alpha Y(\tau)
 + \lambda_{\exp}\int_{-\infty}^{\tau} d\tau' \,
\frac{e^{-(\tau-\tau')/T_h}}{T_h}\,Y(\tau').
```

- Status: Proven. The genuine hereditary candidate is

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

## Why These Models Respect The Scope

1. Status: Proven. The primitive-family content is unchanged, because all three models use the same already admitted reduced scalar coordinate `Y` and add no new primitive tensor family.
2. Status: Proven. The fixed-order local operator side is unchanged, because the hereditary term is introduced only at the monopole-response layer after the local normal-form basis has already been fixed.
3. Status: Proven. The local analytic control leaves both `A3` and `A5` active and therefore serves as the baseline theorem-domain case.
4. Status: Proven. The exponential control drops locality in kernel form but keeps analyticity on the instantaneous side; it is therefore a clean control case for separating nonlocality from nonanalyticity.
5. Status: Proven. The power-law kernel also keeps the instantaneous analytic term `\alpha Y(\tau)` and changes only the history dependence; it is therefore the smallest one-coordinate candidate aimed specifically at `A3`.

## Boundary

- Status: Proven. This lemma constructs the kernel classes only; it does not yet classify which of them are reducible to finitely many local states.
- Status: Proven. That classification is the next step of the `A3` boundary-stress program.
