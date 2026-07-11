# Hereditary Counterexample Program

## Target

- Status: Proven. The target of this boundary-stress program is no longer family-envelope completeness, fixed-order operator finiteness, or analytic monopole collapse.
- Status: Proven. The target is the smallest genuinely hereditary monopole model that keeps the finite primitive-family envelope and fixed-order local operator closure intact but breaks the local finite-family collapse theorem.
- Status: Proven. The program is confined to the parity-even, nonspinning, free-fall MVP theorem domain, with locality `A3` the only intended stress point and instantaneous analyticity kept whenever possible.

## Layer Separation

- Status: Proven. Local analytic collapse is the theorem already closed in [`finite-family-collapse-theorem.md`](finite-family-collapse-theorem.md).
- Status: Proven. The local nonanalytic `A5` escape is the separate one-coordinate smooth-flat counterexample recorded in [`nonanalytic-counterexample-program.md`](nonanalytic-counterexample-program.md).
- Status: Proven. Finite-state Markovianizable memory is not yet a sharp `A3` escape, because it can be re-read as a finite auxiliary-state extension of the local theory.
- Status: Proven. Genuinely hereditary memory is the class in which the monopole response depends on a past-history functional that is not reducible to finitely many local state variables.

## Kernel Classes

- Status: Proven. The local analytic control keeps the already proved theorem domain untouched:

```math
m_A(Y(\tau)) = m_0 + \alpha Y(\tau) + \beta Y(\tau)^2.
```

- Status: Proven. The exponential-memory control keeps the primitive-family and operator content fixed but writes the monopole response as

```math
m_A[Y](\tau)
=
m_0 + \alpha Y(\tau)
+ \lambda_{\exp}\int_{-\infty}^{\tau} d\tau' \,
\frac{e^{-(\tau-\tau')/T_h}}{T_h}\,Y(\tau').
```

- Status: Proven. The genuine hereditary candidate keeps the same one-coordinate scalar content but uses the causal power-law kernel

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

## Current Verdict

- Status: Proven. The exponential kernel is not the sharp `A3` boundary because it is finitely Markovianizable and therefore collapses into a one-state `A4`-type extension.
- Status: Proven. The smallest genuine hereditary counterexample found so far is the one-coordinate causal power-law kernel model.
- Status: Proven. It preserves the finite scalar/vector/STF primitive-family envelope, fixed-order local operator finiteness, and the finite local normal-form quotient.
- Status: Proven. It breaks the reduction of the monopole response to a local function `m_A(Y(\tau))` of finitely many instantaneous normal-form coordinates.
- Status: Proven. The replacement bookkeeping is kernel or spectral data, not Wilson coefficients.
