# Lemma 61: Markovianizable Versus Genuine Hereditary Memory

## Statement

- Status: Proven. Work within the one-coordinate causal convolution class

```math
{\cal M}[Y](\tau) = \int_{-\infty}^{\tau} d\tau'\,K(\tau-\tau')\,Y(\tau').
```

- Status: Proven. Finite sums of exponential kernels are finitely Markovianizable and therefore are not sharp pure `A3` escapes.
- Status: Proven. The causal power-law kernel `K_{\gamma}` with `0<\gamma<1` is not finitely Markovianizable inside that class.
- Status: Proven. Therefore the smallest explicit genuine hereditary class in the current program is the one-coordinate power-law kernel branch, not the exponential branch.

## Exponential Control

1. Status: Proven. Define

```math
\chi(\tau)
=
\int_{-\infty}^{\tau} d\tau' \,
\frac{e^{-(\tau-\tau')/T_h}}{T_h}\,Y(\tau').
```

2. Status: Proven. Differentiating under the integral gives

```math
\dot{\chi}(\tau) = -\frac{1}{T_h}\chi(\tau) + \frac{1}{T_h}Y(\tau).
```

3. Status: Proven. Hence the exponential-memory monopole model is equivalent to a local worldline action with one auxiliary state variable `\chi`.
4. Status: Proven. More generally, any finite linear combination of exponentials yields a finite-dimensional first-order state system.
5. Status: Proven. Therefore the exponential branch is not a pure `A3`-only escape. It collapses back into an `A4`-type finite-state extension.

## Genuine Hereditary Branch

1. Status: Proven. A finite-dimensional linear state-space realization with constant coefficients yields a rational Laplace-domain transfer function in the Laplace variable `s`.
2. Status: Proven. The exponential kernel has Laplace transform

```math
\widehat{K}_{\exp}(s)=\frac{1}{1+sT_h},
```

which is rational.
3. Status: Proven. The power-law kernel has Laplace transform

```math
\widehat{K}_{\gamma}(s) = (T_h s)^{\gamma-1},
\qquad 0<\gamma<1,
```

up to the fixed normalization chosen in Lemma 60.
4. Status: Proven. The function `s^{\gamma-1}` is not rational and carries a branch cut.
5. Status: Proven. Therefore no finite-dimensional constant-coefficient local state-space realization exists for the power-law kernel inside the one-coordinate causal convolution class used here.

## Verdict

- Status: Proven. The exponential-memory control is finitely Markovianizable and therefore should not be counted as the sharp `A3` counterexample.
- Status: Proven. The one-coordinate power-law kernel is the smallest explicit genuine hereditary counterexample class identified so far.
