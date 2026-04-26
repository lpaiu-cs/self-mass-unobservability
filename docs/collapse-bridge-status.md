# Collapse Bridge Status

## Purpose

- Status: Proven. This note records the bridge from irreducible family-envelope closure to the positive finite-family collapse theorem.
- Status: Proven. It is written to keep four distinct layers separate:
  envelope closure,
  operator finiteness,
  monopole jet collapse,
  and sensitivity versus Wilson bookkeeping.

| Layer | Status | Exact dependency | Current outcome |
| --- | --- | --- | --- |
| Irreducible family-envelope closure | Proven | [`irreducible-envelope-theorem.md`](irreducible-envelope-theorem.md), [`../lemmas/50-cartesian-irrep-reduction.md`](../lemmas/50-cartesian-irrep-reduction.md), [`../lemmas/51-trace-descendant-absorption.md`](../lemmas/51-trace-descendant-absorption.md), [`../lemmas/52-mixed-symmetry-family-gate.md`](../lemmas/52-mixed-symmetry-family-gate.md) | Not the bottleneck anymore. |
| Fixed-order candidate operator finiteness | Proven | [`../lemmas/53-fixed-order-operator-finiteness.md`](../lemmas/53-fixed-order-operator-finiteness.md), [`../symbolic/fixed_family_operator_count.py`](../symbolic/fixed_family_operator_count.py), [`power-counting.md`](power-counting.md), `A7`, `A8` | Closed positively inside the current theorem domain. |
| Finite normal-form basis closure | Proven | [`../lemmas/54-normal-form-basis-closure.md`](../lemmas/54-normal-form-basis-closure.md), [`reduction-rules.md`](reduction-rules.md), [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) | Closed positively as a finite-dimensional quotient statement. |
| Analytic monopole jet collapse | Proven | [`../lemmas/55-monopole-jet-collapse.md`](../lemmas/55-monopole-jet-collapse.md), `A5` | Closed positively inside the analytic theorem domain; fails exactly if `A5` is dropped. |
| Sensitivity versus Wilson split | Proven | [`../lemmas/56-sensitivity-vs-wilson-split.md`](../lemmas/56-sensitivity-vs-wilson-split.md), `A3`, `A5`, `A8` | Closed positively inside the local analytic theorem domain; fails exactly if locality or analyticity is dropped. |
| Overall positive finite-family collapse theorem | Proven | [`finite-family-collapse-theorem.md`](finite-family-collapse-theorem.md), [`theorem-A-freefall.md`](theorem-A-freefall.md) | Closed positively inside the current theorem domain. |

## Readout

- Status: Proven. Envelope closure is no longer the bridge bottleneck.
- Status: Proven. Fixed-order operator finiteness and finite-dimensional scalar normal-form closure are now proven.
- Status: Proven. Monopole jet collapse is not being smuggled in; it is closed explicitly through analyticity `A5`.
- Status: Proven. Sensitivity coefficients and Wilson coefficients are no longer conflated.
- Status: Proven. The remaining risks are explicit assumption-drop failures, not a hidden unresolved bridge step inside the present theorem domain.

## First Explicit Assumption-Drop Stress Test

- Status: Proven. The former boundary-stress bottleneck `"construct or refute the smallest local nonanalytic counterexample to analytic monopole jet collapse (A5 boundary)"` is now resolved positively by the one-coordinate smooth flat monopole model `m_A(Y)=m_0+\alpha e^{-1/Y^2}\Theta(Y)`.
- Status: Proven. This stress test leaves irreducible family-envelope closure, finite admitted family content, fixed-order operator finiteness, and finite normal-form closure intact.
- Status: Proven. It breaks only the analytic Taylor-jet step of [`../lemmas/55-monopole-jet-collapse.md`](../lemmas/55-monopole-jet-collapse.md).
- Status: Proven. The replacement bookkeeping is non-Taylor monopole germ data, not Wilson coefficients.

## A3 Boundary-Stress: Hereditary Counterexample Program

- Status: Proven. The former boundary-stress bottleneck `"construct or refute the smallest genuinely hereditary counterexample to the local finite-family collapse theorem (A3 boundary)"` is now resolved positively by the one-coordinate causal power-law kernel model of [`hereditary-counterexample-program.md`](hereditary-counterexample-program.md).
- Status: Proven. The exponential-memory control case is not the sharp `A3` escape because it is finitely Markovianizable and collapses into a finite auxiliary-state `A4`-type extension.
- Status: Proven. The genuine power-law kernel keeps the finite primitive-family envelope, fixed admitted family content, fixed-order local operator finiteness, and finite scalar normal-form quotient intact.
- Status: Proven. The first theorem layer it breaks is the reduction of the monopole response to a local function `m_A(Y(\tau))` of finitely many instantaneous normal-form coordinates, so Lemmas 55 and 56 are no longer applicable in their local form.
- Status: Proven. The replacement bookkeeping is kernel or spectral memory data, not Wilson coefficients.

## A4 Boundary-Stress: Finite Local State-Variable Counterexample Program

- Status: Proven. The former boundary-stress bottleneck `"construct or refute the smallest local finite-state counterexample to the no-state theorem (A4 boundary), and determine whether a finite state-augmented collapse theorem survives"` is now resolved positively by the one-state local analytic model of [`state-variable-counterexample-program.md`](state-variable-counterexample-program.md):

```math
m_A(Y,\chi)=m_0+\alpha Y+\lambda\chi,
\qquad
\dot{\chi}=-(1/T_h)\chi+(\beta/T_h)Y.
```

- Status: Proven. The adiabatic or slaved control case is not the sharp `A4` escape, because when `\chi` is algebraically or adiabatically eliminable the model collapses back into a `Y`-only monopole description.
- Status: Proven. The genuinely dynamical one-state branch keeps locality `A3`, analyticity `A5`, finite primitive-family content, fixed-order operator finiteness, and the finite scalar normal-form quotient intact.
- Status: Proven. The first theorem layer it breaks is the original no-state reduction to a monopole function `m_A(Y)` of instantaneous scalar normal-form coordinates alone; equivalently, Lemma 55 and the Y-only reading of Lemma 56 fail in that branch.
- Status: Proven. A second positive branch survives: once finitely many local state variables `\chi^a` are kept explicitly, a finite state-augmented collapse theorem survives as a finite-dimensional local state-space statement.
- Status: Proven. The replacement bookkeeping is finite augmented state-space data `(Y^I,\chi^a)` plus local state-evolution parameters, and it remains separate from higher-multipole Wilson coefficients.
