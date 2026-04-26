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
