# Failure Ledger

This ledger records exact failure modes, not vague concerns. If a proof step breaks, the failed step and the smallest missing assumption should appear here.

| ID | Status | Failed step or loophole | Minimal violated or missing assumption | Current action |
| --- | --- | --- | --- | --- |
| F1 | Conjectural | The M3 theorem candidate still lacks a proof that every admissible local parity-even `Delta<=4` operator reduces to the target normal forms. | `normal-form completeness modulo total derivatives and lower-order equations of motion` for the exact `Delta<=4` primitive catalog. | Keep this burden in [`theorem-A-freefall.md`](theorem-A-freefall.md) and [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md), not in the assumptions ledger. |
| F2 | Counterexample candidate | A light internal coordinate `chi_A(t)` with orbital-timescale dynamics can enter the free-fall action as an explicit worldline field. | A4 fails. | See [`counterexamples/chi-state/README.md`](../counterexamples/chi-state/README.md). |
| F3 | Counterexample candidate | A threshold response such as `sqrt(Y - Y_c) Theta(Y - Y_c)` breaks the analytic sensitivity jet. | A5 fails. | See [`counterexamples/nonanalytic-activation/README.md`](../counterexamples/nonanalytic-activation/README.md). |
| F4 | Counterexample candidate | Retarded memory kernels can preserve body dependence that is not instantaneous in the local invariants. | A3 fails. | See [`counterexamples/hereditary/README.md`](../counterexamples/hereditary/README.md). |
| F5 | Counterexample candidate | A nonmetric clock coupling can evade the free-fall theorem because it lives outside the MVP observable sector. | The free-fall-only scope assumption fails. | See [`counterexamples/nonmetric-clock/README.md`](../counterexamples/nonmetric-clock/README.md). |
| F6 | Counterexample candidate | If finite external field content fails, infinitely many primitive species can survive at the same fixed order and block basis closure before the conditional collapse lemma applies. | A8 fails. | Keep this loophole explicit rather than smuggling a finite basis into the assumptions. |

## Current Theorem Failure Boundary

- Status: Proven. Abstract monomial finiteness is no longer the unresolved step once a finite primitive catalog with positive weights is fixed.
- Status: Conjectural. The exact remaining burden is `normal-form completeness modulo total derivatives and lower-order equations of motion` for the exact `Delta<=4` primitive catalog.
- Status: Counterexample candidate. If that reduction fails, the repository should extract the smallest explicit loophole model rather than restoring finite-basis closure as an assumption.
