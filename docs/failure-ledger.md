# Failure Ledger

This ledger records exact failure modes, not vague concerns. If a proof step breaks, the failed step and the smallest missing assumption should appear here.

| ID | Status | Failed step or loophole | Minimal violated or missing assumption | Current action |
| --- | --- | --- | --- | --- |
| F1 | Proven | The previous five-element `Delta<=4` target omitted explicit surviving gradient operators. | The omitted operators are `divE2` and `mixedGradE2`. | The corrected normal-form path now includes them explicitly in [`theorem-A-freefall.md`](theorem-A-freefall.md), [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md), and [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py). |
| F2 | Counterexample candidate | A light internal coordinate `chi_A(t)` with orbital-timescale dynamics can enter the free-fall action as an explicit worldline field. | A4 fails. | See [`counterexamples/chi-state/README.md`](../counterexamples/chi-state/README.md). |
| F3 | Counterexample candidate | A threshold response such as `sqrt(Y - Y_c) Theta(Y - Y_c)` breaks the analytic sensitivity jet. | A5 fails. | See [`counterexamples/nonanalytic-activation/README.md`](../counterexamples/nonanalytic-activation/README.md). |
| F4 | Counterexample candidate | Retarded memory kernels can preserve body dependence that is not instantaneous in the local invariants. | A3 fails. | See [`counterexamples/hereditary/README.md`](../counterexamples/hereditary/README.md). |
| F5 | Counterexample candidate | A nonmetric clock coupling can evade the free-fall theorem because it lives outside the MVP observable sector. | The free-fall-only scope assumption fails. | See [`counterexamples/nonmetric-clock/README.md`](../counterexamples/nonmetric-clock/README.md). |
| F6 | Counterexample candidate | If finite external field content fails, infinitely many primitive species can survive at the same fixed order and block basis closure before the conditional collapse lemma applies. | A8 fails. | Keep this loophole explicit rather than smuggling a finite basis into the assumptions. |

## Current Theorem Failure Boundary

- Status: Proven. Abstract monomial finiteness is no longer the unresolved step once a finite primitive catalog with positive weights is fixed.
- Status: Proven. The previous minimal obstruction to the five-element target is explicit: `divE2`, with `mixedGradE2` as a second surviving gradient invariant.
- Status: Conjectural. For the exact current primitive set, the present completeness path is the corrected seven-element normal-form basis rather than the older five-element target.
- Status: Counterexample candidate. If the explicit contraction generator misses a valid `Delta<=4` scalar class, that omitted class would become the next minimal obstruction.
