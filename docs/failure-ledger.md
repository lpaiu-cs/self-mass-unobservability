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
| F7 | Proven | The electric-only exact current primitive set is not yet adequate as a physically justified minimal free-fall sector. | A physically motivated magnetic tidal family `B_ij` is excluded only by working restriction, not by a derived theorem assumption. | [`primitive-set-adequacy.md`](primitive-set-adequacy.md), [`magnetic-family-ordering.md`](magnetic-family-ordering.md), and [`../lemmas/10-magnetic-family-obstruction.md`](../lemmas/10-magnetic-family-obstruction.md) now record `B2` as the smallest explicit adequacy obstruction. |
| F8 | Proven | The corrected exact-current-set `E/B` theorem is not yet adequate as a physically justified minimal free-fall sector either. | A parity-even scalar-like external family is excluded only by working restriction, not by a derived theorem assumption. | [`scalar-family-ordering.md`](scalar-family-ordering.md), [`../lemmas/13-scalar-family-obstruction.md`](../lemmas/13-scalar-family-obstruction.md), and [`../symbolic/es_sector_delta4.py`](../symbolic/es_sector_delta4.py) now record `S` as the smallest explicit scalar-family obstruction. |
| F9 | Proven | Even a derivative-only or shift-symmetric scalar family does not rescue minimal-sector uniqueness. | Removing bare `S` alone is insufficient; derivative invariants still survive unless a stronger suppression principle is added. | [`../lemmas/14-derivative-only-scalar-audit.md`](../lemmas/14-derivative-only-scalar-audit.md) and [`../symbolic/shift_scalar_sector_delta4.py`](../symbolic/shift_scalar_sector_delta4.py) now record the first weight-`4` derivative-only scalar obstructions, with canonical representative `dotS2`. |

## Current Theorem Failure Boundary

- Status: Proven. Abstract monomial finiteness is no longer the unresolved step once a finite primitive catalog with positive weights is fixed.
- Status: Proven. The previous minimal obstruction to the five-element target is explicit: `divE2`, with `mixedGradE2` as a second surviving gradient invariant.
- Status: Proven. For the exact current primitive set, the corrected seven-element normal-form basis is now linearly independent.
- Status: Proven. The live bottleneck for electric-only salvage remains magnetic-family ordering.
- Status: Proven. The raw `E/B` survivor candidate list had one explicit mixed quartic dependence relation.
- Status: Proven. After that correction, the `E/B`-expanded `Delta<=4` sector closes on a finite linearly independent `18`-element basis.
- Status: Proven. The smallest explicit adequacy obstruction found so far is the magnetic-family survivor `B2`.
- Status: Proven. The scalar-family extension contributes the smaller broad-program adequacy obstruction `S` as soon as a rank-0 external family is admitted.
- Status: Proven. The corrected `E/B+scalar` `Delta<=4` sector closes on a finite linearly independent `33`-element basis.
- Status: Proven. The derivative-only scalar audit removes bare `S` but still produces new weight-`4` survivors, including canonical obstruction `dotS2`.
- Status: Proven. Minimal-sector uniqueness is therefore already effectively dead across the audited unsuppressed family extensions.
- Status: Conjectural. The positive finite-family collapse candidate remains plausible because fixed-order closure survives in the corrected `E/B+scalar` audit and in the derivative-only scalar audit.
- Status: Conjectural. The live positive bottleneck is no longer uniqueness; it is the finite-family collapse theorem beyond the explicitly audited catalogs.
- Status: Counterexample candidate. If the sector is enlarged beyond the explicit audited family catalogs by another physically admissible family, the first new survivor in that extension reinforces the uniqueness no-go and only obstructs the positive branch if fixed-order finiteness itself fails.
