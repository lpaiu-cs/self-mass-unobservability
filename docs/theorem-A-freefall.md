# Theorem A Candidate: Free-Fall Fixed-Order Collapse

## Candidate 1: Electric-Only Exact-Current-Set Theorem

- Status: Conjectural. Work in the free-fall MVP sector only: nonspinning, nearly spherical, parity-even, local worldline EFT, fixed truncation order `\Delta \le 4`, and finite external field content.
- Status: Conjectural. Assume the body is quasi-static, self-bound in equilibrium, and carries no orbital-timescale internal state variable in the free-fall sector.
- Status: Conjectural. Assume couplings to the external environment are analytic once the local operator space has been reduced to normal form at fixed order.
- Status: Imported from prior work. Exclude literal internal self-unobservability as an equilibrium mechanism.
- Status: Imported from prior work. Exclude COM self-subtraction as a source of a new monopole force.
- Status: Conjectural. Then, for the exact current electric-only primitive set at `\Delta \le 4`, the admissible local free-fall scalar operator space closes on

```math
\{E2,\ E3,\ E2^2,\ dotE2,\ gradE2,\ divE2,\ mixedGradE2\}
```

modulo total derivatives, lower-order equations of motion, and the explicitly stated algebraic identities.
- Status: Conjectural. Consequently, by [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md), the leading body-dependent free-fall deviation collapses to finitely many sensitivity coordinates plus finitely many higher-multipole Wilson coefficients.

## Candidate 2: Stronger Physically Justified Minimal-Sector Theorem

- Status: Proven. This stronger claim is not currently available.
- Status: Proven. It would require the electric-only primitive set to be justified as the physically correct minimal free-fall sector rather than merely the exact current primitive set.
- Status: Proven. The magnetic-family survivor `B2` blocks that upgrade unless an explicit ordering or background-restriction assumption is added.
- Status: Proven. [`magnetic-family-ordering.md`](magnetic-family-ordering.md) records the allowed ways to make that suppression explicit rather than tacit.

## Candidate 3: Broader Finite-Dimensional Collapse Conjecture

- Status: Proven. The raw `E/B` survivor candidate set is not linearly independent; [`../lemmas/11-eb-survivor-independence-delta4.md`](../lemmas/11-eb-survivor-independence-delta4.md) finds one explicit mixed quartic dependence relation.
- Status: Conjectural. After reducing by that explicit relation, the corrected `E/B`-expanded `\Delta \le 4` basis is

```math
\{E2,\ B2,\ E3,\ EB2,\ E2^2,\ B2^2,\ dotE2,\ dotB2,\ EBDtB,\ E2B2,\ EB\_sq,\ TrE2B2,\ gradE2,\ divE2,\ mixedGradE2,\ gradB2,\ divB2,\ mixedGradB2\}.
```

- Status: Proven. Therefore the magnetic family does not presently falsify finite fixed-order closure by itself.
- Status: Proven. The corrected `E/B` basis is linearly independent.
- Status: Conjectural. [`eb-conditional-collapse.md`](eb-conditional-collapse.md) records the conditional finite-dimensional collapse step for this corrected `E/B` basis.
- Status: Conjectural. The broader finite-dimensional collapse program remains plausible after the magnetic extension, although its exact reduced basis is larger than the electric-only basis.

## Scope Separation

- Status: Proven. Candidate 1 is a theorem candidate for the exact current electric-only primitive set.
- Status: Proven. Candidate 2 is the stronger physically justified minimal-sector theorem and is blocked at M6.
- Status: Conjectural. Candidate 3 is the broader finite-dimensional collapse conjecture after explicit `E/B` enlargement and mixed-quartic correction.

## Proof Route

1. Status: Proven. Use the fixed-order counting rule in [`power-counting.md`](power-counting.md) to define a bounded-weight operator search space.
2. Status: Proven. Finite external field content and positive weights imply a finite set of decorated primitive building blocks at order `\Delta \le 4`.
3. Status: Proven. [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) exhaustively enumerates all parity-even scalar contractions from the exact current electric-only primitive set.
4. Status: Proven. [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) reduces all total-derivative, lower-order-EOM, and single-family Cayley-Hamilton-reducible electric classes explicitly.
5. Status: Proven. [`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py) shows that the corrected seven electric-only survivors are linearly independent as operators.
6. Status: Proven. [`../symbolic/eb_sector_delta4.py`](../symbolic/eb_sector_delta4.py) shows that admitting the magnetic family yields a finite enlarged `E/B` candidate list rather than a fixed-order blowup.
7. Status: Proven. [`../symbolic/eb_survivor_rank_check.py`](../symbolic/eb_survivor_rank_check.py) extracts the first explicit `E/B` dependence relation and the corrected linearly independent `18`-element basis.
8. Status: Conjectural. Invoke the conditional collapse lemma only after choosing whether the theorem target is the exact current electric-only set or the corrected enlarged `E/B` sector.

## Current Verdict

- Status: Proven. The old five-scalar electric target was false because it omitted `divE2` and `mixedGradE2`.
- Status: Proven. The corrected seven-scalar electric-only basis is linearly independent for the exact current primitive set.
- Status: Proven. The live electric-only salvage bottleneck is magnetic-family ordering, not survivor dependence.
- Status: Proven. The electric-only theorem candidate remains valid only as an exact-current-set statement.
- Status: Proven. The stronger physically justified minimal-sector theorem should not presently be treated as established.
- Status: Proven. The raw `E/B` candidate basis needed one explicit mixed quartic correction.
- Status: Conjectural. The broader finite-dimensional collapse conjecture remains plausible after the magnetic extension because the corrected `E/B` sector is still finite and linearly independent at `\Delta \le 4`.

## Dependencies

- Status: Imported from prior work. [`../lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md)
- Status: Imported from prior work. [`../lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md)
- Status: Conjectural. [`../lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md)
- Status: Proven. [`../lemmas/05-finite-basis-closure.md`](../lemmas/05-finite-basis-closure.md)
- Status: Conjectural. [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md)
- Status: Conjectural. [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md)
- Status: Proven. [`../lemmas/08-mixed-time-derivative-audit.md`](../lemmas/08-mixed-time-derivative-audit.md)
- Status: Proven. [`../lemmas/09-survivor-independence-delta4.md`](../lemmas/09-survivor-independence-delta4.md)
- Status: Proven. [`../lemmas/10-magnetic-family-obstruction.md`](../lemmas/10-magnetic-family-obstruction.md)
- Status: Proven. [`../lemmas/11-eb-survivor-independence-delta4.md`](../lemmas/11-eb-survivor-independence-delta4.md)
- Status: Proven. [`../lemmas/12-magnetic-ordering-salvage.md`](../lemmas/12-magnetic-ordering-salvage.md)
- Status: Conjectural. [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md)
- Status: Conjectural. [`eb-conditional-collapse.md`](eb-conditional-collapse.md)
- Status: Conjectural. [`power-counting.md`](power-counting.md)
- Status: Proven. [`primitive-catalog.md`](primitive-catalog.md)
- Status: Proven. [`primitive-set-adequacy.md`](primitive-set-adequacy.md)
- Status: Proven. [`magnetic-family-ordering.md`](magnetic-family-ordering.md)
- Status: Proven. [`reduction-rules.md`](reduction-rules.md)

## Failure Triggers

- Status: Counterexample candidate. An explicit `chi_A` state invalidates the no-state reduction step.
- Status: Counterexample candidate. Nonanalytic activation invalidates the analytic collapse step.
- Status: Counterexample candidate. Hereditary couplings invalidate locality before basis closure is even applied.
- Status: Counterexample candidate. Infinite or uncontrolled external field content invalidates the fixed-order closure theorem candidate.
- Status: Counterexample candidate. Any additional `\Delta \le 4` scalar contraction from the exact current electric-only primitive set that evades the explicit enumeration/reduction pipeline would obstruct Candidate 1.
- Status: Counterexample candidate. Any further physically admissible primitive family that yields new survivors beyond the corrected explicit `E/B` basis would enlarge Candidate 3 or obstruct its present finite-basis path.
