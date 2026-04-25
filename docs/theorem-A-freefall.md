# Theorem A Candidate: Free-Fall Fixed-Order Collapse

## Main Positive Target: Finite-Family Fixed-Order Collapse Candidate

- Status: Conjectural. Work in the free-fall MVP sector only: nonspinning, nearly spherical, parity-even, local worldline EFT, fixed truncation order `\Delta \le 4`, and finite admitted external field content.
- Status: Conjectural. Assume the body is quasi-static, self-bound in equilibrium, and carries no orbital-timescale internal state variable in the free-fall sector.
- Status: Conjectural. Assume couplings to the admitted external-family catalog are analytic once the local operator space has been reduced to normal form at fixed order.
- Status: Imported from prior work. Exclude literal internal self-unobservability as an equilibrium mechanism.
- Status: Imported from prior work. Exclude COM self-subtraction as a source of a new monopole force.
- Status: Conjectural. Then, for an explicit finite admitted primitive-family catalog, the admissible local free-fall scalar operator space should collapse to a finite family-conditioned normal form modulo the explicitly stated total-derivative, lower-order-EOM, and algebraic reduction rules.
- Status: Conjectural. Consequently, by the relevant conditional collapse lemma, the leading body-dependent free-fall deviation should collapse to finitely many sensitivity coordinates plus finitely many higher-multipole Wilson coefficients for that admitted family catalog.

## Exact-Current-Set Specializations

### Electric-Only Exact-Current-Set Theorem

- Status: Conjectural. For the exact current electric-only primitive set at `\Delta \le 4`, the admissible local free-fall scalar operator space closes on

```math
\{E2,\ E3,\ E2^2,\ dotE2,\ gradE2,\ divE2,\ mixedGradE2\}.
```

- Status: Proven. This corrected seven-element electric-only basis is linearly independent.
- Status: Conjectural. This is an exact-current-set theorem candidate, not yet a physically justified minimal-sector theorem.

### Corrected `E/B` Exact-Current-Set Theorem

- Status: Proven. The raw `E/B` survivor candidate set is not linearly independent; one explicit mixed quartic dependence relation reduces it.
- Status: Conjectural. After that correction, the exact-current-set `E/B` basis at `\Delta \le 4` is

```math
\{E2,\ B2,\ E3,\ EB2,\ E2^2,\ B2^2,\ dotE2,\ dotB2,\ EBDtB,\ E2B2,\ EB\_sq,\ TrE2B2,\ gradE2,\ divE2,\ mixedGradE2,\ gradB2,\ divB2,\ mixedGradB2\}.
```

- Status: Proven. This corrected `E/B` basis is linearly independent.
- Status: Conjectural. This is likewise an exact-current-set theorem candidate, not yet a physically justified minimal-sector theorem.

## Negative Branch: Family-Admission No-Go For Minimal-Sector Uniqueness

- Status: Proven. Minimal-sector uniqueness is not the positive theorem target anymore.
- Status: Proven. Across the explicitly audited family classes, unsuppressed family admission always yields a new low-order survivor.
- Status: Proven. For STF rank-2 family admission, the self witness `X2` and the mixed witness `EX` both appear at weight `2`; the audited explicit self instance is `B2`.
- Status: Proven. For unsuppressed rank-0 scalar-family admission, the self witness `S` appears before the first mixed witness `SE2`.
- Status: Proven. For derivative-only or shift-symmetric scalar-family admission, the self witness `dotS2` and the mixed witness `DtS_E2` both first appear at weight `4`.
- Status: Proven. For a genuine parity-even rank-1 vector-family admission, the self witness `V2` appears at weight `2` and the first mixed witness `EVV` appears at weight `3`.
- Status: Proven. For a genuine parity-even rank-3 STF family admission, the self witness `T2` appears at weight `2` and the first mixed witness `ETT` appears at weight `3`.
- Status: Proven. For a genuine parity-even rank-4 STF family admission, the self witness `Q2` appears at weight `2` and the first mixed witness `EQQ` appears at weight `3`.
- Status: Proven. Therefore a stronger physically justified minimal-sector theorem is already on the no-go branch unless explicit suppression, ordering, or background-restriction assumptions are added family by family.

## Scope Separation

- Status: Proven. The exact-current-set electric-only theorem and the corrected exact-current-set `E/B` theorem are special audited instances of the broader positive finite-family target.
- Status: Proven. The stronger physically justified minimal-sector theorem is a separate claim and is not rescued by current audits.
- Status: Proven. The negative branch is the family-admission no-go for minimal-sector uniqueness.
- Status: Conjectural. The positive branch is the finite-family fixed-order collapse candidate.

## Proof Route

1. Status: Proven. Use the fixed-order counting rule in [`power-counting.md`](power-counting.md) to define a bounded-weight operator search space.
2. Status: Proven. Finite admitted family content and positive weights imply a finite set of decorated primitive building blocks at order `\Delta \le 4`.
3. Status: Proven. [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) exhaustively enumerates all parity-even scalar contractions from the exact current electric-only primitive set.
4. Status: Proven. [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) reduces all total-derivative, lower-order-EOM, and single-family Cayley-Hamilton-reducible electric classes explicitly.
5. Status: Proven. [`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py) shows that the corrected seven electric-only survivors are linearly independent as operators.
6. Status: Proven. [`../symbolic/eb_sector_delta4.py`](../symbolic/eb_sector_delta4.py) shows that admitting the magnetic family yields a finite enlarged `E/B` candidate list rather than a fixed-order blowup.
7. Status: Proven. [`../symbolic/eb_survivor_rank_check.py`](../symbolic/eb_survivor_rank_check.py) extracts the first explicit `E/B` dependence relation and the corrected linearly independent `18`-element basis.
8. Status: Proven. [`../symbolic/es_sector_delta4.py`](../symbolic/es_sector_delta4.py) shows that adjoining one scalar-like external family still yields a finite enlarged `E/B+scalar` candidate list rather than a fixed-order blowup.
9. Status: Proven. [`../symbolic/es_survivor_rank_check.py`](../symbolic/es_survivor_rank_check.py) shows that the corrected `E/B+scalar` `33`-element basis is linearly independent.
10. Status: Proven. [`../symbolic/shift_scalar_sector_delta4.py`](../symbolic/shift_scalar_sector_delta4.py) shows that even derivative-only scalar admission still yields new weight-`4` survivors beyond the corrected `E/B` sector.
11. Status: Proven. [`../symbolic/family_witness_map.py`](../symbolic/family_witness_map.py) records the current family-class witness table and the exact theorem layer each audited class obstructs.
12. Status: Proven. [`family-admission-theorem.md`](family-admission-theorem.md) promotes the repeated family-attack pattern into an explicit no-go statement for the audited family classes.
13. Status: Proven. [`suppression-budget-theorem.md`](suppression-budget-theorem.md), [`../lemmas/17-witness-threshold-classification.md`](../lemmas/17-witness-threshold-classification.md), and [`../symbolic/witness_threshold_map.py`](../symbolic/witness_threshold_map.py) promote the audited witnesses into necessary suppression thresholds for minimal-sector uniqueness at `\Delta_{\max} = 4`.
14. Status: Proven. [`mixed-witness-thresholds.md`](mixed-witness-thresholds.md), [`../lemmas/18-rank2-mixed-witness-threshold.md`](../lemmas/18-rank2-mixed-witness-threshold.md), [`../lemmas/19-rank0-mixed-witness-threshold.md`](../lemmas/19-rank0-mixed-witness-threshold.md), and [`../symbolic/mixed_witness_map.py`](../symbolic/mixed_witness_map.py) classify the first self and first mixed witnesses for the audited family classes.
15. Status: Proven. [`threshold-formulas.md`](threshold-formulas.md), [`../lemmas/20-rank2-threshold-formula.md`](../lemmas/20-rank2-threshold-formula.md), [`../lemmas/21-rank0-threshold-formula.md`](../lemmas/21-rank0-threshold-formula.md), and [`../symbolic/threshold_formula_check.py`](../symbolic/threshold_formula_check.py) lift the audited family classes into formula-level threshold statements.
16. Status: Proven. [`budget-composition-theorem.md`](budget-composition-theorem.md), [`../lemmas/22-pairwise-composition-audit.md`](../lemmas/22-pairwise-composition-audit.md), [`../lemmas/23-three-family-composition-audit.md`](../lemmas/23-three-family-composition-audit.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) first closed the old audited-set composition layer for `R2`, `R0a`, and `R0b`.
17. Status: Proven. [`vector-family-ordering.md`](vector-family-ordering.md), [`../lemmas/26-rank1-family-admission.md`](../lemmas/26-rank1-family-admission.md), [`../lemmas/27-rank1-threshold-formula.md`](../lemmas/27-rank1-threshold-formula.md), [`../symbolic/r1_sector_delta4.py`](../symbolic/r1_sector_delta4.py), and [`../symbolic/r1_survivor_rank_check.py`](../symbolic/r1_survivor_rank_check.py) show that a genuine local parity-even rank-1 vector family is a real obstruction class rather than a derivative-generated artifact.
18. Status: Proven. [`../lemmas/28-r1-pairwise-composition-audit.md`](../lemmas/28-r1-pairwise-composition-audit.md), [`../lemmas/29-r1-augmented-composition-audit.md`](../lemmas/29-r1-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, and quadruple composition layer for the enlarged audited family set `{R2, R0a, R0b, R1}`.
19. Status: Proven. [`family-envelope-theorem.md`](family-envelope-theorem.md), [`family-envelope-table.md`](family-envelope-table.md), [`../lemmas/24-family-envelope-census.md`](../lemmas/24-family-envelope-census.md), [`../lemmas/25-next-unaudited-family-gate.md`](../lemmas/25-next-unaudited-family-gate.md), and [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py) separate audited-set sufficiency from envelope completeness and identify the odd-rank tensor gate beyond the resolved vector gate.
20. Status: Proven. [`rank3-family-ordering.md`](rank3-family-ordering.md), [`../lemmas/30-rank3-family-admission.md`](../lemmas/30-rank3-family-admission.md), [`../lemmas/31-rank3-threshold-formula.md`](../lemmas/31-rank3-threshold-formula.md), [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py), and [`../symbolic/r3_survivor_rank_check.py`](../symbolic/r3_survivor_rank_check.py) now show that `Rodd+` is a genuine new obstruction class rather than a derivative-generated artifact.
21. Status: Proven. [`../lemmas/32-r3-pairwise-composition-audit.md`](../lemmas/32-r3-pairwise-composition-audit.md), [`../lemmas/33-r3-augmented-composition-audit.md`](../lemmas/33-r3-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, and five-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}`.
22. Status: Proven. [`rank4-family-ordering.md`](rank4-family-ordering.md), [`../lemmas/34-rank4-family-admission.md`](../lemmas/34-rank4-family-admission.md), [`../lemmas/35-rank4-threshold-formula.md`](../lemmas/35-rank4-threshold-formula.md), [`../symbolic/r4_sector_delta4.py`](../symbolic/r4_sector_delta4.py), and [`../symbolic/r4_survivor_rank_check.py`](../symbolic/r4_survivor_rank_check.py) now show that `Reven4+` is a genuine new obstruction class rather than a trace-descended or derivative-generated artifact.
23. Status: Conjectural. The remaining broader positive burden is no longer family-envelope census alone; it is a fresh enlarged audited-set composition re-close that includes the new rank-4 threshold, and only after that the next envelope gate can move beyond `Reven4+`.

## Current Verdict

- Status: Proven. The old five-scalar electric target was false because it omitted `divE2` and `mixedGradE2`.
- Status: Proven. The corrected seven-scalar electric-only basis is linearly independent for the exact current primitive set.
- Status: Proven. The corrected `E/B` basis is linearly independent for the corrected exact-current-set `E/B` primitive sector.
- Status: Proven. The stronger physically justified minimal-sector theorem should not presently be treated as alive.
- Status: Proven. Minimal-sector uniqueness is already effectively dead across the audited unsuppressed family classes.
- Status: Proven. No audited family class is currently harmless without extra assumptions.
- Status: Proven. The vector-family clarification shows that derivative-generated vector blocks do not absorb a genuine primitive family `V_i`.
- Status: Proven. The smallest explicit vector-family obstruction is `V2`, with first mixed witness `EVV`.
- Status: Proven. The rank-3 clarification shows that derivative-generated rank-3 blocks do not absorb a genuine primitive STF family `T_{ijk}`.
- Status: Proven. The smallest explicit rank-3 obstruction is `T2`, with first mixed witness `ETT`.
- Status: Proven. The rank-4 clarification shows that trace descendants and derivative-generated rank-4 blocks do not absorb a genuine primitive STF family `Q_{ijkl}`.
- Status: Proven. The smallest explicit rank-4 obstruction is `Q2`, with first mixed witness `EQQ`.
- Status: Proven. Necessary suppression budgets are explicit for the audited classes: for rank-2 STF families the self-only lower bound is `w_X \ge 3` but the mixed-aware necessary threshold is `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness; for bare scalar families the self-only sharp threshold is `w_S \ge 5`; for derivative-only scalar families the tied-sharp threshold is `w_D \ge 3`; for genuine vector families the self-only sharp threshold is `w_V \ge 3`; for genuine rank-3 STF families the self-only sharp threshold is `w_T \ge 3`; for genuine rank-4 STF families the self-only sharp threshold is `w_Q \ge 3`.
- Status: Proven. The pairwise, triple, and quadruple composition audits already found no surviving operator beyond the baseline electric sector once the current audited thresholds were imposed across the enlarged audited family set `{R2, R0a, R0b, R1}`.
- Status: Proven. The post-Rodd+ pairwise, triple, quadruple, and five-family composition audits likewise find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}`.
- Status: Proven. Therefore the current thresholds are jointly sufficient for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}` at `\Delta \le 4`, and no explicit post-Rodd+ cross-family obstruction has appeared in the audited set.
- Status: Proven. The raw R3-extended survivor list is not linearly independent; its first exact quartic mixed dependence relation is `-E2T2/2 + 2 E2T2_mixed_1 - E2T2_mixed_2 + E2T2_mixed_3 = 0`.
- Status: Proven. The raw R3-extended `21`-element survivor list is bookkeeping only and should not be read as a corrected basis statement; the composition audit uses the thresholded rank-3 family instead.
- Status: Proven. The raw R4-extended survivor list is not linearly independent; its first exact quartic mixed dependence relation is `-E2Q2/2 + 2 E2Q2_mixed_1 - E2Q2_mixed_2 + E2Q2_mixed_3 = 0`.
- Status: Proven. The raw R4-extended `22`-element survivor list has rank `19` and nullity `3`, and it is bookkeeping only rather than a corrected basis statement.
- Status: Proven. The current MVP assumptions still do not close the family envelope.
- Status: Proven. The family-envelope census now points to `Rodd5+` as the next smallest unaudited class, but the theorem does not move there yet because the `Reven4+` composition re-close is still pending.
- Status: Proven. The current live bottleneck is the `Reven4+` obstruction: a genuine parity-even rank-4 STF family yields smallest witness `Q2`, so the enlarged audited set must be re-closed before the theorem can move to the next envelope gate.
- Status: Conjectural. The broad positive target remains finite-family fixed-order collapse, and it remains plausible because the audited enlarged family catalogs, including the current rank-3 extension, still close finitely at `\Delta \le 4`.

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
- Status: Proven. [`../lemmas/13-scalar-family-obstruction.md`](../lemmas/13-scalar-family-obstruction.md)
- Status: Proven. [`../lemmas/14-derivative-only-scalar-audit.md`](../lemmas/14-derivative-only-scalar-audit.md)
- Status: Proven. [`../lemmas/15-rank2-family-admission.md`](../lemmas/15-rank2-family-admission.md)
- Status: Proven. [`../lemmas/16-rank0-family-admission.md`](../lemmas/16-rank0-family-admission.md)
- Status: Proven. [`../lemmas/18-rank2-mixed-witness-threshold.md`](../lemmas/18-rank2-mixed-witness-threshold.md)
- Status: Proven. [`../lemmas/19-rank0-mixed-witness-threshold.md`](../lemmas/19-rank0-mixed-witness-threshold.md)
- Status: Proven. [`../lemmas/20-rank2-threshold-formula.md`](../lemmas/20-rank2-threshold-formula.md)
- Status: Proven. [`../lemmas/21-rank0-threshold-formula.md`](../lemmas/21-rank0-threshold-formula.md)
- Status: Proven. [`../lemmas/22-pairwise-composition-audit.md`](../lemmas/22-pairwise-composition-audit.md)
- Status: Proven. [`../lemmas/23-three-family-composition-audit.md`](../lemmas/23-three-family-composition-audit.md)
- Status: Proven. [`../lemmas/24-family-envelope-census.md`](../lemmas/24-family-envelope-census.md)
- Status: Proven. [`../lemmas/25-next-unaudited-family-gate.md`](../lemmas/25-next-unaudited-family-gate.md)
- Status: Proven. [`../lemmas/26-rank1-family-admission.md`](../lemmas/26-rank1-family-admission.md)
- Status: Proven. [`../lemmas/27-rank1-threshold-formula.md`](../lemmas/27-rank1-threshold-formula.md)
- Status: Proven. [`../lemmas/28-r1-pairwise-composition-audit.md`](../lemmas/28-r1-pairwise-composition-audit.md)
- Status: Proven. [`../lemmas/29-r1-augmented-composition-audit.md`](../lemmas/29-r1-augmented-composition-audit.md)
- Status: Proven. [`../lemmas/30-rank3-family-admission.md`](../lemmas/30-rank3-family-admission.md)
- Status: Proven. [`../lemmas/31-rank3-threshold-formula.md`](../lemmas/31-rank3-threshold-formula.md)
- Status: Conjectural. [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md)
- Status: Conjectural. [`eb-conditional-collapse.md`](eb-conditional-collapse.md)
- Status: Conjectural. [`power-counting.md`](power-counting.md)
- Status: Proven. [`family-admission-no-go.md`](family-admission-no-go.md)
- Status: Proven. [`family-admission-theorem.md`](family-admission-theorem.md)
- Status: Proven. [`family-envelope-theorem.md`](family-envelope-theorem.md)
- Status: Proven. [`family-envelope-table.md`](family-envelope-table.md)
- Status: Proven. [`rank3-family-ordering.md`](rank3-family-ordering.md)
- Status: Proven. [`vector-family-ordering.md`](vector-family-ordering.md)
- Status: Proven. [`mixed-witness-thresholds.md`](mixed-witness-thresholds.md)
- Status: Proven. [`threshold-formulas.md`](threshold-formulas.md)
- Status: Proven. [`budget-composition-theorem.md`](budget-composition-theorem.md)
- Status: Proven. [`composition-status.md`](composition-status.md)
- Status: Proven. [`suppression-budget-theorem.md`](suppression-budget-theorem.md)
- Status: Proven. [`sharp-threshold-status.md`](sharp-threshold-status.md)
- Status: Proven. [`broad-collapse-reformulation.md`](broad-collapse-reformulation.md)
- Status: Proven. [`family-class-table.md`](family-class-table.md)
- Status: Proven. [`primitive-catalog.md`](primitive-catalog.md)
- Status: Proven. [`primitive-set-adequacy.md`](primitive-set-adequacy.md)
- Status: Proven. [`magnetic-family-ordering.md`](magnetic-family-ordering.md)
- Status: Proven. [`scalar-family-ordering.md`](scalar-family-ordering.md)
- Status: Proven. [`reduction-rules.md`](reduction-rules.md)
- Status: Proven. [`../symbolic/witness_threshold_map.py`](../symbolic/witness_threshold_map.py)
- Status: Proven. [`../symbolic/mixed_witness_map.py`](../symbolic/mixed_witness_map.py)
- Status: Proven. [`../symbolic/threshold_formula_check.py`](../symbolic/threshold_formula_check.py)
- Status: Proven. [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py)
- Status: Proven. [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py)
- Status: Proven. [`../symbolic/r1_sector_delta4.py`](../symbolic/r1_sector_delta4.py)
- Status: Proven. [`../symbolic/r1_survivor_rank_check.py`](../symbolic/r1_survivor_rank_check.py)
- Status: Proven. [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py)
- Status: Proven. [`../symbolic/r3_survivor_rank_check.py`](../symbolic/r3_survivor_rank_check.py)

## Failure Triggers

- Status: Counterexample candidate. An explicit `chi_A` state invalidates the no-state reduction step.
- Status: Counterexample candidate. Nonanalytic activation invalidates the analytic collapse step.
- Status: Counterexample candidate. Hereditary couplings invalidate locality before basis closure is even applied.
- Status: Counterexample candidate. Infinite or uncontrolled external field content invalidates the fixed-order closure theorem candidate.
- Status: Counterexample candidate. Any additional `\Delta \le 4` scalar contraction from the exact current electric-only primitive set that evades the explicit enumeration/reduction pipeline obstructs the electric-only exact-current-set theorem.
- Status: Counterexample candidate. Any further physically admissible primitive family that yields new survivors beyond the corrected explicit audited family catalogs reinforces the uniqueness no-go branch.
- Status: Counterexample candidate. Any further physically admissible primitive family obstructs the positive finite-family collapse branch only if it also destroys fixed-order finiteness or the family-conditioned reduction path.
