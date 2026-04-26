# Theorem A Candidate: Free-Fall Fixed-Order Collapse

## Main Positive Target: Finite-Family Fixed-Order Collapse Theorem

- Status: Proven. Work in the free-fall MVP sector only: nonspinning, nearly spherical, parity-even, local worldline EFT, fixed truncation order `\Delta \le 4`, and finite admitted external field content.
- Status: Proven. Within that theorem domain, the irreducible primitive-family envelope is already closed on the audited scalar/vector/STF classes.
- Status: Proven. For any explicit finite admitted primitive-family catalog drawn from that closed envelope, the candidate parity-even local free-fall scalar operator space at `\Delta \le 4` is finite before reduction.
- Status: Proven. Modulo the explicit family-envelope, total-derivative, lower-order-EOM, algebraic, and extracted linear-dependence relations already fixed in this repo, the reduced scalar operator space is finite-dimensional and admits a finite normal-form basis `Y^I`.
- Status: Proven. Under the locality and analyticity assumptions `A3` and `A5`, the monopole term `m_A(Y)` therefore has a finite fixed-order Taylor jet in finitely many scalar normal-form coordinates.
- Status: Proven. The residual higher-multipole sector is then carried by finitely many Wilson coefficients `C_{A,\alpha}` that are separate from the monopole sensitivity coefficients.
- Status: Proven. Therefore the positive finite-family collapse theorem closes inside the current theorem domain.

## Positive Bridge Layers

- Status: Proven. Layer 1 is irreducible family-envelope closure; this is already closed positively on the audited scalar/vector/STF classes.
- Status: Proven. Layer 2 is fixed-order candidate operator finiteness; this now follows from finite admitted family content, positive weights, and fixed `\Delta \le 4`.
- Status: Proven. Layer 3 is finite normal-form basis closure; this now follows because the reduced operator space is a quotient of a finite-dimensional candidate space.
- Status: Proven. Layer 4 is analytic monopole jet collapse; this now follows conditionally on `A5` once the finite scalar basis `Y^I` is fixed.
- Status: Proven. Layer 5 is the sensitivity/Wilson split; this now follows conditionally on `A3`, `A5`, and `A8`.

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
- Status: Proven. Rank `2` remains a separate mixed-aware STF special case because the mixed quadratic witness `EX` exists.
- Status: Proven. For every genuinely new parity-even STF primitive family with rank `L \ge 3`, there is no linear scalar witness, no mixed quadratic `E Y_L` witness, and the quadratic norm `Y2` always exists.
- Status: Proven. Therefore the higher-rank STF branch obeys a universal self-only sharp threshold `w_Y \ge 3` at `\Delta \le 4`.
- Status: Proven. The stronger mixed-pattern theorem still fails at audited rank `L = 4`, because the cubic mixed class `EEQ` survives alongside `EQQ`.
- Status: Proven. The audited mixed-pattern split is therefore separate from the threshold theorem: rank-2 mixed-aware special case, rank-4 mixed exception, and ranks `3`, `5`, `6` supporting the `EYY`-first audited pattern.
- Status: Proven. Therefore a stronger physically justified minimal-sector theorem is already on the no-go branch unless explicit suppression, ordering, or background-restriction assumptions are added family by family.

## Scope Separation

- Status: Proven. The exact-current-set electric-only theorem and the corrected exact-current-set `E/B` theorem are special audited instances of the broader positive finite-family target.
- Status: Proven. The stronger physically justified minimal-sector theorem is a separate claim and is not rescued by current audits.
- Status: Proven. The negative branch is the family-admission no-go for minimal-sector uniqueness.
- Status: Proven. The positive branch is the finite-family fixed-order collapse theorem inside the current theorem domain.

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
19. Status: Proven. [`family-envelope-theorem.md`](family-envelope-theorem.md), [`family-envelope-table.md`](family-envelope-table.md), [`../lemmas/24-family-envelope-census.md`](../lemmas/24-family-envelope-census.md), [`../lemmas/25-next-unaudited-family-gate.md`](../lemmas/25-next-unaudited-family-gate.md), and [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py) first separated audited-set sufficiency from envelope completeness in the earlier family-envelope census stage before the later irreducible-envelope closure superseded the rank-by-rank gate march.
20. Status: Proven. [`rank3-family-ordering.md`](rank3-family-ordering.md), [`../lemmas/30-rank3-family-admission.md`](../lemmas/30-rank3-family-admission.md), [`../lemmas/31-rank3-threshold-formula.md`](../lemmas/31-rank3-threshold-formula.md), [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py), and [`../symbolic/r3_survivor_rank_check.py`](../symbolic/r3_survivor_rank_check.py) now show that `Rodd+` is a genuine new obstruction class rather than a derivative-generated artifact.
21. Status: Proven. [`../lemmas/32-r3-pairwise-composition-audit.md`](../lemmas/32-r3-pairwise-composition-audit.md), [`../lemmas/33-r3-augmented-composition-audit.md`](../lemmas/33-r3-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, and five-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}`.
22. Status: Proven. [`rank4-family-ordering.md`](rank4-family-ordering.md), [`../lemmas/34-rank4-family-admission.md`](../lemmas/34-rank4-family-admission.md), [`../lemmas/35-rank4-threshold-formula.md`](../lemmas/35-rank4-threshold-formula.md), [`../symbolic/r4_sector_delta4.py`](../symbolic/r4_sector_delta4.py), and [`../symbolic/r4_survivor_rank_check.py`](../symbolic/r4_survivor_rank_check.py) now show that `Reven4+` is a genuine new obstruction class rather than a trace-descended or derivative-generated artifact.
23. Status: Proven. [`../lemmas/36-r4-pairwise-composition-audit.md`](../lemmas/36-r4-pairwise-composition-audit.md), [`../lemmas/37-r4-augmented-composition-audit.md`](../lemmas/37-r4-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, quintuple, and six-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+}`.
24. Status: Proven. [`rank5-family-ordering.md`](rank5-family-ordering.md), [`../lemmas/38-rank5-family-admission.md`](../lemmas/38-rank5-family-admission.md), [`../lemmas/39-rank5-threshold-formula.md`](../lemmas/39-rank5-threshold-formula.md), [`../symbolic/r5_sector_delta4.py`](../symbolic/r5_sector_delta4.py), and [`../symbolic/r5_survivor_rank_check.py`](../symbolic/r5_survivor_rank_check.py) now show that `Rodd5+` is a genuine new obstruction class rather than a trace-descended or derivative-generated artifact.
25. Status: Proven. [`../lemmas/40-r5-pairwise-composition-audit.md`](../lemmas/40-r5-pairwise-composition-audit.md), [`../lemmas/41-r5-augmented-composition-audit.md`](../lemmas/41-r5-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, quintuple, six-family, and seven-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
26. Status: Proven. [`high-rank-audit-methodology.md`](high-rank-audit-methodology.md), [`../lemmas/42-rank3-exhaustiveness-check.md`](../lemmas/42-rank3-exhaustiveness-check.md), [`../lemmas/43-rank4-exhaustiveness-check.md`](../lemmas/43-rank4-exhaustiveness-check.md), [`../lemmas/44-rank5-exhaustiveness-check.md`](../lemmas/44-rank5-exhaustiveness-check.md), [`../symbolic/high_rank_family_enumerator.py`](../symbolic/high_rank_family_enumerator.py), and [`../symbolic/high_rank_diff_report.py`](../symbolic/high_rank_diff_report.py) now separate manual high-rank survivor bookkeeping from exhaustive contraction generation.
27. Status: Proven. The exhaustive high-rank audit validates the current manual rank-3 and rank-5 candidate-survivor lists, but it finds omitted rank-4 classes `\{EEQ,\ Q3,\ EEDtQ,\ EEEQ,\ EQQQ_1,\ EQQQ_2,\ GradEGradQ\}`.
28. Status: Proven. [`stf-tower-theorem.md`](stf-tower-theorem.md), [`../lemmas/42-stf-rankL-admission.md`](../lemmas/42-stf-rankL-admission.md), [`../lemmas/43-stf-rankL-threshold-formula.md`](../lemmas/43-stf-rankL-threshold-formula.md), [`stf-family-class-table.md`](stf-family-class-table.md), and [`../symbolic/stf_rankL_pattern_check.py`](../symbolic/stf_rankL_pattern_check.py) show that the attempted single STF-tower theorem for all genuine parity-even STF ranks `L \ge 3` fails at `L = 4`.
29. Status: Proven. [`rank6-family-ordering.md`](rank6-family-ordering.md), [`../lemmas/44-rank6-family-admission.md`](../lemmas/44-rank6-family-admission.md), [`../lemmas/45-rank6-threshold-formula.md`](../lemmas/45-rank6-threshold-formula.md), [`../symbolic/r6_sector_delta4.py`](../symbolic/r6_sector_delta4.py), and [`../symbolic/r6_survivor_rank_check.py`](../symbolic/r6_survivor_rank_check.py) audit `Reven6+` exhaustively and show that the first mixed layer is `EZZ`, not a new `EEY`-type even-rank exception.
30. Status: Proven. [`../lemmas/46-r6-pairwise-composition-audit.md`](../lemmas/46-r6-pairwise-composition-audit.md), [`../lemmas/47-r6-augmented-composition-audit.md`](../lemmas/47-r6-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, quintuple, six-family, seven-family, and eight-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}`.
31. Status: Proven. [`stf-self-witness-theorem.md`](stf-self-witness-theorem.md), [`../lemmas/48-stf-self-witness-theorem.md`](../lemmas/48-stf-self-witness-theorem.md), [`../lemmas/49-stf-threshold-theorem.md`](../lemmas/49-stf-threshold-theorem.md), and [`../symbolic/stf_self_witness_check.py`](../symbolic/stf_self_witness_check.py) promote the repeated audited higher-rank STF cases into a single universal self-witness threshold theorem for genuine parity-even STF primitive families of rank `L \ge 3`.
32. Status: Proven. [`stf-mixed-pattern-classification.md`](stf-mixed-pattern-classification.md), [`stf-family-class-table.md`](stf-family-class-table.md), and [`stf-tower-theorem.md`](stf-tower-theorem.md) now separate the closed threshold theorem from the still-split mixed-pattern story.
33. Status: Proven. [`irreducible-envelope-theorem.md`](irreducible-envelope-theorem.md), [`../lemmas/50-cartesian-irrep-reduction.md`](../lemmas/50-cartesian-irrep-reduction.md), and [`../lemmas/51-trace-descendant-absorption.md`](../lemmas/51-trace-descendant-absorption.md) replace the old rank-by-rank family march with a three-dimensional irreducible-envelope reduction.
34. Status: Proven. [`../lemmas/52-mixed-symmetry-family-gate.md`](../lemmas/52-mixed-symmetry-family-gate.md), [`mixed-symmetry-risk-register.md`](mixed-symmetry-risk-register.md), [`../symbolic/irrep_family_census.py`](../symbolic/irrep_family_census.py), and [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py) show that no explicit mixed-symmetry or otherwise non-STF primitive family survives inside the current theorem domain.
35. Status: Proven. [`finite-family-collapse-theorem.md`](finite-family-collapse-theorem.md), [`../lemmas/53-fixed-order-operator-finiteness.md`](../lemmas/53-fixed-order-operator-finiteness.md), and [`../symbolic/fixed_family_operator_count.py`](../symbolic/fixed_family_operator_count.py) now close the pre-reduction finite-candidate-operator step for any explicit finite admitted family catalog inside the current theorem domain.
36. Status: Proven. [`../lemmas/54-normal-form-basis-closure.md`](../lemmas/54-normal-form-basis-closure.md), [`reduction-rules.md`](reduction-rules.md), and [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) now close the finite-dimensional scalar normal-form quotient step.
37. Status: Proven. [`../lemmas/55-monopole-jet-collapse.md`](../lemmas/55-monopole-jet-collapse.md) now closes the finite analytic monopole jet step once `A5` is imposed.
38. Status: Proven. [`../lemmas/56-sensitivity-vs-wilson-split.md`](../lemmas/56-sensitivity-vs-wilson-split.md) now closes the finite sensitivity-versus-Wilson split once `A3`, `A5`, and `A8` are imposed.
39. Status: Proven. [`collapse-bridge-status.md`](collapse-bridge-status.md) makes the bridge status explicit and records that the former positive-branch bottleneck is now resolved positively.
40. Status: Proven. [`nonanalytic-counterexample-program.md`](nonanalytic-counterexample-program.md), [`../lemmas/57-local-finite-but-nonanalytic.md`](../lemmas/57-local-finite-but-nonanalytic.md), and [`../lemmas/58-nonanalytic-jet-failure.md`](../lemmas/58-nonanalytic-jet-failure.md) isolate the `A5` boundary and show that locality plus finite operator closure survive even when the analytic monopole jet fails.
41. Status: Proven. [`../lemmas/59-nonanalytic-sensitivity-replacement.md`](../lemmas/59-nonanalytic-sensitivity-replacement.md) separates non-Taylor monopole germ data from higher-multipole Wilson coefficients.
42. Status: Proven. [`nonanalytic-model-table.md`](nonanalytic-model-table.md) and [`../symbolic/nonanalytic_jet_demo.py`](../symbolic/nonanalytic_jet_demo.py) record the explicit analytic control, smooth-flat, and threshold models side by side.

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
- Status: Proven. The rank-4 clarification still shows that trace descendants and derivative-generated rank-4 blocks do not absorb a genuine primitive STF family `Q_{ijkl}`.
- Status: Proven. The smallest explicit rank-4 obstruction remains `Q2`, while the exhaustive audit now shows that the mixed cubic layer at weight `3` includes both the previously recorded `EQQ` and the omitted `EEQ`.
- Status: Proven. The rank-5 clarification shows that trace descendants and derivative-generated rank-5 blocks do not absorb a genuine primitive STF family `U_{ijklm}`.
- Status: Proven. The smallest explicit rank-5 obstruction is `U2`, with first mixed witness `EUU`.
- Status: Proven. The stronger mixed-pattern STF-tower theorem for all genuine parity-even ranks `L \ge 3` is false at the first audited even-rank member `L = 4`.
- Status: Proven. The weaker self-witness threshold theorem for all genuine parity-even STF primitive families of rank `L \ge 3` is now proven.
- Status: Proven. Necessary suppression budgets are explicit for the audited classes: for rank-2 STF families the self-only lower bound is `w_X \ge 3` but the mixed-aware necessary threshold is `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness; for bare scalar families the self-only sharp threshold is `w_S \ge 5`; for derivative-only scalar families the tied-sharp threshold is `w_D \ge 3`; for genuine vector families the self-only sharp threshold is `w_V \ge 3`; for genuine rank-3 STF families the self-only sharp threshold is `w_T \ge 3`; for genuine rank-4 STF families the self-only sharp threshold is `w_Q \ge 3`; for genuine rank-5 STF families the self-only sharp threshold is `w_U \ge 3`; for genuine rank-6 STF families the self-only sharp threshold is `w_Z \ge 3`.
- Status: Proven. The pairwise, triple, and quadruple composition audits already found no surviving operator beyond the baseline electric sector once the current audited thresholds were imposed across the enlarged audited family set `{R2, R0a, R0b, R1}`.
- Status: Proven. The post-Rodd+ pairwise, triple, quadruple, and five-family composition audits likewise find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}`.
- Status: Proven. Therefore the current thresholds are jointly sufficient for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}` at `\Delta \le 4`, and no explicit post-Rodd+ cross-family obstruction has appeared in the audited set.
- Status: Proven. The raw R3-extended survivor list is not linearly independent; its first exact quartic mixed dependence relation is `-E2T2/2 + 2 E2T2_mixed_1 - E2T2_mixed_2 + E2T2_mixed_3 = 0`.
- Status: Proven. The raw R3-extended `21`-element survivor list is bookkeeping only and should not be read as a corrected basis statement; the composition audit uses the thresholded rank-3 family instead.
- Status: Proven. The old raw R4-extended survivor bookkeeping is not linearly independent; its first exact quartic mixed dependence relation is `-E2Q2/2 + 2 E2Q2_mixed_1 - E2Q2_mixed_2 + E2Q2_mixed_3 = 0`.
- Status: Proven. The old raw R4 bookkeeping with `22` labels, rank `19`, and nullity `3` now applies only to the pre-exhaustiveness manual subset, not to an exhaustive corrected rank-4 basis.
- Status: Conjectural. The raw R5-extended survivor list is not linearly independent; its first sample-extracted quartic mixed dependence relation is `E2U2 - 4 E2U2_mixed_1 + 2 E2U2_mixed_2 - 2 E2U2_mixed_3 = 0`.
- Status: Conjectural. The raw R5-extended `23`-element survivor list has sample-stable rank `19` and nullity `4`, and it is bookkeeping only rather than a corrected basis statement.
- Status: Conjectural. The raw R6-extended survivor list is not linearly independent; its first sample-revalidated mixed quartic dependence relation is `E2Z2 + 2 E2Z2_mixed_1 - 2 E2Z2_mixed_2 - 4 E2Z2_mixed_3 = 0`.
- Status: Conjectural. The raw R6-extended `22`-label family list has sample-stable rank `16` and nullity `6`, and it is bookkeeping only rather than a corrected basis statement.
- Status: Proven. The irreducible family-envelope theorem now closes the parity-even nonspinning local MVP family domain on the audited scalar/vector/STF classes.
- Status: Proven. The post-Reven4+ pairwise, triple, quadruple, quintuple, and six-family composition audits find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+}`.
- Status: Proven. The post-Rodd5+ pairwise, triple, quadruple, quintuple, six-family, and seven-family composition audits likewise find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
- Status: Proven. Therefore the current thresholds are jointly sufficient for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}` at `\Delta \le 4`, and no explicit post-Rodd5+ cross-family obstruction has appeared in the audited set.
- Status: Proven. The `Reven6+` audit now resolves the next even-rank gate as a genuine obstruction class with smallest witness `Z2` and first mixed witness `EZZ`.
- Status: Proven. The post-Reven6+ pairwise, triple, quadruple, quintuple, six-family, seven-family, and eight-family composition audits find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}`.
- Status: Proven. Therefore the current thresholds are jointly sufficient for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}` at `\Delta \le 4`, and no explicit post-Reven6+ cross-family obstruction has appeared in the audited set.
- Status: Proven. The former live bottleneck `"prove or refute irreducible family-envelope closure beyond the audited scalar/vector/STF classes"` is now resolved positively.
- Status: Proven. The former live bottleneck `"prove or refute the positive finite-family collapse theorem once the irreducible scalar/vector/STF envelope is fixed"` is now resolved positively.
- Status: Proven. Fixed-order operator-space non-closure is no longer a live failure mode inside the current theorem domain once finite admitted family content and irreducible-envelope closure are imposed.
- Status: Proven. Finite normal-form basis closure is no longer a live failure mode inside the current theorem domain once the explicit quotient rules are imposed.
- Status: Proven. The positive finite-family collapse theorem now closes inside the current theorem domain.
- Status: Proven. The former boundary-stress bottleneck `"construct or refute the smallest local nonanalytic counterexample to analytic monopole jet collapse (A5 boundary)"` is now resolved positively by the one-coordinate smooth flat monopole model `m_A(Y)=m_0+\alpha e^{-1/Y^2}\Theta(Y)`.
- Status: Proven. That model preserves locality, finite admitted family content, fixed-order operator finiteness, and finite scalar normal-form closure, but it destroys the finite analytic Taylor jet at the reference point.
- Status: Proven. The replacement bookkeeping is non-Taylor monopole germ data and remains separate from the residual Wilson coefficients.
- Status: Counterexample candidate. The remaining explicit risks are assumption-drop failures: nonanalytic monopole dependence at `A5`, locality failure of the sensitivity/Wilson split at `A3`, failure of finite admitted family content at `A8`, or escape beyond the present theorem domain.

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
- Status: Proven. [`../lemmas/32-r3-pairwise-composition-audit.md`](../lemmas/32-r3-pairwise-composition-audit.md)
- Status: Proven. [`../lemmas/33-r3-augmented-composition-audit.md`](../lemmas/33-r3-augmented-composition-audit.md)
- Status: Proven. [`../lemmas/34-rank4-family-admission.md`](../lemmas/34-rank4-family-admission.md)
- Status: Proven. [`../lemmas/35-rank4-threshold-formula.md`](../lemmas/35-rank4-threshold-formula.md)
- Status: Proven. [`../lemmas/36-r4-pairwise-composition-audit.md`](../lemmas/36-r4-pairwise-composition-audit.md)
- Status: Proven. [`../lemmas/37-r4-augmented-composition-audit.md`](../lemmas/37-r4-augmented-composition-audit.md)
- Status: Proven. [`../lemmas/38-rank5-family-admission.md`](../lemmas/38-rank5-family-admission.md)
- Status: Proven. [`../lemmas/39-rank5-threshold-formula.md`](../lemmas/39-rank5-threshold-formula.md)
- Status: Proven. [`../lemmas/40-r5-pairwise-composition-audit.md`](../lemmas/40-r5-pairwise-composition-audit.md)
- Status: Proven. [`../lemmas/41-r5-augmented-composition-audit.md`](../lemmas/41-r5-augmented-composition-audit.md)
- Status: Proven. [`../lemmas/48-stf-self-witness-theorem.md`](../lemmas/48-stf-self-witness-theorem.md)
- Status: Proven. [`../lemmas/49-stf-threshold-theorem.md`](../lemmas/49-stf-threshold-theorem.md)
- Status: Proven. [`../lemmas/50-cartesian-irrep-reduction.md`](../lemmas/50-cartesian-irrep-reduction.md)
- Status: Proven. [`../lemmas/51-trace-descendant-absorption.md`](../lemmas/51-trace-descendant-absorption.md)
- Status: Proven. [`../lemmas/52-mixed-symmetry-family-gate.md`](../lemmas/52-mixed-symmetry-family-gate.md)
- Status: Proven. [`../lemmas/53-fixed-order-operator-finiteness.md`](../lemmas/53-fixed-order-operator-finiteness.md)
- Status: Proven. [`../lemmas/54-normal-form-basis-closure.md`](../lemmas/54-normal-form-basis-closure.md)
- Status: Proven. [`../lemmas/55-monopole-jet-collapse.md`](../lemmas/55-monopole-jet-collapse.md)
- Status: Proven. [`../lemmas/56-sensitivity-vs-wilson-split.md`](../lemmas/56-sensitivity-vs-wilson-split.md)
- Status: Proven. [`../lemmas/57-local-finite-but-nonanalytic.md`](../lemmas/57-local-finite-but-nonanalytic.md)
- Status: Proven. [`../lemmas/58-nonanalytic-jet-failure.md`](../lemmas/58-nonanalytic-jet-failure.md)
- Status: Proven. [`../lemmas/59-nonanalytic-sensitivity-replacement.md`](../lemmas/59-nonanalytic-sensitivity-replacement.md)
- Status: Conjectural. [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md)
- Status: Conjectural. [`eb-conditional-collapse.md`](eb-conditional-collapse.md)
- Status: Proven. [`finite-family-collapse-theorem.md`](finite-family-collapse-theorem.md)
- Status: Proven. [`collapse-bridge-status.md`](collapse-bridge-status.md)
- Status: Proven. [`nonanalytic-counterexample-program.md`](nonanalytic-counterexample-program.md)
- Status: Proven. [`nonanalytic-model-table.md`](nonanalytic-model-table.md)
- Status: Conjectural. [`power-counting.md`](power-counting.md)
- Status: Proven. [`family-admission-no-go.md`](family-admission-no-go.md)
- Status: Proven. [`family-admission-theorem.md`](family-admission-theorem.md)
- Status: Proven. [`family-envelope-theorem.md`](family-envelope-theorem.md)
- Status: Proven. [`family-envelope-table.md`](family-envelope-table.md)
- Status: Proven. [`irreducible-envelope-theorem.md`](irreducible-envelope-theorem.md)
- Status: Proven. [`mixed-symmetry-risk-register.md`](mixed-symmetry-risk-register.md)
- Status: Proven. [`rank4-family-ordering.md`](rank4-family-ordering.md)
- Status: Proven. [`rank5-family-ordering.md`](rank5-family-ordering.md)
- Status: Proven. [`rank6-family-ordering.md`](rank6-family-ordering.md)
- Status: Proven. [`stf-self-witness-theorem.md`](stf-self-witness-theorem.md)
- Status: Proven. [`stf-mixed-pattern-classification.md`](stf-mixed-pattern-classification.md)
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
- Status: Proven. [`../symbolic/irrep_family_census.py`](../symbolic/irrep_family_census.py)
- Status: Proven. [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py)
- Status: Proven. [`../symbolic/fixed_family_operator_count.py`](../symbolic/fixed_family_operator_count.py)
- Status: Proven. [`../symbolic/nonanalytic_jet_demo.py`](../symbolic/nonanalytic_jet_demo.py)
- Status: Proven. [`../symbolic/r1_sector_delta4.py`](../symbolic/r1_sector_delta4.py)
- Status: Proven. [`../symbolic/r1_survivor_rank_check.py`](../symbolic/r1_survivor_rank_check.py)
- Status: Proven. [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py)
- Status: Proven. [`../symbolic/r3_survivor_rank_check.py`](../symbolic/r3_survivor_rank_check.py)
- Status: Proven. [`../symbolic/r4_sector_delta4.py`](../symbolic/r4_sector_delta4.py)
- Status: Proven. [`../symbolic/r4_survivor_rank_check.py`](../symbolic/r4_survivor_rank_check.py)
- Status: Proven. [`../symbolic/r5_sector_delta4.py`](../symbolic/r5_sector_delta4.py)
- Status: Proven. [`../symbolic/r5_survivor_rank_check.py`](../symbolic/r5_survivor_rank_check.py)
- Status: Proven. [`../symbolic/r6_sector_delta4.py`](../symbolic/r6_sector_delta4.py)
- Status: Proven. [`../symbolic/r6_survivor_rank_check.py`](../symbolic/r6_survivor_rank_check.py)
- Status: Proven. [`../symbolic/stf_self_witness_check.py`](../symbolic/stf_self_witness_check.py)

## Failure Triggers

- Status: Counterexample candidate. An explicit `chi_A` state invalidates the no-state reduction step.
- Status: Counterexample candidate. Nonanalytic activation invalidates the analytic collapse step.
- Status: Counterexample candidate. Hereditary couplings invalidate locality before basis closure is even applied.
- Status: Counterexample candidate. Infinite or uncontrolled external field content invalidates the fixed-order closure theorem candidate.
- Status: Counterexample candidate. Any additional `\Delta \le 4` scalar contraction from the exact current electric-only primitive set that evades the explicit enumeration/reduction pipeline obstructs the electric-only exact-current-set theorem.
- Status: Counterexample candidate. Any further physically admissible primitive family that yields new survivors beyond the corrected explicit audited family catalogs reinforces the uniqueness no-go branch.
- Status: Counterexample candidate. Any further physically admissible primitive family obstructs the positive finite-family collapse branch only if it also destroys fixed-order finiteness or the family-conditioned reduction path.
