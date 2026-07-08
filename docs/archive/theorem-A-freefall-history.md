# Theorem A Candidate: Free-Fall Fixed-Order Collapse

## Main Positive Target: Finite-Family Fixed-Order Collapse Theorem

- Status: Note. Work in the free-fall MVP sector only: nonspinning, nearly spherical, parity-even, local worldline EFT, fixed truncation order `\Delta \le 4`, and finite admitted external field content.
- Status: Note. Within that theorem domain, the irreducible primitive-family envelope is already closed on the audited scalar/vector/STF classes.
- Status: Note. For any admitted primitive-family spectrum drawn from that closed envelope which is locally finite below `\Delta \le 4`, the candidate parity-even local free-fall scalar operator space at `\Delta \le 4` is finite before reduction.
- Status: Note. Modulo the explicit family-envelope, total-derivative, lower-order-EOM, algebraic, and extracted linear-dependence relations already fixed in this repo, the reduced scalar operator space is finite-dimensional and admits a finite normal-form basis `Y^I`.
- Status: Note. Under the locality and analyticity assumptions `A3` and `A5`, the monopole term `m_A(Y)` therefore has a finite fixed-order Taylor jet in finitely many scalar normal-form coordinates.
- Status: Note. The residual higher-multipole sector is then carried by finitely many Wilson coefficients `C_{A,\alpha}` that are separate from the monopole sensitivity coefficients.
- Status: Note. Therefore the positive finite-family collapse theorem closes inside the current theorem domain.

## Positive Bridge Layers

- Status: Note. Layer 1 is irreducible family-envelope closure; this is already closed positively on the audited scalar/vector/STF classes.
- Status: Note. Layer 2 is fixed-order candidate operator finiteness; this now follows from local weight-spectrum finiteness below the cutoff, positive weights, and fixed `\Delta \le 4`.
- Status: Note. Layer 3 is finite normal-form basis closure; this now follows because the reduced operator space is a quotient of a finite-dimensional candidate space.
- Status: Note. Layer 4 is analytic monopole jet collapse; this now follows conditionally on `A5` once the finite scalar basis `Y^I` is fixed.
- Status: Note. Layer 5 is the sensitivity/Wilson split; this now follows conditionally on `A3`, `A5`, and the sharpened `A8`.

## Exact-Current-Set Specializations

### Electric-Only Exact-Current-Set Theorem

- Status: Conjectural. For the exact current electric-only primitive set at `\Delta \le 4`, the admissible local free-fall scalar operator space closes on

```math
\{E2,\ E3,\ E2^2,\ dotE2,\ gradE2,\ divE2,\ mixedGradE2\}.
```

- Status: Note. This corrected seven-element electric-only basis is linearly independent.
- Status: Conjectural. This is an exact-current-set theorem candidate, not yet a physically justified minimal-sector theorem.

### Corrected `E/B` Exact-Current-Set Theorem

- Status: Note. The raw `E/B` survivor candidate set is not linearly independent; one explicit mixed quartic dependence relation reduces it.
- Status: Conjectural. After that correction, the exact-current-set `E/B` basis at `\Delta \le 4` is

```math
\{E2,\ B2,\ E3,\ EB2,\ E2^2,\ B2^2,\ dotE2,\ dotB2,\ EBDtB,\ E2B2,\ EB\_sq,\ TrE2B2,\ gradE2,\ divE2,\ mixedGradE2,\ gradB2,\ divB2,\ mixedGradB2\}.
```

- Status: Note. This corrected `E/B` basis is linearly independent.
- Status: Conjectural. This is likewise an exact-current-set theorem candidate, not yet a physically justified minimal-sector theorem.

## Negative Branch: Family-Admission No-Go For Minimal-Sector Uniqueness

- Status: Note. Minimal-sector uniqueness is not the positive theorem target anymore.
- Status: Note. Across the explicitly audited family classes, unsuppressed family admission always yields a new low-order survivor.
- Status: Note. For STF rank-2 family admission, the self witness `X2` and the mixed witness `EX` both appear at weight `2`; the audited explicit self instance is `B2`.
- Status: Note. For unsuppressed rank-0 scalar-family admission, the self witness `S` appears before the first mixed witness `SE2`.
- Status: Note. For derivative-only or shift-symmetric scalar-family admission, the self witness `dotS2` and the mixed witness `DtS_E2` both first appear at weight `4`.
- Status: Note. For a genuine parity-even rank-1 vector-family admission, the self witness `V2` appears at weight `2` and the first mixed witness `EVV` appears at weight `3`.
- Status: Note. Rank `2` remains a separate mixed-aware STF special case because the mixed quadratic witness `EX` exists.
- Status: Note. For every genuinely new parity-even STF primitive family with rank `L \ge 3`, there is no linear scalar witness, no mixed quadratic `E Y_L` witness, and the quadratic norm `Y2` always exists.
- Status: Note. Therefore the higher-rank STF branch obeys a universal self-only sharp threshold `w_Y \ge 3` at `\Delta \le 4`.
- Status: Note. The stronger mixed-pattern theorem still fails at audited rank `L = 4`, because the cubic mixed class `EEQ` survives alongside `EQQ`.
- Status: Note. The audited mixed-pattern split is therefore separate from the threshold theorem: rank-2 mixed-aware special case, rank-4 mixed exception, and ranks `3`, `5`, `6` supporting the `EYY`-first audited pattern.
- Status: Note. Therefore a stronger physically justified minimal-sector theorem is already on the no-go branch unless explicit suppression, ordering, or background-restriction assumptions are added family by family.

## Scope Separation

- Status: Note. The exact-current-set electric-only theorem and the corrected exact-current-set `E/B` theorem are special audited instances of the broader positive finite-family target.
- Status: Note. The stronger physically justified minimal-sector theorem is a separate claim and is not rescued by current audits.
- Status: Note. The negative branch is the family-admission no-go for minimal-sector uniqueness.
- Status: Note. The positive branch is the finite-family fixed-order collapse theorem inside the current theorem domain.

## Proof Route

1. Status: Note. Use the fixed-order counting rule in [`power-counting.md`](power-counting.md) to define a bounded-weight operator search space.
2. Status: Note. Local weight-spectrum finiteness below the cutoff and positive weights imply a finite set of decorated primitive building blocks at order `\Delta \le 4`.
3. Status: Note. [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) exhaustively enumerates all parity-even scalar contractions from the exact current electric-only primitive set.
4. Status: Note. [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) reduces all total-derivative, lower-order-EOM, and single-family Cayley-Hamilton-reducible electric classes explicitly.
5. Status: Note. [`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py) shows that the corrected seven electric-only survivors are linearly independent as operators.
6. Status: Note. [`../symbolic/eb_sector_delta4.py`](../symbolic/eb_sector_delta4.py) shows that admitting the magnetic family yields a finite enlarged `E/B` candidate list rather than a fixed-order blowup.
7. Status: Note. [`../symbolic/eb_survivor_rank_check.py`](../symbolic/eb_survivor_rank_check.py) extracts the first explicit `E/B` dependence relation and the corrected linearly independent `18`-element basis.
8. Status: Note. [`../symbolic/es_sector_delta4.py`](../symbolic/es_sector_delta4.py) shows that adjoining one scalar-like external family still yields a finite enlarged `E/B+scalar` candidate list rather than a fixed-order blowup.
9. Status: Note. [`../symbolic/es_survivor_rank_check.py`](../symbolic/es_survivor_rank_check.py) shows that the corrected `E/B+scalar` `33`-element basis is linearly independent.
10. Status: Note. [`../symbolic/shift_scalar_sector_delta4.py`](../symbolic/shift_scalar_sector_delta4.py) shows that even derivative-only scalar admission still yields new weight-`4` survivors beyond the corrected `E/B` sector.
11. Status: Note. [`../symbolic/family_witness_map.py`](../symbolic/family_witness_map.py) records the current family-class witness table and the exact theorem layer each audited class obstructs.
12. Status: Note. [`family-admission-theorem.md`](family-admission-theorem.md) promotes the repeated family-attack pattern into an explicit no-go statement for the audited family classes.
13. Status: Note. [`suppression-budget-theorem.md`](suppression-budget-theorem.md), [`../lemmas/17-witness-threshold-classification.md`](../lemmas/17-witness-threshold-classification.md), and [`../symbolic/witness_threshold_map.py`](../symbolic/witness_threshold_map.py) promote the audited witnesses into necessary suppression thresholds for minimal-sector uniqueness at `\Delta_{\max} = 4`.
14. Status: Note. [`mixed-witness-thresholds.md`](mixed-witness-thresholds.md), [`../lemmas/18-rank2-mixed-witness-threshold.md`](../lemmas/18-rank2-mixed-witness-threshold.md), [`../lemmas/19-rank0-mixed-witness-threshold.md`](../lemmas/19-rank0-mixed-witness-threshold.md), and [`../symbolic/mixed_witness_map.py`](../symbolic/mixed_witness_map.py) classify the first self and first mixed witnesses for the audited family classes.
15. Status: Note. [`threshold-formulas.md`](threshold-formulas.md), [`../lemmas/20-rank2-threshold-formula.md`](../lemmas/20-rank2-threshold-formula.md), [`../lemmas/21-rank0-threshold-formula.md`](../lemmas/21-rank0-threshold-formula.md), and [`../symbolic/threshold_formula_check.py`](../symbolic/threshold_formula_check.py) lift the audited family classes into formula-level threshold statements.
16. Status: Note. [`budget-composition-theorem.md`](budget-composition-theorem.md), [`../lemmas/22-pairwise-composition-audit.md`](../lemmas/22-pairwise-composition-audit.md), [`../lemmas/23-three-family-composition-audit.md`](../lemmas/23-three-family-composition-audit.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) first closed the old audited-set composition layer for `R2`, `R0a`, and `R0b`.
17. Status: Note. [`vector-family-ordering.md`](vector-family-ordering.md), [`../lemmas/26-rank1-family-admission.md`](../lemmas/26-rank1-family-admission.md), [`../lemmas/27-rank1-threshold-formula.md`](../lemmas/27-rank1-threshold-formula.md), [`../symbolic/r1_sector_delta4.py`](../symbolic/r1_sector_delta4.py), and [`../symbolic/r1_survivor_rank_check.py`](../symbolic/r1_survivor_rank_check.py) show that a genuine local parity-even rank-1 vector family is a real obstruction class rather than a derivative-generated artifact.
18. Status: Note. [`../lemmas/28-r1-pairwise-composition-audit.md`](../lemmas/28-r1-pairwise-composition-audit.md), [`../lemmas/29-r1-augmented-composition-audit.md`](../lemmas/29-r1-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, and quadruple composition layer for the enlarged audited family set `{R2, R0a, R0b, R1}`.
19. Status: Note. [`family-envelope-theorem.md`](family-envelope-theorem.md), [`family-envelope-table.md`](family-envelope-table.md), [`../lemmas/24-family-envelope-census.md`](../lemmas/24-family-envelope-census.md), [`../lemmas/25-next-unaudited-family-gate.md`](../lemmas/25-next-unaudited-family-gate.md), and [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py) first separated audited-set sufficiency from envelope completeness in the earlier family-envelope census stage before the later irreducible-envelope closure superseded the rank-by-rank gate march.
20. Status: Note. [`rank3-family-ordering.md`](rank3-family-ordering.md), [`../lemmas/30-rank3-family-admission.md`](../lemmas/30-rank3-family-admission.md), [`../lemmas/31-rank3-threshold-formula.md`](../lemmas/31-rank3-threshold-formula.md), [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py), and [`../symbolic/r3_survivor_rank_check.py`](../symbolic/r3_survivor_rank_check.py) now show that `Rodd+` is a genuine new obstruction class rather than a derivative-generated artifact.
21. Status: Note. [`../lemmas/32-r3-pairwise-composition-audit.md`](../lemmas/32-r3-pairwise-composition-audit.md), [`../lemmas/33-r3-augmented-composition-audit.md`](../lemmas/33-r3-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, and five-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}`.
22. Status: Note. [`rank4-family-ordering.md`](rank4-family-ordering.md), [`../lemmas/34-rank4-family-admission.md`](../lemmas/34-rank4-family-admission.md), [`../lemmas/35-rank4-threshold-formula.md`](../lemmas/35-rank4-threshold-formula.md), [`../symbolic/r4_sector_delta4.py`](../symbolic/r4_sector_delta4.py), and [`../symbolic/r4_survivor_rank_check.py`](../symbolic/r4_survivor_rank_check.py) now show that `Reven4+` is a genuine new obstruction class rather than a trace-descended or derivative-generated artifact.
23. Status: Note. [`../lemmas/36-r4-pairwise-composition-audit.md`](../lemmas/36-r4-pairwise-composition-audit.md), [`../lemmas/37-r4-augmented-composition-audit.md`](../lemmas/37-r4-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, quintuple, and six-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+}`.
24. Status: Note. [`rank5-family-ordering.md`](rank5-family-ordering.md), [`../lemmas/38-rank5-family-admission.md`](../lemmas/38-rank5-family-admission.md), [`../lemmas/39-rank5-threshold-formula.md`](../lemmas/39-rank5-threshold-formula.md), [`../symbolic/r5_sector_delta4.py`](../symbolic/r5_sector_delta4.py), and [`../symbolic/r5_survivor_rank_check.py`](../symbolic/r5_survivor_rank_check.py) now show that `Rodd5+` is a genuine new obstruction class rather than a trace-descended or derivative-generated artifact.
25. Status: Note. [`../lemmas/40-r5-pairwise-composition-audit.md`](../lemmas/40-r5-pairwise-composition-audit.md), [`../lemmas/41-r5-augmented-composition-audit.md`](../lemmas/41-r5-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, quintuple, six-family, and seven-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
26. Status: Note. [`high-rank-audit-methodology.md`](high-rank-audit-methodology.md), [`../lemmas/42-rank3-exhaustiveness-check.md`](../lemmas/42-rank3-exhaustiveness-check.md), [`../lemmas/43-rank4-exhaustiveness-check.md`](../lemmas/43-rank4-exhaustiveness-check.md), [`../lemmas/44-rank5-exhaustiveness-check.md`](../lemmas/44-rank5-exhaustiveness-check.md), [`../symbolic/high_rank_family_enumerator.py`](../symbolic/high_rank_family_enumerator.py), and [`../symbolic/high_rank_diff_report.py`](../symbolic/high_rank_diff_report.py) now separate manual high-rank survivor bookkeeping from exhaustive contraction generation.
27. Status: Note. The exhaustive high-rank audit validates the current manual rank-3 and rank-5 candidate-survivor lists, but it finds omitted rank-4 classes `\{EEQ,\ Q3,\ EEDtQ,\ EEEQ,\ EQQQ_1,\ EQQQ_2,\ GradEGradQ\}`.
28. Status: Note. [`stf-tower-theorem.md`](stf-tower-theorem.md), [`../lemmas/42-stf-rankL-admission.md`](../lemmas/42-stf-rankL-admission.md), [`../lemmas/43-stf-rankL-threshold-formula.md`](../lemmas/43-stf-rankL-threshold-formula.md), [`stf-family-class-table.md`](stf-family-class-table.md), and [`../symbolic/stf_rankL_pattern_check.py`](../symbolic/stf_rankL_pattern_check.py) show that the attempted single STF-tower theorem for all genuine parity-even STF ranks `L \ge 3` fails at `L = 4`.
29. Status: Note. [`rank6-family-ordering.md`](rank6-family-ordering.md), [`../lemmas/44-rank6-family-admission.md`](../lemmas/44-rank6-family-admission.md), [`../lemmas/45-rank6-threshold-formula.md`](../lemmas/45-rank6-threshold-formula.md), [`../symbolic/r6_sector_delta4.py`](../symbolic/r6_sector_delta4.py), and [`../symbolic/r6_survivor_rank_check.py`](../symbolic/r6_survivor_rank_check.py) audit `Reven6+` exhaustively and show that the first mixed layer is `EZZ`, not a new `EEY`-type even-rank exception.
30. Status: Note. [`../lemmas/46-r6-pairwise-composition-audit.md`](../lemmas/46-r6-pairwise-composition-audit.md), [`../lemmas/47-r6-augmented-composition-audit.md`](../lemmas/47-r6-augmented-composition-audit.md), [`composition-status.md`](composition-status.md), and [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) re-close the pairwise, triple, quadruple, quintuple, six-family, seven-family, and eight-family composition layer for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}`.
31. Status: Note. [`stf-self-witness-theorem.md`](stf-self-witness-theorem.md), [`../lemmas/48-stf-self-witness-theorem.md`](../lemmas/48-stf-self-witness-theorem.md), [`../lemmas/49-stf-threshold-theorem.md`](../lemmas/49-stf-threshold-theorem.md), and [`../symbolic/stf_self_witness_check.py`](../symbolic/stf_self_witness_check.py) promote the repeated audited higher-rank STF cases into a single universal self-witness threshold theorem for genuine parity-even STF primitive families of rank `L \ge 3`.
32. Status: Note. [`stf-mixed-pattern-classification.md`](stf-mixed-pattern-classification.md), [`stf-family-class-table.md`](stf-family-class-table.md), and [`stf-tower-theorem.md`](stf-tower-theorem.md) now separate the closed threshold theorem from the still-split mixed-pattern story.
33. Status: Note. [`irreducible-envelope-theorem.md`](irreducible-envelope-theorem.md), [`../lemmas/50-cartesian-irrep-reduction.md`](../lemmas/50-cartesian-irrep-reduction.md), and [`../lemmas/51-trace-descendant-absorption.md`](../lemmas/51-trace-descendant-absorption.md) replace the old rank-by-rank family march with a three-dimensional irreducible-envelope reduction.
34. Status: Note. [`../lemmas/52-mixed-symmetry-family-gate.md`](../lemmas/52-mixed-symmetry-family-gate.md), [`mixed-symmetry-risk-register.md`](mixed-symmetry-risk-register.md), [`../symbolic/irrep_family_census.py`](../symbolic/irrep_family_census.py), and [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py) show that no explicit mixed-symmetry or otherwise non-STF primitive family survives inside the current theorem domain.
35. Status: Note. [`finite-family-collapse-theorem.md`](finite-family-collapse-theorem.md), [`../lemmas/53-fixed-order-operator-finiteness.md`](../lemmas/53-fixed-order-operator-finiteness.md), and [`../symbolic/fixed_family_operator_count.py`](../symbolic/fixed_family_operator_count.py) now close the pre-reduction finite-candidate-operator step for any explicit finite admitted family catalog inside the current theorem domain.
36. Status: Note. [`../lemmas/54-normal-form-basis-closure.md`](../lemmas/54-normal-form-basis-closure.md), [`reduction-rules.md`](reduction-rules.md), and [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) now close the finite-dimensional scalar normal-form quotient step.
37. Status: Note. [`../lemmas/55-monopole-jet-collapse.md`](../lemmas/55-monopole-jet-collapse.md) now closes the finite analytic monopole jet step once `A5` is imposed.
38. Status: Note. [`../lemmas/56-sensitivity-vs-wilson-split.md`](../lemmas/56-sensitivity-vs-wilson-split.md) now closes the finite sensitivity-versus-Wilson split once `A3`, `A5`, and the sharpened `A8` are imposed.
39. Status: Note. [`weight-spectrum-theorem.md`](weight-spectrum-theorem.md), [`../lemmas/69-local-weight-spectrum-finiteness.md`](../lemmas/69-local-weight-spectrum-finiteness.md), [`../lemmas/70-infinite-low-weight-tower-counterexample.md`](../lemmas/70-infinite-low-weight-tower-counterexample.md), [`../lemmas/71-a8-sharpness.md`](../lemmas/71-a8-sharpness.md), and [`../symbolic/weight_spectrum_demo.py`](../symbolic/weight_spectrum_demo.py) sharpen the old finite-family-content assumption to local weight-spectrum finiteness below the cutoff and record the explicit infinite low-weight tower failure when it is dropped.
40. Status: Note. [`collapse-bridge-status.md`](collapse-bridge-status.md) makes the bridge status explicit and records that the former positive-branch bottleneck is now resolved positively.
41. Status: Note. [`nonanalytic-counterexample-program.md`](nonanalytic-counterexample-program.md), [`../lemmas/57-local-finite-but-nonanalytic.md`](../lemmas/57-local-finite-but-nonanalytic.md), and [`../lemmas/58-nonanalytic-jet-failure.md`](../lemmas/58-nonanalytic-jet-failure.md) isolate the `A5` boundary and show that locality plus finite operator closure survive even when the analytic monopole jet fails.
42. Status: Note. [`../lemmas/59-nonanalytic-sensitivity-replacement.md`](../lemmas/59-nonanalytic-sensitivity-replacement.md) separates non-Taylor monopole germ data from higher-multipole Wilson coefficients.
43. Status: Note. [`nonanalytic-model-table.md`](nonanalytic-model-table.md) and [`../symbolic/nonanalytic_jet_demo.py`](../symbolic/nonanalytic_jet_demo.py) record the explicit analytic control, smooth-flat, and threshold models side by side.

## Current Verdict

- Status: Note. The old five-scalar electric target was false because it omitted `divE2` and `mixedGradE2`.
- Status: Note. The corrected seven-scalar electric-only basis is linearly independent for the exact current primitive set.
- Status: Note. The corrected `E/B` basis is linearly independent for the corrected exact-current-set `E/B` primitive sector.
- Status: Note. The stronger physically justified minimal-sector theorem should not presently be treated as alive.
- Status: Note. Minimal-sector uniqueness is already effectively dead across the audited unsuppressed family classes.
- Status: Note. No audited family class is currently harmless without extra assumptions.
- Status: Note. The vector-family clarification shows that derivative-generated vector blocks do not absorb a genuine primitive family `V_i`.
- Status: Note. The smallest explicit vector-family obstruction is `V2`, with first mixed witness `EVV`.
- Status: Note. The rank-3 clarification shows that derivative-generated rank-3 blocks do not absorb a genuine primitive STF family `T_{ijk}`.
- Status: Note. The smallest explicit rank-3 obstruction is `T2`, with first mixed witness `ETT`.
- Status: Note. The rank-4 clarification still shows that trace descendants and derivative-generated rank-4 blocks do not absorb a genuine primitive STF family `Q_{ijkl}`.
- Status: Note. The smallest explicit rank-4 obstruction remains `Q2`, while the exhaustive audit now shows that the mixed cubic layer at weight `3` includes both the previously recorded `EQQ` and the omitted `EEQ`.
- Status: Note. The rank-5 clarification shows that trace descendants and derivative-generated rank-5 blocks do not absorb a genuine primitive STF family `U_{ijklm}`.
- Status: Note. The smallest explicit rank-5 obstruction is `U2`, with first mixed witness `EUU`.
- Status: Note. The stronger mixed-pattern STF-tower theorem for all genuine parity-even ranks `L \ge 3` is false at the first audited even-rank member `L = 4`.
- Status: Note. The weaker self-witness threshold theorem for all genuine parity-even STF primitive families of rank `L \ge 3` is now proven.
- Status: Note. Necessary suppression budgets are explicit for the audited classes: for rank-2 STF families the self-only lower bound is `w_X \ge 3` but the mixed-aware necessary threshold is `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness; for bare scalar families the self-only sharp threshold is `w_S \ge 5`; for derivative-only scalar families the tied-sharp threshold is `w_D \ge 3`; for genuine vector families the self-only sharp threshold is `w_V \ge 3`; for genuine rank-3 STF families the self-only sharp threshold is `w_T \ge 3`; for genuine rank-4 STF families the self-only sharp threshold is `w_Q \ge 3`; for genuine rank-5 STF families the self-only sharp threshold is `w_U \ge 3`; for genuine rank-6 STF families the self-only sharp threshold is `w_Z \ge 3`.
- Status: Note. The pairwise, triple, and quadruple composition audits already found no surviving operator beyond the baseline electric sector once the current audited thresholds were imposed across the enlarged audited family set `{R2, R0a, R0b, R1}`.
- Status: Note. The post-Rodd+ pairwise, triple, quadruple, and five-family composition audits likewise find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}`.
- Status: Note. Therefore the current thresholds are jointly sufficient for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+}` at `\Delta \le 4`, and no explicit post-Rodd+ cross-family obstruction has appeared in the audited set.
- Status: Note. The raw R3-extended survivor list is not linearly independent; its first exact quartic mixed dependence relation is `-E2T2/2 + 2 E2T2_mixed_1 - E2T2_mixed_2 + E2T2_mixed_3 = 0`.
- Status: Note. The raw R3-extended `21`-element survivor list is bookkeeping only and should not be read as a corrected basis statement; the composition audit uses the thresholded rank-3 family instead.
- Status: Note. The old raw R4-extended survivor bookkeeping is not linearly independent; its first exact quartic mixed dependence relation is `-E2Q2/2 + 2 E2Q2_mixed_1 - E2Q2_mixed_2 + E2Q2_mixed_3 = 0`.
- Status: Note. The old raw R4 bookkeeping with `22` labels, rank `19`, and nullity `3` now applies only to the pre-exhaustiveness manual subset, not to an exhaustive corrected rank-4 basis. (The exhaustive corrected rank-4 basis has survivor dimension `25`; see `../../symbolic/r4_survivor_rank_check.py` and `../../verification/rederive_rank4.py`.)
- Status: Conjectural. The raw R5-extended survivor list is not linearly independent; its first sample-extracted quartic mixed dependence relation is `E2U2 - 4 E2U2_mixed_1 + 2 E2U2_mixed_2 - 2 E2U2_mixed_3 = 0`.
- Status: Conjectural. The raw R5-extended `23`-element survivor list has sample-stable rank `19` and nullity `4`, and it is bookkeeping only rather than a corrected basis statement.
- Status: Conjectural. The raw R6-extended survivor list is not linearly independent; its first sample-revalidated mixed quartic dependence relation is `E2Z2 + 2 E2Z2_mixed_1 - 2 E2Z2_mixed_2 - 4 E2Z2_mixed_3 = 0`.
- Status: Conjectural. The raw R6-extended `22`-label family list has sample-stable rank `16` and nullity `6`, and it is bookkeeping only rather than a corrected basis statement.
- Status: Note. The irreducible family-envelope theorem now closes the parity-even nonspinning local MVP family domain on the audited scalar/vector/STF classes.
- Status: Note. The post-Reven4+ pairwise, triple, quadruple, quintuple, and six-family composition audits find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+}`.
- Status: Note. The post-Rodd5+ pairwise, triple, quadruple, quintuple, six-family, and seven-family composition audits likewise find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
- Status: Note. Therefore the current thresholds are jointly sufficient for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}` at `\Delta \le 4`, and no explicit post-Rodd5+ cross-family obstruction has appeared in the audited set.
- Status: Note. The `Reven6+` audit now resolves the next even-rank gate as a genuine obstruction class with smallest witness `Z2` and first mixed witness `EZZ`.
- Status: Note. The post-Reven6+ pairwise, triple, quadruple, quintuple, six-family, seven-family, and eight-family composition audits find no surviving operator beyond the baseline electric sector once the current audited thresholds are imposed across the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}`.
- Status: Note. Therefore the current thresholds are jointly sufficient for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}` at `\Delta \le 4`, and no explicit post-Reven6+ cross-family obstruction has appeared in the audited set.
- Status: Note. The former live bottleneck `"prove or refute irreducible family-envelope closure beyond the audited scalar/vector/STF classes"` is now resolved positively.
- Status: Note. The former live bottleneck `"prove or refute the positive finite-family collapse theorem once the irreducible scalar/vector/STF envelope is fixed"` is now resolved positively.
- Status: Note. The former live bottleneck `"replace A8 by a sharp local weight-spectrum finiteness condition"` is now resolved positively.
- Status: Note. Fixed-order operator-space non-closure is no longer a live failure mode inside the current theorem domain once local weight-spectrum finiteness below the cutoff and irreducible-envelope closure are imposed.
- Status: Note. Finite normal-form basis closure is no longer a live failure mode inside the current theorem domain once the explicit quotient rules are imposed.
- Status: Note. The positive finite-family collapse theorem now closes inside the current theorem domain.
- Status: Note. The explicit infinite low-weight STF tower is the sharp `A8` failure mode when the local spectrum condition is dropped.
- Status: Note. The former boundary-stress bottleneck `"construct or refute the smallest local nonanalytic counterexample to analytic monopole jet collapse (A5 boundary)"` is now resolved positively by the one-coordinate smooth flat monopole model `m_A(Y)=m_0+\alpha e^{-1/Y^2}\Theta(Y)`.
- Status: Note. That model preserves locality, local weight-spectrum finiteness below the cutoff, fixed-order operator finiteness, and finite scalar normal-form closure, but it destroys the finite analytic Taylor jet at the reference point.
- Status: Note. The replacement bookkeeping is non-Taylor monopole germ data and remains separate from the residual Wilson coefficients.
- Status: Note. The former boundary-stress bottleneck `"construct or refute the smallest genuinely hereditary counterexample to the local finite-family collapse theorem (A3 boundary)"` is now resolved positively by the one-coordinate causal power-law kernel model.
- Status: Note. The exponential-memory control is not the sharp `A3` escape because it collapses into a finite auxiliary-state `A4`-type extension.
- Status: Note. The exact broken layer for the genuine hereditary branch is the reduction of the monopole response to a local function `m_A(Y(\tau))`; therefore Lemmas 55 and 56 no longer apply in the local theorem form.
- Status: Note. The former boundary-stress bottleneck `"construct or refute the smallest local finite-state counterexample to the no-state theorem (A4 boundary), and determine whether a finite state-augmented collapse theorem survives"` is now resolved positively by the one-state local analytic `\chi` model of [`state-variable-counterexample-program.md`](state-variable-counterexample-program.md).
- Status: Note. The exact broken layer in that genuinely dynamical local-state branch is the original no-state reduction to `m_A(Y)` alone; Lemma 55 fails in Y-only form, and Lemma 56 fails in Y-only form as well.
- Status: Note. The adiabatic or slaved local-state control does not count as the sharp `A4` escape, because it is eliminable back into a `Y`-only monopole model.
- Status: Note. A second positive branch survives beyond the original theorem: once finitely many local state variables `\chi^a` are kept explicitly, a finite state-augmented collapse theorem survives as a finite-dimensional local state-space statement with state data kept separate from Wilson coefficients.
- Status: Counterexample candidate. The remaining explicit risks are now explicit assumption-drop or domain-escape failures: nonanalytic monopole dependence at `A5`, genuinely hereditary nonlocal monopole dependence at `A3`, failure of local weight-spectrum finiteness at `A8` through an infinite low-weight tower, mixed `A4+A5` or `A3+A4` escape beyond the sharp one-state local branch, or escape beyond the present theorem domain.

## Dependencies

- Status: Imported from prior work. [`../lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md)
- Status: Imported from prior work. [`../lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md)
- Status: Conjectural. [`../lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md)
- Status: Note. [`../lemmas/05-finite-basis-closure.md`](../lemmas/05-finite-basis-closure.md)
- Status: Conjectural. [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md)
- Status: Conjectural. [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md)
- Status: Note. [`../lemmas/08-mixed-time-derivative-audit.md`](../lemmas/08-mixed-time-derivative-audit.md)
- Status: Note. [`../lemmas/09-survivor-independence-delta4.md`](../lemmas/09-survivor-independence-delta4.md)
- Status: Note. [`../lemmas/10-magnetic-family-obstruction.md`](../lemmas/10-magnetic-family-obstruction.md)
- Status: Note. [`../lemmas/11-eb-survivor-independence-delta4.md`](../lemmas/11-eb-survivor-independence-delta4.md)
- Status: Note. [`../lemmas/12-magnetic-ordering-salvage.md`](../lemmas/12-magnetic-ordering-salvage.md)
- Status: Note. [`../lemmas/13-scalar-family-obstruction.md`](../lemmas/13-scalar-family-obstruction.md)
- Status: Note. [`../lemmas/14-derivative-only-scalar-audit.md`](../lemmas/14-derivative-only-scalar-audit.md)
- Status: Note. [`../lemmas/15-rank2-family-admission.md`](../lemmas/15-rank2-family-admission.md)
- Status: Note. [`../lemmas/16-rank0-family-admission.md`](../lemmas/16-rank0-family-admission.md)
- Status: Note. [`../lemmas/18-rank2-mixed-witness-threshold.md`](../lemmas/18-rank2-mixed-witness-threshold.md)
- Status: Note. [`../lemmas/19-rank0-mixed-witness-threshold.md`](../lemmas/19-rank0-mixed-witness-threshold.md)
- Status: Note. [`../lemmas/20-rank2-threshold-formula.md`](../lemmas/20-rank2-threshold-formula.md)
- Status: Note. [`../lemmas/21-rank0-threshold-formula.md`](../lemmas/21-rank0-threshold-formula.md)
- Status: Note. [`../lemmas/22-pairwise-composition-audit.md`](../lemmas/22-pairwise-composition-audit.md)
- Status: Note. [`../lemmas/23-three-family-composition-audit.md`](../lemmas/23-three-family-composition-audit.md)
- Status: Note. [`../lemmas/24-family-envelope-census.md`](../lemmas/24-family-envelope-census.md)
- Status: Note. [`../lemmas/25-next-unaudited-family-gate.md`](../lemmas/25-next-unaudited-family-gate.md)
- Status: Note. [`../lemmas/26-rank1-family-admission.md`](../lemmas/26-rank1-family-admission.md)
- Status: Note. [`../lemmas/27-rank1-threshold-formula.md`](../lemmas/27-rank1-threshold-formula.md)
- Status: Note. [`../lemmas/28-r1-pairwise-composition-audit.md`](../lemmas/28-r1-pairwise-composition-audit.md)
- Status: Note. [`../lemmas/29-r1-augmented-composition-audit.md`](../lemmas/29-r1-augmented-composition-audit.md)
- Status: Note. [`../lemmas/30-rank3-family-admission.md`](../lemmas/30-rank3-family-admission.md)
- Status: Note. [`../lemmas/31-rank3-threshold-formula.md`](../lemmas/31-rank3-threshold-formula.md)
- Status: Note. [`../lemmas/32-r3-pairwise-composition-audit.md`](../lemmas/32-r3-pairwise-composition-audit.md)
- Status: Note. [`../lemmas/33-r3-augmented-composition-audit.md`](../lemmas/33-r3-augmented-composition-audit.md)
- Status: Note. [`../lemmas/34-rank4-family-admission.md`](../lemmas/34-rank4-family-admission.md)
- Status: Note. [`../lemmas/35-rank4-threshold-formula.md`](../lemmas/35-rank4-threshold-formula.md)
- Status: Note. [`../lemmas/36-r4-pairwise-composition-audit.md`](../lemmas/36-r4-pairwise-composition-audit.md)
- Status: Note. [`../lemmas/37-r4-augmented-composition-audit.md`](../lemmas/37-r4-augmented-composition-audit.md)
- Status: Note. [`../lemmas/38-rank5-family-admission.md`](../lemmas/38-rank5-family-admission.md)
- Status: Note. [`../lemmas/39-rank5-threshold-formula.md`](../lemmas/39-rank5-threshold-formula.md)
- Status: Note. [`../lemmas/40-r5-pairwise-composition-audit.md`](../lemmas/40-r5-pairwise-composition-audit.md)
- Status: Note. [`../lemmas/41-r5-augmented-composition-audit.md`](../lemmas/41-r5-augmented-composition-audit.md)
- Status: Note. [`../lemmas/48-stf-self-witness-theorem.md`](../lemmas/48-stf-self-witness-theorem.md)
- Status: Note. [`../lemmas/49-stf-threshold-theorem.md`](../lemmas/49-stf-threshold-theorem.md)
- Status: Note. [`../lemmas/50-cartesian-irrep-reduction.md`](../lemmas/50-cartesian-irrep-reduction.md)
- Status: Note. [`../lemmas/51-trace-descendant-absorption.md`](../lemmas/51-trace-descendant-absorption.md)
- Status: Note. [`../lemmas/52-mixed-symmetry-family-gate.md`](../lemmas/52-mixed-symmetry-family-gate.md)
- Status: Note. [`../lemmas/53-fixed-order-operator-finiteness.md`](../lemmas/53-fixed-order-operator-finiteness.md)
- Status: Note. [`../lemmas/54-normal-form-basis-closure.md`](../lemmas/54-normal-form-basis-closure.md)
- Status: Note. [`../lemmas/55-monopole-jet-collapse.md`](../lemmas/55-monopole-jet-collapse.md)
- Status: Note. [`../lemmas/56-sensitivity-vs-wilson-split.md`](../lemmas/56-sensitivity-vs-wilson-split.md)
- Status: Note. [`../lemmas/57-local-finite-but-nonanalytic.md`](../lemmas/57-local-finite-but-nonanalytic.md)
- Status: Note. [`../lemmas/58-nonanalytic-jet-failure.md`](../lemmas/58-nonanalytic-jet-failure.md)
- Status: Note. [`../lemmas/59-nonanalytic-sensitivity-replacement.md`](../lemmas/59-nonanalytic-sensitivity-replacement.md)
- Status: Note. [`../lemmas/60-hereditary-kernel-construction.md`](../lemmas/60-hereditary-kernel-construction.md)
- Status: Note. [`../lemmas/61-markovianizable-vs-genuine-hereditary.md`](../lemmas/61-markovianizable-vs-genuine-hereditary.md)
- Status: Note. [`../lemmas/62-hereditary-collapse-failure.md`](../lemmas/62-hereditary-collapse-failure.md)
- Status: Note. [`../lemmas/63-hereditary-replacement-bookkeeping.md`](../lemmas/63-hereditary-replacement-bookkeeping.md)
- Status: Note. [`../lemmas/64-local-finite-state-construction.md`](../lemmas/64-local-finite-state-construction.md)
- Status: Note. [`../lemmas/65-adiabatic-vs-dynamical-state-boundary.md`](../lemmas/65-adiabatic-vs-dynamical-state-boundary.md)
- Status: Note. [`../lemmas/66-a4-collapse-failure.md`](../lemmas/66-a4-collapse-failure.md)
- Status: Note. [`../lemmas/67-state-augmented-collapse.md`](../lemmas/67-state-augmented-collapse.md)
- Status: Note. [`../lemmas/68-state-vs-wilson-split.md`](../lemmas/68-state-vs-wilson-split.md)
- Status: Note. [`../lemmas/69-local-weight-spectrum-finiteness.md`](../lemmas/69-local-weight-spectrum-finiteness.md)
- Status: Note. [`../lemmas/70-infinite-low-weight-tower-counterexample.md`](../lemmas/70-infinite-low-weight-tower-counterexample.md)
- Status: Note. [`../lemmas/71-a8-sharpness.md`](../lemmas/71-a8-sharpness.md)
- Status: Conjectural. [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md)
- Status: Conjectural. [`eb-conditional-collapse.md`](eb-conditional-collapse.md)
- Status: Note. [`finite-family-collapse-theorem.md`](finite-family-collapse-theorem.md)
- Status: Note. [`collapse-bridge-status.md`](collapse-bridge-status.md)
- Status: Note. [`nonanalytic-counterexample-program.md`](nonanalytic-counterexample-program.md)
- Status: Note. [`nonanalytic-model-table.md`](nonanalytic-model-table.md)
- Status: Conjectural. [`power-counting.md`](power-counting.md)
- Status: Note. [`family-admission-no-go.md`](family-admission-no-go.md)
- Status: Note. [`family-admission-theorem.md`](family-admission-theorem.md)
- Status: Note. [`family-envelope-theorem.md`](family-envelope-theorem.md)
- Status: Note. [`family-envelope-table.md`](family-envelope-table.md)
- Status: Note. [`irreducible-envelope-theorem.md`](irreducible-envelope-theorem.md)
- Status: Note. [`mixed-symmetry-risk-register.md`](mixed-symmetry-risk-register.md)
- Status: Note. [`rank4-family-ordering.md`](rank4-family-ordering.md)
- Status: Note. [`rank5-family-ordering.md`](rank5-family-ordering.md)
- Status: Note. [`rank6-family-ordering.md`](rank6-family-ordering.md)
- Status: Note. [`stf-self-witness-theorem.md`](stf-self-witness-theorem.md)
- Status: Note. [`stf-mixed-pattern-classification.md`](stf-mixed-pattern-classification.md)
- Status: Note. [`rank3-family-ordering.md`](rank3-family-ordering.md)
- Status: Note. [`vector-family-ordering.md`](vector-family-ordering.md)
- Status: Note. [`mixed-witness-thresholds.md`](mixed-witness-thresholds.md)
- Status: Note. [`threshold-formulas.md`](threshold-formulas.md)
- Status: Note. [`budget-composition-theorem.md`](budget-composition-theorem.md)
- Status: Note. [`composition-status.md`](composition-status.md)
- Status: Note. [`suppression-budget-theorem.md`](suppression-budget-theorem.md)
- Status: Note. [`sharp-threshold-status.md`](sharp-threshold-status.md)
- Status: Note. [`broad-collapse-reformulation.md`](broad-collapse-reformulation.md)
- Status: Note. [`family-class-table.md`](family-class-table.md)
- Status: Note. [`primitive-catalog.md`](primitive-catalog.md)
- Status: Note. [`primitive-set-adequacy.md`](primitive-set-adequacy.md)
- Status: Note. [`magnetic-family-ordering.md`](magnetic-family-ordering.md)
- Status: Note. [`scalar-family-ordering.md`](scalar-family-ordering.md)
- Status: Note. [`reduction-rules.md`](reduction-rules.md)
- Status: Note. [`../symbolic/witness_threshold_map.py`](../symbolic/witness_threshold_map.py)
- Status: Note. [`../symbolic/mixed_witness_map.py`](../symbolic/mixed_witness_map.py)
- Status: Note. [`../symbolic/threshold_formula_check.py`](../symbolic/threshold_formula_check.py)
- Status: Note. [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py)
- Status: Note. [`../symbolic/irrep_family_census.py`](../symbolic/irrep_family_census.py)
- Status: Note. [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py)
- Status: Note. [`../symbolic/fixed_family_operator_count.py`](../symbolic/fixed_family_operator_count.py)
- Status: Note. [`../symbolic/hereditary_kernel_demo.py`](../symbolic/hereditary_kernel_demo.py)
- Status: Note. [`../symbolic/nonanalytic_jet_demo.py`](../symbolic/nonanalytic_jet_demo.py)
- Status: Note. [`../symbolic/stateful_counterexample_demo.py`](../symbolic/stateful_counterexample_demo.py)
- Status: Note. [`../symbolic/weight_spectrum_demo.py`](../symbolic/weight_spectrum_demo.py)
- Status: Note. [`../symbolic/r1_sector_delta4.py`](../symbolic/r1_sector_delta4.py)
- Status: Note. [`../symbolic/r1_survivor_rank_check.py`](../symbolic/r1_survivor_rank_check.py)
- Status: Note. [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py)
- Status: Note. [`../symbolic/r3_survivor_rank_check.py`](../symbolic/r3_survivor_rank_check.py)
- Status: Note. [`../symbolic/r4_sector_delta4.py`](../symbolic/r4_sector_delta4.py)
- Status: Note. [`../symbolic/r4_survivor_rank_check.py`](../symbolic/r4_survivor_rank_check.py)
- Status: Note. [`../symbolic/r5_sector_delta4.py`](../symbolic/r5_sector_delta4.py)
- Status: Note. [`../symbolic/r5_survivor_rank_check.py`](../symbolic/r5_survivor_rank_check.py)
- Status: Note. [`../symbolic/r6_sector_delta4.py`](../symbolic/r6_sector_delta4.py)
- Status: Note. [`../symbolic/r6_survivor_rank_check.py`](../symbolic/r6_survivor_rank_check.py)
- Status: Note. [`../symbolic/stf_self_witness_check.py`](../symbolic/stf_self_witness_check.py)

## Failure Triggers

- Status: Note. An explicit orbital-timescale `chi_A` state invalidates the original no-state reduction step; the canonical one-state local analytic model in [`state-variable-counterexample-program.md`](state-variable-counterexample-program.md) now realizes this trigger explicitly.
- Status: Counterexample candidate. Nonanalytic activation invalidates the analytic collapse step.
- Status: Counterexample candidate. Genuinely hereditary kernels invalidate the reduction of the monopole response to a local function of the instantaneous normal-form coordinates even when the local operator basis is already finite.
- Status: Note. An infinite low-weight family tower invalidates the fixed-order closure theorem candidate already at the candidate-operator-finiteness layer.
- Status: Counterexample candidate. Any additional `\Delta \le 4` scalar contraction from the exact current electric-only primitive set that evades the explicit enumeration/reduction pipeline obstructs the electric-only exact-current-set theorem.
- Status: Counterexample candidate. Any further physically admissible primitive family that yields new survivors beyond the corrected explicit audited family catalogs reinforces the uniqueness no-go branch.
- Status: Counterexample candidate. Any further physically admissible primitive family obstructs the positive finite-family collapse branch only if it also destroys fixed-order finiteness or the family-conditioned reduction path.
