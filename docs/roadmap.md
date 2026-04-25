# Roadmap

Status: Conjectural. The repository roadmap is organized around theorem progress first and loophole progress second.

## M1: Free-Fall Sensitivity-Collapse Scaffold

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M1.1 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) | Precise theorem statement and proof skeleton exist. |
| M1.2 | Imported from prior work | [`lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md) | Imported obstruction to literal internal self-unobservability is recorded. |
| M1.3 | Imported from prior work | [`lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md) | Imported COM decoupling lemma is recorded. |
| M1.4 | Conjectural | [`lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md) | Generic local worldline action note exists. |
| M1.5 | Counterexample candidate | [`counterexamples/chi-state/README.md`](../counterexamples/chi-state/README.md) | At least one explicit loophole model is written down. |
| M1.6 | Proven | [`symbolic/worldline_expand.py`](../symbolic/worldline_expand.py) and [`symbolic/checks/test_symbolic.py`](../symbolic/checks/test_symbolic.py) | Symbolic scaffold runs without error. |

## M2: Non-Circular Fixed-Order Theorem Candidate

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M2.1 | Conjectural | [`docs/conditional-collapse-lemma.md`](conditional-collapse-lemma.md) | Conditional collapse is isolated from basis closure. |
| M2.2 | Conjectural | [`docs/power-counting.md`](power-counting.md) | Fixed-order counting rule is explicit and does not assume finite basis closure. |
| M2.3 | Proven | [`lemmas/05-finite-basis-closure.md`](../lemmas/05-finite-basis-closure.md) | Abstract candidate-set finiteness is separated from the physical completeness burden. |
| M2.4 | Proven | [`../symbolic/enumerate_basis.py`](../symbolic/enumerate_basis.py) | The chosen-order candidate monomial list enumerates without error. |
| M2.5 | Counterexample candidate | [`../counterexamples/chi-state/README.md`](../counterexamples/chi-state/README.md) | The smallest loophole model includes an adiabatic-collapse analysis. |

## M3: Delta<=4 Normal-Form Completeness Attack

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M3.1 | Proven | [`docs/primitive-catalog.md`](primitive-catalog.md) | Exact `Delta<=4` primitive and candidate operator content is fixed. |
| M3.2 | Proven | [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) | Total-derivative, lower-order EOM, and algebraic reductions run without error. |
| M3.3 | Proven | [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md) | The old five-element target is recognized as incomplete rather than silently assumed. |

## M4: Contraction-Level Exhaustiveness Audit

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M4.1 | Proven | [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) | All parity-even scalar contractions from the exact current primitive blocks through `Delta<=4` are enumerated. |
| M4.2 | Proven | [`docs/reduction-rules.md`](reduction-rules.md) | Every allowed reduction rule is listed explicitly. |
| M4.3 | Conjectural | [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md) | The gradient sector is audited and the surviving basis is identified explicitly. |
| M4.4 | Proven | [`../lemmas/08-mixed-time-derivative-audit.md`](../lemmas/08-mixed-time-derivative-audit.md) | Mixed `E-E-D_tau E` terms are audited explicitly. |
| M4.5 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) | The corrected seven-element `Delta<=4` normal-form path is stated explicitly. |

## M5: Independence And Primitive-Set Adequacy

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M5.1 | Proven | [`../lemmas/09-survivor-independence-delta4.md`](../lemmas/09-survivor-independence-delta4.md) | The corrected seven survivors are checked for linear independence rather than inferred from non-reducibility. |
| M5.2 | Proven | [`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py) | Exact polynomial coefficient rank confirms the current seven-scalar basis is linearly independent. |
| M5.3 | Proven | [`primitive-set-adequacy.md`](primitive-set-adequacy.md) | The exact-current-set theorem is separated sharply from any physically justified minimal-sector claim. |
| M5.4 | Proven | [`../symbolic/primitive_family_attack.py`](../symbolic/primitive_family_attack.py) | One physically motivated primitive-family extension is audited explicitly. |
| M5.5 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) and [`failure-ledger.md`](failure-ledger.md) | The main remaining risk is written uniformly as primitive-set adequacy, not survivor dependence. |

## M6: Magnetic-Family Verdict

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M6.1 | Proven | [`docs/magnetic-family-ordering.md`](magnetic-family-ordering.md) | The magnetic-ordering question and suppression options are written explicitly rather than tacitly. |
| M6.2 | Proven | [`../lemmas/10-magnetic-family-obstruction.md`](../lemmas/10-magnetic-family-obstruction.md) | The exact object obstructed by `B2` is stated sharply. |
| M6.3 | Proven | [`../symbolic/eb_sector_delta4.py`](../symbolic/eb_sector_delta4.py) | The `E/B`-expanded `Delta<=4` survivor list is computed under the current rules. |
| M6.4 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) | Electric-only, stronger minimal-sector, and broader collapse statements are separated cleanly. |
| M6.5 | Proven | [`failure-ledger.md`](failure-ledger.md) | The live bottleneck is identified as magnetic-family ordering rather than hidden basis failure. |

## M7: `E/B`-Expanded Collapse Split

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M7.1 | Proven | [`../symbolic/eb_survivor_rank_check.py`](../symbolic/eb_survivor_rank_check.py) | The raw `E/B` survivor candidate set is checked for independence and the first exact relation is extracted. |
| M7.2 | Proven | [`../lemmas/11-eb-survivor-independence-delta4.md`](../lemmas/11-eb-survivor-independence-delta4.md) | The corrected `E/B` basis is identified and shown linearly independent. |
| M7.3 | Conjectural | [`eb-conditional-collapse.md`](eb-conditional-collapse.md) | The conditional collapse step is restated on the corrected `E/B` basis. |
| M7.4 | Proven | [`../lemmas/12-magnetic-ordering-salvage.md`](../lemmas/12-magnetic-ordering-salvage.md) | The exact strength and cost of electric-only salvage assumptions are written explicitly. |
| M7.5 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) and [`failure-ledger.md`](failure-ledger.md) | The repo now splits electric-only salvage from the corrected `E/B` theorem route. |

## M8: Family-Admission No-Go And Necessary Budgets

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M8.1 | Proven | [`family-admission-theorem.md`](family-admission-theorem.md) | The negative uniqueness branch is promoted from isolated audits into a class-limited family-admission theorem. |
| M8.2 | Proven | [`../lemmas/15-rank2-family-admission.md`](../lemmas/15-rank2-family-admission.md) and [`../lemmas/16-rank0-family-admission.md`](../lemmas/16-rank0-family-admission.md) | Rank-2 STF and rank-0 scalar family classes are abstracted from the explicit audits. |
| M8.3 | Proven | [`../symbolic/family_witness_map.py`](../symbolic/family_witness_map.py) | The current audited family classes and their canonical witnesses are summarized deterministically. |
| M8.4 | Proven | [`suppression-budget-theorem.md`](suppression-budget-theorem.md), [`family-class-table.md`](family-class-table.md), and [`../symbolic/witness_threshold_map.py`](../symbolic/witness_threshold_map.py) | Necessary class-limited suppression budgets are explicit at `\Delta_{\max} = 4`. |
| M8.5 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) and [`failure-ledger.md`](failure-ledger.md) | The main positive theorem target is finite-family fixed-order collapse rather than minimal-sector uniqueness. |

## M9: Mixed-Witness Threshold Classification

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M9.1 | Proven | [`mixed-witness-thresholds.md`](mixed-witness-thresholds.md) | Self witnesses, mixed witnesses, and `W_{\min}` are defined explicitly for the audited family classes. |
| M9.2 | Proven | [`../lemmas/18-rank2-mixed-witness-threshold.md`](../lemmas/18-rank2-mixed-witness-threshold.md) | The audited rank-2 STF class is identified as tied in the unsuppressed audit but mixed-aware at the threshold level. |
| M9.3 | Proven | [`../lemmas/19-rank0-mixed-witness-threshold.md`](../lemmas/19-rank0-mixed-witness-threshold.md) | The audited rank-0 bare-scalar and derivative-only subclasses are split into self-dominated versus tied sharp thresholds. |
| M9.4 | Proven | [`../symbolic/mixed_witness_map.py`](../symbolic/mixed_witness_map.py) | A deterministic machine-readable mixed-witness threshold table runs without error. |
| M9.5 | Proven | [`sharp-threshold-status.md`](sharp-threshold-status.md) and [`family-class-table.md`](family-class-table.md) | Each audited family class is labeled as self-only, mixed-aware, or tied-sharp. |
| M9.6 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) and [`failure-ledger.md`](failure-ledger.md) | The live negative-branch bottleneck is stated as sharp mixed-aware threshold classification beyond mere witness existence. |

## M10: Formula-Level Threshold Classification

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M10.1 | Proven | [`threshold-formulas.md`](threshold-formulas.md) | Audited family classes are lifted into formula-level `W_{\mathrm{self}}`, `W_{\mathrm{mix}}`, and `W_{\min}` notation. |
| M10.2 | Proven | [`../lemmas/20-rank2-threshold-formula.md`](../lemmas/20-rank2-threshold-formula.md) | The rank-2 STF class now distinguishes the self-only lower bound from the mixed-aware necessary threshold. |
| M10.3 | Proven | [`../lemmas/21-rank0-threshold-formula.md`](../lemmas/21-rank0-threshold-formula.md) | The bare-scalar and derivative-only subclasses are classified as self-only versus tied-sharp at `\Delta_{\max}=4`. |
| M10.4 | Proven | [`../symbolic/threshold_formula_check.py`](../symbolic/threshold_formula_check.py) | The audited formula layer reproduces the current `\Delta \le 4` witness map deterministically. |
| M10.5 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) and [`family-admission-theorem.md`](family-admission-theorem.md) | The live negative-branch bottleneck is written uniformly as sharp mixed-aware threshold classification beyond mere witness existence. |

## M11: Audited-Set Budget Composition

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M11.1 | Proven | [`budget-composition-theorem.md`](budget-composition-theorem.md) | The joint-sufficiency question is stated explicitly for the already audited family classes. |
| M11.2 | Proven | [`../lemmas/22-pairwise-composition-audit.md`](../lemmas/22-pairwise-composition-audit.md) | All three pairwise combinations are audited separately. |
| M11.3 | Proven | [`../lemmas/23-three-family-composition-audit.md`](../lemmas/23-three-family-composition-audit.md) | The simultaneous three-family admission is audited separately from the pairwise cases. |
| M11.4 | Proven | [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) | The current audited thresholds are tested jointly against all relevant `\Delta \le 4` scalar candidates. |
| M11.5 | Proven | [`composition-status.md`](composition-status.md), [`docs/theorem-A-freefall.md`](theorem-A-freefall.md), and [`failure-ledger.md`](failure-ledger.md) | The repo records whether the audited thresholds are jointly sufficient or whether an explicit cross-family obstruction survives. |

## M12: Family-Envelope Census

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M12.1 | Proven | [`family-envelope-theorem.md`](family-envelope-theorem.md) | Audited-set sufficiency is separated from MVP-envelope completeness. |
| M12.2 | Proven | [`family-envelope-table.md`](family-envelope-table.md) | Covered, excluded, and still unaudited family classes are listed explicitly. |
| M12.3 | Proven | [`../lemmas/24-family-envelope-census.md`](../lemmas/24-family-envelope-census.md) | The current MVP family-class census is derived under the active assumptions. |
| M12.4 | Proven | [`../symbolic/family_envelope_census.py`](../symbolic/family_envelope_census.py) | The family-envelope census runs and reports the current next unaudited class deterministically. |
| M12.5 | Proven | [`../lemmas/25-next-unaudited-family-gate.md`](../lemmas/25-next-unaudited-family-gate.md), [`theorem-A-freefall.md`](theorem-A-freefall.md), and [`failure-ledger.md`](failure-ledger.md) | The repo states explicitly whether the audited-set result upgrades to envelope sufficiency or remains audited-set only. |

## M13: Rank-1 Vector Family Gate

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M13.1 | Proven | [`vector-family-ordering.md`](vector-family-ordering.md) | A genuine primitive vector family is distinguished from derivative-generated vector blocks already belonging to audited families. |
| M13.2 | Proven | [`../lemmas/26-rank1-family-admission.md`](../lemmas/26-rank1-family-admission.md) | The first self and mixed witnesses for a genuine rank-1 vector family are identified explicitly. |
| M13.3 | Proven | [`../symbolic/r1_sector_delta4.py`](../symbolic/r1_sector_delta4.py) | The R1-extended `\Delta \le 4` scalar sector enumerates without error and reports the smallest new witness. |
| M13.4 | Proven | [`../symbolic/r1_survivor_rank_check.py`](../symbolic/r1_survivor_rank_check.py) | The corrected R1-extended survivor list is checked for linear dependence whenever more than one new survivor appears. |
| M13.5 | Proven | [`../lemmas/27-rank1-threshold-formula.md`](../lemmas/27-rank1-threshold-formula.md) | The vector-family threshold is classified as self-only, mixed-aware, or tied-sharp. |
| M13.6 | Proven | [`family-envelope-table.md`](family-envelope-table.md), [`theorem-A-freefall.md`](theorem-A-freefall.md), and [`failure-ledger.md`](failure-ledger.md) | The R1 row is resolved explicitly and the next live envelope gate is stated uniformly. |

## M14: Rank-3 Family Gate

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M14.1 | Proven | [`rank3-family-ordering.md`](rank3-family-ordering.md) | A genuine primitive rank-3 family is distinguished from derivative-generated rank-3 descendants already attached to audited families. |
| M14.2 | Proven | [`../lemmas/30-rank3-family-admission.md`](../lemmas/30-rank3-family-admission.md) | The first self and mixed witnesses for the genuine rank-3 family are identified explicitly. |
| M14.3 | Proven | [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py) | The Rodd+-extended `\Delta \le 4` scalar sector enumerates without error and reports the smallest new witness. |
| M14.4 | Proven | [`../symbolic/r3_survivor_rank_check.py`](../symbolic/r3_survivor_rank_check.py) | The raw Rodd+-extended survivor list is checked for linear dependence whenever multiple new survivors appear. |
| M14.5 | Proven | [`../lemmas/31-rank3-threshold-formula.md`](../lemmas/31-rank3-threshold-formula.md) | The rank-3 family threshold is classified as self-only, mixed-aware, or tied-sharp. |
| M14.6 | Proven | [`family-envelope-table.md`](family-envelope-table.md), [`theorem-A-freefall.md`](theorem-A-freefall.md), and [`failure-ledger.md`](failure-ledger.md) | The Rodd+ row is resolved explicitly and the post-audit live bottleneck is stated uniformly. |

## Near-Term Sequence

1. Status: Proven. The enlarged audited-set composition question is now closed positively for the audited family classes `{R2, R0a, R0b, R1}`.
2. Status: Proven. The rank-3 family gate `Rodd+` is now resolved as a genuine new obstruction class rather than an absorbed derivative block.
3. Status: Proven. Re-close audited-set composition for the further enlarged family set that now includes `Rodd+`.
4. Status: Proven. The rank-4 family gate `Reven4+` is now resolved as a genuine new obstruction class rather than an absorbed trace descendant or derivative block.
5. Status: Proven. Audited-set composition is now re-closed for the further enlarged family set that includes `Reven4+`.
6. Status: Conjectural. The live envelope gate can now move to `Rodd5+`, the next smallest unaudited family class.
7. Status: Conjectural. Push the positive branch from audited finite examples toward an abstract finite-family fixed-order collapse theorem.
8. Status: Counterexample candidate. Promote the dominant loophole class only if it produces a genuine failure of the corrected finite-family collapse candidate rather than only another uniqueness obstruction.
9. Status: Conjectural. Revisit the scalar `s_A` corollary only after the higher-dimensional sensitivity manifold statement is sharp.

## Explicit Non-Goals For M1

- Status: Proven. No new LLR, MLRS, PEP, Nutimo, or TOA pipeline work belongs in M1.
- Status: Proven. No clock-sector theorem is attempted in M1.
- Status: Proven. No build or runtime chores count as theorem progress.

## Classification Rule

- `theorem progress`: theorem draft, lemma draft, assumption cleanup, symbolic reduction check
- `loophole progress`: explicit counterexample model, exact failed proof step, minimal violated assumption
