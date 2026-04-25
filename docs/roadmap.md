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

## Near-Term Sequence

1. Status: Conjectural. Decide whether to pursue the corrected `E/B` theorem route as the main line or to impose a stronger magnetic-ordering salvage assumption for an electric-only theorem.
2. Status: Conjectural. If the corrected `E/B` route is accepted, attack the next primitive-family extension and the conditional collapse step rather than returning to electric-only salvage by default.
3. Status: Conjectural. Decide whether any further mixed-trace or gradient identity should be introduced explicitly, rather than silently collapsing enlarged-sector operators.
4. Status: Counterexample candidate. Promote the dominant loophole class only if it produces a genuine failure of the corrected `E/B` basis or of the next primitive-family attack.
5. Status: Conjectural. Revisit the scalar `s_A` corollary only after the higher-dimensional sensitivity manifold statement is sharp.

## Dynamic-Loophole Addendum

Status: Counterexample candidate. This addendum is separate from the numbered
static-family milestones above. It attacks A4 directly rather than adding more
static primitive families.

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| D1 | Counterexample candidate | [`../counterexamples/chi-state/README.md`](../counterexamples/chi-state/README.md) | A one-state `chi` loophole is written as an explicit A4 violation. |
| D2 | Proven | [`nonlinear-comparator-audit.md`](nonlinear-comparator-audit.md) | Single-sideband static nonlinear mimicry and the shared-pole obstruction are separated. |
| D3 | Proven | [`sideband-observable-targets.md`](sideband-observable-targets.md) | Sideband existence, shared phase law, ratio target, projection collapse, and nuisance collapse are classified. |
| D4 | Proven | [`../symbolic/nonlinear_comparator_audit.py`](../symbolic/nonlinear_comparator_audit.py) | The shared `tau_chi` phase law and finite-polynomial interpolation residual are checked symbolically. |
| D5 | Counterexample candidate | [`shared-tau-ratio-theorem.md`](shared-tau-ratio-theorem.md) | The generated-line and linear-line ratio law is stated as a theorem/no-go target. |
| D6 | Proven | [`shared-tau-observable-targets.md`](shared-tau-observable-targets.md) | Ratio success and failure conditions are separated from sideband existence. |
| D7 | Proven | [`../symbolic/shared_tau_ratio_audit.py`](../symbolic/shared_tau_ratio_audit.py) | The exact shared-tau diagnostic residual and finite interpolation budgets are checked symbolically. |
| D8 | Proven | [`sample-budget-theorem.md`](sample-budget-theorem.md) | The minimum sample count is stated in terms of linear samples, sideband pairs, and projection nuisance parameters. |
| D9 | Proven | [`orbital-budget-case.md`](orbital-budget-case.md) and [`two-tone-budget-case.md`](two-tone-budget-case.md) | The current orbital and two-tone dictionaries are classified against realistic comparator classes. |
| D10 | Proven | [`../symbolic/sample_budget_audit.py`](../symbolic/sample_budget_audit.py) | The sample-budget theorem and case classifications are checked symbolically. |
| D11 | Counterexample candidate | [`richer-orbital-forcing-budget-audit.md`](richer-orbital-forcing-budget-audit.md) | The bounded exact-in-e/higher-harmonic rescue route is classified by `(H,R)` harmonic budgets. |
| D12 | Proven | [`../symbolic/orbital_harmonic_budget_audit.py`](../symbolic/orbital_harmonic_budget_audit.py) | The minimum richer-orbital harmonic counts are computed for the comparator grid. |
| D13 | Counterexample candidate | [`second-order-internal-mode.md`](second-order-internal-mode.md) | The next dynamic loophole model is defined with resonance, phase wrapping, and first-order collapse limits. |
| D14 | Proven | [`../symbolic/second_order_mode_response.py`](../symbolic/second_order_mode_response.py) | The second-order transfer, low-frequency expansion, and resonance condition are checked symbolically. |
| D15 | Counterexample candidate | [`second-order-projection-audit.md`](second-order-projection-audit.md) | Acceleration-like and range-like readouts are audited for resonance survival under finite shared projection nuisance. |
| D16 | Counterexample candidate | [`resonant-comparator-theorem.md`](resonant-comparator-theorem.md) | The second-order line-shape budget theorem is stated with resonance bracketing and projection nuisance counts. |
| D17 | Proven | [`../symbolic/second_order_projection_audit.py`](../symbolic/second_order_projection_audit.py) and [`../symbolic/resonant_comparator_audit.py`](../symbolic/resonant_comparator_audit.py) | Projection residuals and resonance-aware comparator budgets are checked symbolically. |
| D18 | Counterexample candidate | [`exact-in-e-resonant-forcing-budget-audit.md`](exact-in-e-resonant-forcing-budget-audit.md) | Exact-in-e orbital harmonics are audited as a sample-provisioning dictionary for the second-order resonance target. |
| D19 | Proven | [`../symbolic/exact_in_e_resonant_forcing_audit.py`](../symbolic/exact_in_e_resonant_forcing_audit.py) | The minimum usable harmonic cutoff `H_min` is computed from comparator budget and resonance bracketing. |
| D20 | Counterexample candidate | [`amplitude-weighted-resonant-design-theorem.md`](amplitude-weighted-resonant-design-theorem.md) | Exact-in-e harmonics are weighted by forcing coefficient, second-order transfer, projection, and a usable-sample cutoff. |
| D21 | Proven | [`../symbolic/amplitude_weighted_resonant_design.py`](../symbolic/amplitude_weighted_resonant_design.py) | The amplitude-weighted usable-sample count is checked for low- and moderate-eccentricity designs. |
| D22 | Counterexample candidate | [`physical-detectability-map.md`](physical-detectability-map.md) | The relative amplitude cutoff is replaced by a parameterized SNR criterion and minimum physical scale boundary. |
| D23 | Proven | [`../symbolic/physical_detectability_map.py`](../symbolic/physical_detectability_map.py) | The SNR-budget classifications and minimum physical scale examples are checked symbolically. |
| D24 | Counterexample candidate | [`nonlinear-second-order-mode.md`](nonlinear-second-order-mode.md) | Minimal nonlinear second-order drive/readout is defined and the shared quadratic-denominator target is separated from sideband existence. |
| D25 | Proven | [`../symbolic/nonlinear_second_order_mode.py`](../symbolic/nonlinear_second_order_mode.py) | The nonlinear second-order transfers, shared-denominator residual, and generated-sideband interpolation boundary are checked symbolically. |
| D26 | Counterexample candidate | [`nonlinear-second-order-detectability.md`](nonlinear-second-order-detectability.md) | Exact-in-e nonlinear generated lines are audited against relative and SNR usable-sample criteria. |
| D27 | Proven | [`../symbolic/nonlinear_second_order_detectability.py`](../symbolic/nonlinear_second_order_detectability.py) | Generated-line sample count, amplitude threshold, projection nuisance, and minimum generated-amplitude scales are checked symbolically. |
| D28 | Counterexample candidate | [`component-separability-theorem.md`](component-separability-theorem.md) | Overlapping linear and nonlinear generated harmonic components are tested for local identifiability by joint-model rank. |
| D29 | Proven | [`../symbolic/component_separability_audit.py`](../symbolic/component_separability_audit.py) | The `Q_beta` rank contribution, full shared-parameter rank, and budget-ready minimum cutoffs are checked numerically. |
| D30 | Counterexample candidate | [`nonlinear-robustness-map.md`](nonlinear-robustness-map.md) | The nonlinear second-order positive region is stress-tested over a finite dimensionless design grid. |
| D31 | Proven | [`../symbolic/nonlinear_robustness_map.py`](../symbolic/nonlinear_robustness_map.py) | Robust-positive, underdesigned, and degenerate grid regions are counted by projection mode. |
| D32 | Counterexample candidate | [`analytic-phase-boundary-theorem.md`](analytic-phase-boundary-theorem.md) | The nonlinear positive region is rewritten as count, bracket, amplitude-usability, and rank-boundary conditions. |
| D33 | Proven | [`../symbolic/analytic_phase_boundary.py`](../symbolic/analytic_phase_boundary.py) | The analytic boundary rules reproduce the default robustness-map verdicts and record the `Q_beta` rank minor. |

## Explicit Non-Goals For M1

- Status: Proven. No new LLR, MLRS, PEP, Nutimo, or TOA pipeline work belongs in M1.
- Status: Proven. No clock-sector theorem is attempted in M1.
- Status: Proven. No build or runtime chores count as theorem progress.

## Classification Rule

- `theorem progress`: theorem draft, lemma draft, assumption cleanup, symbolic reduction check
- `loophole progress`: explicit counterexample model, exact failed proof step, minimal violated assumption
