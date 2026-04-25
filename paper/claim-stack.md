# Dynamic Loophole Claim Stack

This file maps every main paper claim to an existing repository artifact. It
does not introduce new science.

## Strongest Thesis

One-state relaxation is a control/no-go branch under the current forcing and
projection dictionaries. Nonlinear second-order internal visibility is the
live counterexample candidate: exact-in-e eccentric forcing plus finite shared
projection nuisance can yield budget-breaking, locally separable
generated-line observables with robustness regions and analytic phase
boundaries. This is a parameterized design theorem, not an instrument
forecast.

Primary source:
[`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md)

## Proven Claims

| Claim | Supporting artifact | Script check | Paper placement |
| --- | --- | --- | --- |
| One-frequency phase lag or one sideband is not sufficient novelty once local derivative or static nonlinear comparators are admitted. | [`../docs/nonlinear-comparator-audit.md`](../docs/nonlinear-comparator-audit.md), [`../docs/sideband-observable-targets.md`](../docs/sideband-observable-targets.md) | [`../symbolic/nonlinear_comparator_audit.py`](../symbolic/nonlinear_comparator_audit.py) | Sections 3 and 5 |
| Finite dynamic samples must exceed the finite shared comparator and projection nuisance budgets. | [`../docs/sample-budget-theorem.md`](../docs/sample-budget-theorem.md) | [`../symbolic/sample_budget_audit.py`](../symbolic/sample_budget_audit.py) | Section 2 |
| The current one-state orbital and two-tone dictionaries are underbudget against realistic comparator classes. | [`../docs/orbital-budget-case.md`](../docs/orbital-budget-case.md), [`../docs/two-tone-budget-case.md`](../docs/two-tone-budget-case.md) | [`../symbolic/sample_budget_audit.py`](../symbolic/sample_budget_audit.py) | Section 3 |
| A resonance or phase-wrap observation is not sufficient unless sample count and bracketing conditions are both satisfied. | [`../docs/resonant-comparator-theorem.md`](../docs/resonant-comparator-theorem.md) | [`../symbolic/resonant_comparator_audit.py`](../symbolic/resonant_comparator_audit.py) | Section 4 |
| Generated-line novelty requires at least `M+2+K_Lambda` usable generated samples for a degree-`M` static nonlinear comparator with shared projection nuisance. | [`../docs/nonlinear-second-order-detectability.md`](../docs/nonlinear-second-order-detectability.md) | [`../symbolic/nonlinear_second_order_detectability.py`](../symbolic/nonlinear_second_order_detectability.py) | Section 7 |
| Component separability is necessary: the generated amplitude column must add local Jacobian rank. | [`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md) | [`../symbolic/component_separability_audit.py`](../symbolic/component_separability_audit.py) | Section 8 |
| The analytic phase-boundary rules reproduce all `162/162` default robustness-map verdicts. | [`../docs/analytic-phase-boundary-theorem.md`](../docs/analytic-phase-boundary-theorem.md) | [`../symbolic/analytic_phase_boundary.py`](../symbolic/analytic_phase_boundary.py) | Section 8 |
| The dynamic paper claim is not an instrument forecast or real-data detection claim. | [`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md), [`../docs/failure-ledger.md`](../docs/failure-ledger.md) | Not applicable | Discussion and conclusion |

## Counterexample Candidate Claims

| Claim | Supporting artifact | Script check | Paper placement |
| --- | --- | --- | --- |
| The second-order internal mode is the live A4-violating dynamic counterexample candidate. | [`../docs/second-order-internal-mode.md`](../docs/second-order-internal-mode.md) | [`../symbolic/second_order_mode_response.py`](../symbolic/second_order_mode_response.py) | Section 4 |
| Finite shared acceleration-like or range-like projection can preserve the internal quadratic denominator. | [`../docs/second-order-projection-audit.md`](../docs/second-order-projection-audit.md) | [`../symbolic/second_order_projection_audit.py`](../symbolic/second_order_projection_audit.py) | Section 6 |
| Exact-in-e eccentric forcing can provide enough harmonics to support resonance-aware budget-breaking designs in some parameter regions. | [`../docs/exact-in-e-resonant-forcing-budget-audit.md`](../docs/exact-in-e-resonant-forcing-budget-audit.md), [`../docs/amplitude-weighted-resonant-design-theorem.md`](../docs/amplitude-weighted-resonant-design-theorem.md), [`../docs/physical-detectability-map.md`](../docs/physical-detectability-map.md) | [`../symbolic/exact_in_e_resonant_forcing_audit.py`](../symbolic/exact_in_e_resonant_forcing_audit.py), [`../symbolic/amplitude_weighted_resonant_design.py`](../symbolic/amplitude_weighted_resonant_design.py), [`../symbolic/physical_detectability_map.py`](../symbolic/physical_detectability_map.py) | Section 6 |
| Minimal nonlinear second-order drive/readout can generate lines that carry the same quadratic denominator as the linear response. | [`../docs/nonlinear-second-order-mode.md`](../docs/nonlinear-second-order-mode.md) | [`../symbolic/nonlinear_second_order_mode.py`](../symbolic/nonlinear_second_order_mode.py) | Section 7 |
| Moderate-eccentricity examples can be generated-budget-breaking under relative cutoff or parameterized SNR criteria. | [`../docs/nonlinear-second-order-detectability.md`](../docs/nonlinear-second-order-detectability.md) | [`../symbolic/nonlinear_second_order_detectability.py`](../symbolic/nonlinear_second_order_detectability.py) | Section 7 |
| Moderate-eccentricity designs can be component-separable-and-budget-ready. | [`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md) | [`../symbolic/component_separability_audit.py`](../symbolic/component_separability_audit.py) | Section 8 |
| The nonlinear second-order positive region is not isolated inside the default dimensionless grid. | [`../docs/nonlinear-robustness-map.md`](../docs/nonlinear-robustness-map.md) | [`../symbolic/nonlinear_robustness_map.py`](../symbolic/nonlinear_robustness_map.py) | Section 8 |

## Conjectural Support Or Future Routes

| Statement | Supporting artifact | Boundary |
| --- | --- | --- |
| A system-anchored detectability note could translate design variables into a concrete channel example. | [`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md), [`dynamic-loophole-outline.md`](dynamic-loophole-outline.md) | Future support material, not part of the current theorem claim. |
| The dynamic paper is best separated from the static theorem/classification paper. | [`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md), [`dynamic-loophole-outline.md`](dynamic-loophole-outline.md) | Packaging recommendation, not a scientific theorem. |
| A merged manuscript would need to put the A4 dynamic loophole in the main text. | [`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md) | Editorial route only. |

## Unsupported Or Forbidden In This Package

| Statement | Verdict |
| --- | --- |
| The dynamic branch predicts a specific real instrument detection. | Unsupported; forbidden by current package scope. |
| Arbitrary per-frequency projection nuisance still leaves the signal identifiable. | False under the current artifacts. |
| Sideband existence alone is a dynamic signature. | False under the current comparator audit. |
| The robustness-grid fraction is a physical population probability. | Unsupported; explicitly disclaimed. |
| The one-state branch is the main positive result. | False for the current forcing/projection dictionaries. |
