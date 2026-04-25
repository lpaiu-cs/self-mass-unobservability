# Dynamic Loophole Non-Claims

This file records what the dynamic paper package does not claim. These
boundaries must stay visible in the submission draft and reviewer notes.

## Required Non-Claims

| Non-claim | Repository basis | How to phrase it |
| --- | --- | --- |
| This is not an instrument forecast. | [`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md), [`../docs/physical-detectability-map.md`](../docs/physical-detectability-map.md), [`../docs/failure-ledger.md`](../docs/failure-ledger.md) | The result is a parameterized design theorem; `Q_0`, `Q_beta`, `sigma_0`, `Gamma`, `kappa`, cadence, and extraction errors are not fixed to an instrument. |
| This is not a real-data detectability claim. | [`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md) | The package identifies conditions under which a design can be budget-breaking; it does not analyze a data set. |
| The model is not identifiable under arbitrary per-frequency projection nuisance. | [`../docs/second-order-projection-audit.md`](../docs/second-order-projection-audit.md), [`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md), [`../docs/analytic-phase-boundary-theorem.md`](../docs/analytic-phase-boundary-theorem.md) | Projection must be known, calibrated, or finite-dimensional with shared nuisance parameters. |
| Sideband existence alone is not a dynamic signature. | [`../docs/nonlinear-comparator-audit.md`](../docs/nonlinear-comparator-audit.md), [`../docs/nonlinear-second-order-detectability.md`](../docs/nonlinear-second-order-detectability.md) | Generated lines matter only after static nonlinear comparator budgets are exceeded. |
| A single generated line is not enough. | [`../docs/nonlinear-second-order-detectability.md`](../docs/nonlinear-second-order-detectability.md) | Generated-line novelty requires at least `M+2+K_Lambda` usable generated samples. |
| Component separability alone is not enough. | [`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md), [`../docs/nonlinear-robustness-map.md`](../docs/nonlinear-robustness-map.md) | A locally identifiable generated component is still underdesigned if the generated sample count is inside the comparator budget. |
| Robustness-map fractions are not population measures. | [`../docs/nonlinear-robustness-map.md`](../docs/nonlinear-robustness-map.md), [`../docs/analytic-phase-boundary-theorem.md`](../docs/analytic-phase-boundary-theorem.md) | The finite grid is a stress test for non-isolation, not a physical prior. |
| The one-state branch is not the positive branch. | [`../docs/sample-budget-theorem.md`](../docs/sample-budget-theorem.md), [`../docs/orbital-budget-case.md`](../docs/orbital-budget-case.md), [`../docs/two-tone-budget-case.md`](../docs/two-tone-budget-case.md) | Present one-state relaxation as control/no-go under the current dictionaries. |

## Wording To Avoid

| Avoid | Use instead |
| --- | --- |
| "The effect is detectable." | "The design is budget-breaking under the stated parameterized SNR or relative-cutoff criterion." |
| "Robust-positive cases are common in nature." | "The positive region is not isolated in the chosen dimensionless grid." |
| "Sidebands prove internal dynamics." | "Sidebands only become useful when budget and separability conditions are both met." |
| "Projection cannot hide the signal." | "Finite shared projection nuisance does not hide the signal in the audited channels." |
| "The one-state model works." | "The one-state model is a control/no-go branch for the current dictionaries." |
| "The theorem is refuted observationally." | "A conditional A4-violating counterexample candidate has a parameterized design region." |
