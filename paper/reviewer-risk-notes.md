# Reviewer Risk Notes

This file lists likely reviewer objections and the exact repository answer.
It should be updated before any submission-style prose is polished further.

## Objection 1: This is only a parameterized design theorem, not a forecast.

Answer: Correct. The paper should say this explicitly. The detectability
variables `Q_0`, `Q_beta`, `sigma_0`, `Gamma`, `kappa`, cadence, and harmonic
extraction uncertainty are design inputs. The package does not bind them to a
real instrument or data set.

Repository answer:
[`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md),
[`../docs/physical-detectability-map.md`](../docs/physical-detectability-map.md),
[`non-claims.md`](non-claims.md).

## Objection 2: Static nonlinear comparators can mimic single sidebands.

Answer: Correct. That is why the paper does not claim novelty from sideband
existence. Generated lines matter only when they exceed the static nonlinear
comparator budget and carry the shared denominator relation.

Repository answer:
[`../docs/nonlinear-comparator-audit.md`](../docs/nonlinear-comparator-audit.md),
[`../docs/nonlinear-second-order-detectability.md`](../docs/nonlinear-second-order-detectability.md),
[`claim-stack.md`](claim-stack.md).

## Objection 3: Generated-line novelty only matters once budget and separability are both satisfied.

Answer: Correct. The generated branch requires both `|U| >= M+2+K_Lambda`
and local component separability by Jacobian rank. The robustness map uses
separate labels for budget-underdesigned and component-degenerate designs.

Repository answer:
[`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md),
[`../docs/nonlinear-robustness-map.md`](../docs/nonlinear-robustness-map.md),
[`../docs/analytic-phase-boundary-theorem.md`](../docs/analytic-phase-boundary-theorem.md).

## Objection 4: The one-state branch failed.

Answer: The one-state branch is intentionally presented as a control/no-go
branch under the current forcing and projection dictionaries. Its role is to
establish the budget discipline and prevent overclaiming from phase lag or
single sidebands.

Repository answer:
[`../docs/sample-budget-theorem.md`](../docs/sample-budget-theorem.md),
[`../docs/orbital-budget-case.md`](../docs/orbital-budget-case.md),
[`../docs/two-tone-budget-case.md`](../docs/two-tone-budget-case.md),
[`result-ledger.md`](result-ledger.md).

## Objection 5: Why is this distinct from the static theorem/classification paper?

Answer: The static paper asks which instantaneous finite-family sensitivity
collapse statements survive under explicit primitive-family assumptions. The
dynamic paper attacks A4 directly by allowing an orbital-timescale internal
state and tracking whether its transfer structure survives finite comparator
and projection budgets.

Repository answer:
[`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md),
[`dynamic-loophole-submission-draft.md`](dynamic-loophole-submission-draft.md).

## Objection 6: The robustness grid could be arbitrary.

Answer: The grid is a stress test, not a physical prior. Its role is to show
that the positive branch is not isolated in the chosen dimensionless design
space. The analytic phase-boundary theorem then explains the grid by count,
bracket, amplitude, and rank boundaries.

Repository answer:
[`../docs/nonlinear-robustness-map.md`](../docs/nonlinear-robustness-map.md),
[`../docs/analytic-phase-boundary-theorem.md`](../docs/analytic-phase-boundary-theorem.md).

## Objection 7: Projection nuisance can absorb the signal.

Answer: Yes, if projection is arbitrary per frequency. The paper only claims
survival under known, calibrated, or finite shared projection nuisance.

Repository answer:
[`../docs/second-order-projection-audit.md`](../docs/second-order-projection-audit.md),
[`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md),
[`non-claims.md`](non-claims.md).
