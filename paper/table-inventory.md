# Table Inventory

This file defines the minimal table package for the dynamic loophole paper.

## Table 1: Theorem / No-Go / Counterexample Claim Stack

- Scientific role: Put the main claims and statuses in front of the reader.
- Source artifacts:
  [`claim-stack.md`](claim-stack.md),
  [`result-ledger.md`](result-ledger.md),
  [`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md).
- Proposed columns: Claim, status class, supporting artifact, paper section,
  caveat.
- Caption: The dynamic paper separates proven comparator/no-go statements
  from the nonlinear second-order counterexample candidate.
- What this does not claim: It does not promote counterexample candidates into
  detected effects.

## Table 2: Collapse Boundaries And Non-Claims

- Scientific role: Keep the paper's scientific caution explicit.
- Source artifacts:
  [`non-claims.md`](non-claims.md),
  [`../docs/failure-ledger.md`](../docs/failure-ledger.md),
  [`../docs/assumptions-ledger.md`](../docs/assumptions-ledger.md).
- Proposed columns: Boundary, failure mode, affected branch, repository
  artifact.
- Caption: The dynamic loophole collapses under arbitrary projection nuisance,
  insufficient sample budget, insufficient generated amplitude, failed
  component separability, or adiabatic/first-order limits.
- What this does not claim: It does not provide an exhaustive list of every
  possible future dynamic loophole.

## Table 3: Key Parameterized Positive-Region Examples

- Scientific role: Provide concrete dimensionless examples without turning
  them into instrument claims.
- Source artifacts:
  [`../docs/physical-detectability-map.md`](../docs/physical-detectability-map.md),
  [`../docs/nonlinear-second-order-detectability.md`](../docs/nonlinear-second-order-detectability.md),
  [`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md),
  [`../docs/nonlinear-robustness-map.md`](../docs/nonlinear-robustness-map.md).
- Proposed rows:
  - Linear second-order, `p=2`, `e=0.3`, `rho=3/2`, `delta=0.2`, acceleration,
    minimum `Q_0=0.909151`.
  - Linear second-order, same dimensionless design, range, minimum
    `Q_0=0.288909`.
  - Nonlinear generated-line branch, `e=0.3`, acceleration, minimum
    `Q_beta=0.429654`.
  - Nonlinear generated-line branch, `e=0.3`, range, minimum
    `Q_beta=0.118468`.
  - Robustness grid, `107/162` robust-positive, with the explicit note that
    this is not a population measure.
- Caption: Parameterized examples show where the theorem-level design
  conditions are satisfied.
- What this does not claim: The numerical scales are not tied to a real
  instrument unless a future system-anchored note supplies that mapping.
