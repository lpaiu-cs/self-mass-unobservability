# Figure Inventory

This file defines the minimal four-figure package for the dynamic loophole
paper. No figure should imply a stronger claim than the artifacts support.

## Figure 1: One-State Control/No-Go Flow

- Scientific role: Show why the first-order relaxation branch is a control
  and no-go branch under the current dictionaries.
- Source artifacts:
  [`../docs/sample-budget-theorem.md`](../docs/sample-budget-theorem.md),
  [`../docs/orbital-budget-case.md`](../docs/orbital-budget-case.md),
  [`../docs/two-tone-budget-case.md`](../docs/two-tone-budget-case.md).
- Caption: First-order relaxation produces a shared pole and phase lag, but
  the current orbital and two-tone forcing dictionaries remain inside finite
  static comparator budgets.
- What this does not claim: It does not claim the one-state branch is the main
  positive dynamic signature.

## Figure 2: Second-Order Detectability Region

- Scientific role: Show how resonance-aware counting, bracketing, and
  amplitude/SNR usability create a conditional positive design region for the
  linear second-order branch.
- Source artifacts:
  [`../docs/second-order-internal-mode.md`](../docs/second-order-internal-mode.md),
  [`../docs/exact-in-e-resonant-forcing-budget-audit.md`](../docs/exact-in-e-resonant-forcing-budget-audit.md),
  [`../docs/amplitude-weighted-resonant-design-theorem.md`](../docs/amplitude-weighted-resonant-design-theorem.md),
  [`../docs/physical-detectability-map.md`](../docs/physical-detectability-map.md).
- Caption: Exact-in-e eccentric forcing can supply enough usable harmonics to
  bracket the second-order denominator and exceed finite shared comparator
  budgets in parameterized moderate-eccentricity designs.
- What this does not claim: It does not attach `Q_0`, `sigma_0`, `Gamma`, or
  `kappa` to a real instrument.

## Figure 3: Nonlinear Robustness Map

- Scientific role: Show that the nonlinear second-order positive region is
  not isolated inside the chosen dimensionless audit grid.
- Source artifacts:
  [`../docs/nonlinear-second-order-detectability.md`](../docs/nonlinear-second-order-detectability.md),
  [`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md),
  [`../docs/nonlinear-robustness-map.md`](../docs/nonlinear-robustness-map.md),
  [`../symbolic/nonlinear_robustness_map.py`](../symbolic/nonlinear_robustness_map.py).
- Caption: In the default finite grid, robust-positive, underdesigned, and
  generated-component-degenerate regions separate by generated sample budget,
  component rank, and finite projection nuisance.
- What this does not claim: The `107/162` robust-positive fraction is not a
  physical population probability.

## Figure 4: Analytic Phase-Boundary Schematic

- Scientific role: Replace the finite robustness grid with explicit boundary
  logic.
- Source artifacts:
  [`../docs/analytic-phase-boundary-theorem.md`](../docs/analytic-phase-boundary-theorem.md),
  [`../symbolic/analytic_phase_boundary.py`](../symbolic/analytic_phase_boundary.py).
- Caption: The nonlinear second-order verdict is governed by count,
  bracketing, amplitude-usability, and rank boundaries; the audit reproduces
  all default robustness-map classifications.
- What this does not claim: Passing the boundaries is still a parameterized
  design result, not a real-data detection.
