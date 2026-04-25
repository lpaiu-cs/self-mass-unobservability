# Dynamic Paper Result Ledger

This ledger gives the paper-facing order of results. It is a packaging index,
not a new theorem note.

| Step | Result | Status in paper | Primary artifacts | Collapse or caveat |
| --- | --- | --- | --- | --- |
| R1 | A4 is the active assumption attacked by the dynamic paper. | Counterexample candidate framing | [`../docs/dynamic-loophole-paper-synthesis.md`](../docs/dynamic-loophole-paper-synthesis.md), [`../docs/assumptions-ledger.md`](../docs/assumptions-ledger.md) | The paper does not reopen static primitive-family audits. |
| R2 | One-state relaxation creates phase lag and a shared pole. | Control setup | [`../counterexamples/chi-state/README.md`](../counterexamples/chi-state/README.md), [`../docs/shared-tau-ratio-theorem.md`](../docs/shared-tau-ratio-theorem.md) | One-frequency quadrature can be derivative-EFT-degenerate. |
| R3 | Static nonlinear comparators can mimic finite sideband data inside budget. | Proven no-go discipline | [`../docs/nonlinear-comparator-audit.md`](../docs/nonlinear-comparator-audit.md), [`../docs/sample-budget-theorem.md`](../docs/sample-budget-theorem.md) | Sideband existence alone is not a signature. |
| R4 | Current one-state orbital and two-tone dictionaries are underbudget. | Proven control/no-go | [`../docs/orbital-budget-case.md`](../docs/orbital-budget-case.md), [`../docs/two-tone-budget-case.md`](../docs/two-tone-budget-case.md) | Do not present one-state as positive. |
| R5 | Second-order internal visibility supplies a shared quadratic denominator. | Counterexample candidate | [`../docs/second-order-internal-mode.md`](../docs/second-order-internal-mode.md) | Collapses in first-order or adiabatic limits. |
| R6 | Finite shared projection can preserve the denominator in acceleration/range channels. | Counterexample candidate | [`../docs/second-order-projection-audit.md`](../docs/second-order-projection-audit.md) | Arbitrary per-frequency projection collapses identifiability. |
| R7 | Exact-in-e forcing provides a resonance-aware harmonic sample dictionary. | Counterexample candidate | [`../docs/exact-in-e-resonant-forcing-budget-audit.md`](../docs/exact-in-e-resonant-forcing-budget-audit.md) | Single-system bracketing fails for `rho<=1`. |
| R8 | Parameterized relative/SNR usability can produce budget-breaking second-order designs. | Counterexample candidate | [`../docs/amplitude-weighted-resonant-design-theorem.md`](../docs/amplitude-weighted-resonant-design-theorem.md), [`../docs/physical-detectability-map.md`](../docs/physical-detectability-map.md) | Not an instrument forecast. |
| R9 | Minimal nonlinear second-order generated lines can carry the same `D_2` denominator. | Counterexample candidate | [`../docs/nonlinear-second-order-mode.md`](../docs/nonlinear-second-order-mode.md) | A single generated line is still static-nonlinear-degenerate. |
| R10 | Generated-line detectability requires usable generated samples beyond `M+2+K_Lambda`. | Proven budget rule and candidate positive examples | [`../docs/nonlinear-second-order-detectability.md`](../docs/nonlinear-second-order-detectability.md) | Low-eccentricity cases can remain generated-underbudget. |
| R11 | Generated contributions can be locally separable from overlapping linear harmonics. | Counterexample candidate | [`../docs/component-separability-theorem.md`](../docs/component-separability-theorem.md) | Full rank and budget must both hold. |
| R12 | The nonlinear second-order positive region persists over the default finite grid. | Counterexample candidate | [`../docs/nonlinear-robustness-map.md`](../docs/nonlinear-robustness-map.md) | Grid fraction is not a population measure. |
| R13 | Count, bracket, amplitude-usability, and rank boundaries explain the grid. | Proven phase-boundary audit | [`../docs/analytic-phase-boundary-theorem.md`](../docs/analytic-phase-boundary-theorem.md), [`../symbolic/analytic_phase_boundary.py`](../symbolic/analytic_phase_boundary.py) | Passing all boundaries is still parameterized, not empirical. |

## Paper-Ready Thesis

The repository is ready for human prose polishing because the package has a
closed control/no-go branch and a conditional positive branch. It should not
be expanded by adding new internal-state models before the present result is
written.
