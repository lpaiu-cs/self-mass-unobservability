# Dynamic Loophole Paper Synthesis

Status: Counterexample candidate. This note packages the dynamic-visibility
side track after the analytic phase-boundary theorem. It does not introduce a
new loophole model. Its purpose is to identify the paper-level claim stack,
the no-go/control branch, and the remaining non-claims.

## Thesis

Status: Counterexample candidate. The dynamic branch attacks A4 directly:
the assumption that no orbital-timescale internal state variable appears in
the free-fall sector.

Status: Counterexample candidate. The strongest honest thesis is:

```text
One-state relaxation is a control/no-go branch under the current
forcing/projection dictionaries, while nonlinear second-order internal
visibility supplies a parameterized, non-fine-tuned counterexample candidate:
exact-in-e eccentric forcing plus finite shared projection nuisance can yield
budget-breaking, locally separable generated-line observables with analytic
phase boundaries.
```

Status: Proven. This is not an instrument forecast. The quantities `Q_0`,
`Q_beta`, `sigma_0`, `Gamma`, `kappa`, `H`, `eta`, and `SNR_min` remain
design inputs unless a later system-anchored note fixes them to a concrete
observable channel.

## Result Ledger

| Branch | Status | Main artifact | Result | Caveat |
| --- | --- | --- | --- | --- |
| One-state relaxation | Proven | [`sample-budget-theorem.md`](sample-budget-theorem.md), [`orbital-budget-case.md`](orbital-budget-case.md), [`two-tone-budget-case.md`](two-tone-budget-case.md) | Shared one-pole transfer and nonlinear sidebands are underbudget for the current orbital and two-tone dictionaries. | Useful as control/no-go, not as the final positive claim. |
| Linear second-order mode | Counterexample candidate | [`second-order-internal-mode.md`](second-order-internal-mode.md), [`physical-detectability-map.md`](physical-detectability-map.md) | A quadratic denominator gives resonance-aware budget-breaking regions under parameterized SNR criteria. | Still conditional on finite shared projection nuisance and usable harmonics. |
| Nonlinear second-order mode | Counterexample candidate | [`nonlinear-second-order-mode.md`](nonlinear-second-order-mode.md), [`nonlinear-second-order-detectability.md`](nonlinear-second-order-detectability.md) | Generated lines can carry the same `D_2(z)` denominator as linear lines and can exceed finite static nonlinear comparator budgets. | A single generated line is still static-nonlinear-degenerate. |
| Component separability | Counterexample candidate | [`component-separability-theorem.md`](component-separability-theorem.md) | The generated component is locally separable in moderate-eccentricity designs by Jacobian rank. | Low-eccentricity and free-`kappa` cases can remain degenerate or underdesigned. |
| Robustness map | Counterexample candidate | [`nonlinear-robustness-map.md`](nonlinear-robustness-map.md) | The default dimensionless grid has `107/162` robust-positive points. | The grid fraction is not a physical prior or population measure. |
| Analytic phase boundary | Proven | [`analytic-phase-boundary-theorem.md`](analytic-phase-boundary-theorem.md), [`../symbolic/analytic_phase_boundary.py`](../symbolic/analytic_phase_boundary.py) | Count, bracket, amplitude-usability, and rank boundaries reproduce `162/162` robustness verdicts. | Passing all boundaries is still a parameterized design theorem, not a detection claim. |

## Paper-Level Claim Stack

Status: Imported from prior work. The static theorem program already showed
that minimal-sector uniqueness is weak once magnetic and scalar primitive
families are admitted, while broader fixed-order finite-family collapse
remains the positive theorem target.

Status: Proven. The dynamic side track changes the attacked assumption from
static primitive-family adequacy to A4. More static bookkeeping is not the
mainline of this result.

Status: Proven. The one-state branch establishes the right comparator
discipline: phase lag, a sideband, or a finite set of ratios does not count
unless the finite shared comparator and projection nuisance budgets are
exceeded.

Status: Counterexample candidate. The second-order branch supplies the live
counterexample class because a shared quadratic denominator can survive finite
shared projection and organize both linear and generated-line responses.

Status: Counterexample candidate. The nonlinear second-order branch is the
current strongest positive result: it has budget-breaking generated lines,
local component separability, finite-grid robustness, and analytic
phase-boundary conditions.

## Non-Claims

Status: Proven. The branch does not claim real-data detectability.

Status: Proven. The branch does not claim identifiability under arbitrary
per-frequency projection nuisance.

Status: Proven. The branch does not claim that sideband existence alone is a
dynamic signature once static nonlinear comparators are admitted.

Status: Proven. The robustness-map positive fraction is not a population
measure. It only shows that the positive region is not isolated inside the
chosen dimensionless audit grid.

Status: Conjectural. A future system-anchored detectability note could attach
the design inputs to one concrete free-fall-style channel, but that is support
material rather than the main theorem result.

## Suggested Paper Split

Status: Conjectural. The cleanest packaging is a two-paper split:

| Paper | Scope | Message |
| --- | --- | --- |
| A: theorem/classification | Static sensitivity-collapse theorem, family-admission obstructions, finite-family fixed-order collapse target. | Minimal-sector uniqueness is fragile; the positive theorem target is explicit finite-family collapse under stated assumptions. |
| B: dynamic loophole | One-state no-go/control, linear second-order detectability, nonlinear generated-line detectability, separability, robustness, analytic boundaries. | Violating A4 with a second-order internal mode yields a conditional but non-fine-tuned dynamic counterexample class. |

Status: Conjectural. If the material is merged into one manuscript, the A4
dynamic loophole should be a main section rather than an appendix-only aside,
because it supplies the smallest explicit route by which the theorem
assumptions can fail observationally.

## Dynamic Paper Outline

Status: Conjectural. A dynamic-loophole paper can be organized as follows:

1. Status: Imported from prior work. Motivation: sensitivity collapse and the
   role of A4.
2. Status: Proven. Comparator discipline: finite samples do not defeat finite
   shared comparators unless their budgets are exceeded.
3. Status: Proven. One-state relaxation as control/no-go: shared pole, phase
   lag, sideband mimicry, and sample-budget failure for current dictionaries.
4. Status: Counterexample candidate. Second-order internal mode: transfer,
   phase wrapping, first-order collapse limit, and projection survival.
5. Status: Counterexample candidate. Exact-in-e eccentric forcing:
   resonance-aware sample provisioning and parameterized physical
   detectability.
6. Status: Counterexample candidate. Nonlinear second-order response:
   generated lines carrying the shared `D_2` denominator.
7. Status: Counterexample candidate. Component separability by joint harmonic
   Jacobian rank.
8. Status: Counterexample candidate. Robustness map over dimensionless design
   variables.
9. Status: Proven. Analytic phase-boundary theorem and collapse/no-go regions.
10. Status: Proven. Limits: no instrument forecast, no arbitrary projection
    nuisance, no sideband-only novelty claim.

## Optional Last Support Note

Status: Conjectural. The only high-value optional follow-up before manuscript
work is a short system-anchored detectability note. It would choose one
free-fall-style observable channel and translate `Q_0`, `Q_beta`, `sigma_0`,
`Gamma`, `kappa`, cadence, and harmonic extraction error into a concrete
design example.

Status: Proven. Opening a new internal-state family is lower priority than
packaging the current result, because the branch already has a one-state
control/no-go and a nonlinear second-order conditional positive theorem.
