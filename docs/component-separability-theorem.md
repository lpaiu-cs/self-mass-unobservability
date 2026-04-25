# Component Separability Theorem

Status: Counterexample candidate. This note addresses the M14 bottleneck left
by the nonlinear second-order detectability audit: generated nonlinear
response can land on the same exact-in-e harmonic as the linear response, so
budget-breaking generated lines are not enough unless the generated component
is locally identifiable.

## Joint Harmonic Model

Status: Counterexample candidate. For each usable harmonic `k`, model the
deprojected or finitely projected observable as

```math
\hat O_k
=
\Lambda_k
\left[
Q_Y A_k(p,e)
+Q_L A_k(p,e)H_2(k)
+Q_\beta B_k(p,e)D_2(k)^{-1}
\right],
```

where

```math
B_k(p,e)=A_k(2p,e),
\qquad
H_2(k)=\frac{\rho^2}{D_2(k)},
\qquad
D_2(k)=\rho^2-k^2+i\delta k .
```

Status: Conjectural. The amplitudes `Q_Y`, `Q_L`, and `Q_beta` are shared
fit parameters. The denominator parameters `rho` and `delta` are also shared.
For the range channel, a finite shared `kappa` projection nuisance may be
included.

## Rank Criterion

Status: Proven. A generated component is locally separable only if the
Jacobian column for `Q_beta` adds rank relative to the model without
`Q_beta`.

Status: Proven. A full local separability design additionally requires the
real Jacobian rank to equal the number of shared parameters included in the
joint model.

Status: Proven. A theorem-level nonlinear observable still requires the
generated-sideband budget from
[`nonlinear-second-order-detectability.md`](nonlinear-second-order-detectability.md):

```math
N_{\rm usable}\ge M+2+K_\Lambda .
```

Thus there are three distinct outcomes:

- Status: Proven. `generated-component-degenerate`: the `Q_beta` column does
  not add rank.
- Status: Proven. `component-separable-but-budget-underdesigned`: `Q_beta` is
  locally separable, but generated samples remain below the static nonlinear
  comparator budget.
- Status: Counterexample candidate. `component-separable-and-budget-ready`:
  the joint rank is full, `Q_beta` adds rank, and the generated sample budget
  is exceeded.

## Default Audit

Status: Counterexample candidate. The deterministic audit uses `p=2`,
`rho=3/2`, `delta=0.2`, relative generated-line cutoff `eta=10^-3`, and
harmonics through `H=6`.

Status: Counterexample candidate. At `e=0.1`, acceleration and fixed-kappa
range are locally component-separable but budget-underdesigned: the usable
set is `{1,2,3}`, while the generated comparator budgets require `4` and `5`
samples respectively. With free range `kappa`, the same `e=0.1` design is
generated-component-degenerate.

Status: Counterexample candidate. At `e=0.3`, the acceleration channel is
component-separable-and-budget-ready with usable harmonics `{1,2,3,4,5}` and
minimum cutoff `H=4`. The fixed-kappa and free-kappa range channels are also
component-separable-and-budget-ready, with minimum cutoff `H=5`.

Status: Counterexample candidate. At `e=0.6`, all default acceleration/range
designs through `H=6` are component-separable-and-budget-ready.

## Collapse Boundaries

Status: Proven. If `B_k(p,e)` is proportional to `A_k(p,e)` over the usable
harmonics, the generated amplitude column is degenerate with the linear
amplitude column.

Status: Proven. If the usable generated harmonic count is below
`M+2+K_Lambda`, the design remains a static nonlinear mimicry candidate even
when `Q_beta` is locally separable.

Status: Proven. If projection nuisance is arbitrary per frequency, the joint
rank theorem is not an observable claim; every harmonic can be absorbed by a
frequency-local nuisance.

Status: Counterexample candidate. The current positive branch therefore
requires moderate eccentricity, enough usable harmonics, finite shared
projection nuisance, and full joint rank including the generated amplitude.
