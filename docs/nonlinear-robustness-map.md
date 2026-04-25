# Nonlinear Robustness Map

Status: Counterexample candidate. This note asks whether the M14 positive
region is an isolated point or a broader parameterized region. It does not add
a new loophole model; it stress-tests the nonlinear second-order branch.

## Grid

Status: Counterexample candidate. The default deterministic grid is

```text
p in {1.5, 2, 3}
e in {0.1, 0.3, 0.6}
rho in {4/3, 3/2, 2}
delta in {0.1, 0.2}
projection in {acceleration, range-fixed-kappa, range-free-kappa}
H = 6
eta = 10^-3
M = 1
```

Status: Conjectural. This is a finite robustness map, not a measure on the
physical parameter space. It tests whether the positive branch persists under
nearby dimensionless design changes.

## Verdict Rule

Status: Proven. A grid point is `robust-positive` only when all of the
following hold:

- generated harmonics exceed the static nonlinear comparator budget;
- usable harmonics bracket `rho`;
- the joint harmonic Jacobian is full rank;
- the generated amplitude column `Q_beta` adds rank;
- projection nuisance is finite and shared.

Status: Proven. If `Q_beta` adds rank but the generated sample count is below
budget, the verdict is
`component-separable-but-budget-underdesigned`, not positive.

Status: Proven. If `Q_beta` does not add rank after nuisance parameters are
included, the verdict is `generated-component-degenerate`.

## Default Result

Status: Counterexample candidate. The default grid has `162` points. The
deterministic audit finds:

| Verdict | Count |
| --- | ---: |
| `robust-positive` | 107 |
| `component-separable-but-budget-underdesigned` | 40 |
| `generated-component-degenerate` | 15 |

Status: Counterexample candidate. The finite-grid robust-positive fraction is
`107/162 = 0.660`. This indicates that the M14 positive region is not a single
fine-tuned point in this dimensionless grid.

Status: Counterexample candidate. By projection channel:

| Projection mode | Positive | Underdesigned | Degenerate |
| --- | ---: | ---: | ---: |
| acceleration | 37 | 17 | 0 |
| range-fixed-kappa | 35 | 19 | 0 |
| range-free-kappa | 35 | 4 | 15 |

Status: Counterexample candidate. The strongest remaining fragility in this
grid is not acceleration/range projection itself, but low-eccentricity sample
underdesign and free-`kappa` degeneracy in part of the low-information range
sector.

## Interpretation

Status: Counterexample candidate. The nonlinear second-order positive branch
is broader than the single `p=2,e=0.3,rho=3/2,delta=0.2` example. It survives
across multiple `p`, `rho`, `delta`, and projection choices in the default
grid.

Status: Proven. This robustness map is still conditional. It collapses if
projection nuisance becomes arbitrary per frequency, if generated components
cannot be separated at all, if the generated sample budget is not exceeded, or
if the harmonic amplitudes fall below the usable threshold.

Status: Counterexample candidate. The next paper-level move is to present the
one-state branch as control/no-go, the linear second-order branch as a
parameterized design theorem, and the nonlinear second-order branch as the
current strongest conditional positive result with an explicit robustness
map.

## Analytic Boundary

Status: Counterexample candidate. The analytic count/bracket/rank boundary
that explains this finite grid is recorded in
[`analytic-phase-boundary-theorem.md`](analytic-phase-boundary-theorem.md).
The associated audit reproduces all `162/162` default grid verdicts.
