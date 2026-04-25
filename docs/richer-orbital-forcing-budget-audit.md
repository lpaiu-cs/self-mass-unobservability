# Richer Orbital Forcing Budget Audit

Status: Counterexample candidate. This note is the bounded rescue audit for the
one-state `chi` program after the current `n,2n -> 3n` dictionary failed the
sample-budget theorem.

## Counting Model

Status: Conjectural. For generic eccentric forcing `F(Y(t))=(a/r(t))^p`, assume
that the positive harmonics `n,2n,...,Hn` are usable as linear samples after
calibration and that no accidental coefficient zero removes them.

Status: Conjectural. Count clean nonlinear sum-sideband pairs `(i,j)->i+j`
only when

```math
1\le i\le j\le H,
\qquad
H < i+j \le R .
```

Here `H` is the highest retained linear harmonic and `R` is the highest
generated output harmonic kept in the nonlinear sideband audit.

Status: Proven. This clean-count rule reproduces the current small-e case:
`H=2`, `R=3` gives exactly one sideband pair, `(1,2)->3`.

## Minimum Designs

| Comparator class | Minimum `(H,R)` | Clean sideband pairs | Status | Verdict |
| --- | ---: | ---: | --- | --- |
| `N=0, M=0, K_Lambda=0` | `(2,4)` | `2` | Proven | budget-breaking only for the trivial comparator |
| `N=1, M=0, K_Lambda=0` | `(3,4)` | `2` | Proven | budget-breaking for constant sideband comparator |
| `N=1, M=1, K_Lambda=0` | `(3,6)` | `4` | Proven | budget-breaking minimum for first-derivative linear plus linear sideband comparator |
| `N=1, M=1, K_Lambda=2` | `(4,6)` | `4` | Proven | budget-breaking by adding one extra linear sample to pay the projection budget |
| `N=1, M=2, K_Lambda=0` | `(5,8)` | `7` | Proven | budget-breaking for quadratic bivariate sideband comparator |

## Exact-In-E Boundary

Status: Proven. Exact-in-e forcing does not automatically strengthen the clean
new-line argument. If every eccentric harmonic is admitted as a linear input
line, then no generated output frequency is absent from the linear dictionary.

Status: Counterexample candidate. Exact-in-e can still help only if the
observable program can separate the nonlinear sideband contribution from the
linear contribution at the same harmonic. That is a stronger decomposition
assumption than the clean new-line sideband test.

## Verdict

Status: Proven. The current `H=2,R=3` orbital dictionary remains underbudget.

Status: Counterexample candidate. A bounded orbital rescue is possible in the
counting model, but it requires substantially richer harmonic support. The
realistic `N=1,M=1,K_Lambda=0` comparator needs at least `H=3` linear
harmonics and clean sidebands through `R=6`.

Status: Counterexample candidate. With two shared projection nuisances, the
minimum rises to `H=4,R=6`; the fourth linear harmonic pays the combined
budget overhead even though the clean sideband-pair count remains `4`.

Status: Counterexample candidate. Therefore the one-state orbital program is
not dead in principle, but the current small-e dictionary is not close to a
budget-breaking design.
