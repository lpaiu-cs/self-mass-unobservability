# Sample-Budget Theorem

Status: Counterexample candidate. This note turns the shared-tau sideband
program into an explicit sample-budget theorem/no-go statement.

## Observational Question

Status: Counterexample candidate. The data are deprojected and
forcing-normalized complex samples:

```math
L_i=L(z_i),
\qquad
S_j=S(u_j,v_j).
```

Here `L_i` are linear-frequency samples and `S_j` are generated sideband-pair
samples. Let

```math
N_L = \#\{z_i\},
\qquad
N_S = \#\{(u_j,v_j)\}.
```

Status: Proven. A finite-data novelty claim is not allowed if a finite shared
static comparator plus projection nuisance can interpolate the same samples.

## Comparator Budgets

Status: Conjectural. The conservative complex static comparator is

```math
L_{\rm stat}(z)=P_N(z),
\qquad
S_{\rm stat}(u,v)=Q_M(u,v),
```

where `P_N` is a degree-`N` univariate polynomial and `Q_M` is a bivariate
polynomial of total degree `M`.

Status: Proven. The shared complex coefficient budgets are

```math
D_L=N+1,
\qquad
D_S=\frac{(M+1)(M+2)}{2},
\qquad
D_{\rm tot}=D_L+D_S+K_\Lambda ,
```

where `K_Lambda` is the number of shared projection nuisance parameters.

## Budget-Breaking Criterion

Status: Proven. A design does not beat the linear shared polynomial budget
unless

```math
N_L > D_L=N+1 .
```

Status: Proven. A design does not beat the static nonlinear sideband budget
unless

```math
N_S > D_S=\frac{(M+1)(M+2)}{2}.
```

Status: Proven. A design does not beat the combined shared nuisance budget
unless

```math
N_L+N_S > D_{\rm tot}=N+1+\frac{(M+1)(M+2)}{2}+K_\Lambda .
```

Status: Counterexample candidate. The shared-tau observable is budget-breaking
only when all three inequalities hold.

## Minimum Design

Status: Proven. One canonical minimum design is

```math
N_L^{\min}=N+2,
```

and

```math
N_S^{\min}
=
\max\left[
\frac{(M+1)(M+2)}{2}+1,
\frac{(M+1)(M+2)}{2}+K_\Lambda
\right].
```

Status: Proven. The corresponding total minimum is

```math
N_{\rm tot}^{\min}
=
N_L^{\min}+N_S^{\min}.
```

This allocation puts the projection-nuisance overhead into the sideband-pair
count. Other allocations are equivalent if they satisfy the three inequalities
above.

## Theorem Statement

Status: Counterexample candidate. For fixed finite `N`, `M`, and finite
`K_Lambda`, a shared-tau sideband program can produce a distinguishability
claim only from sample designs satisfying

```math
N_L\ge N+2,
\qquad
N_S\ge \frac{(M+1)(M+2)}{2}+1,
\qquad
N_L+N_S\ge N+2+\frac{(M+1)(M+2)}{2}+K_\Lambda .
```

Status: Proven. If any inequality fails, the result is an underbudget no-go
against the stated comparator class, not a novelty claim.

Status: Proven. If the projection nuisance is arbitrary per frequency or per
sideband pair, then `K_Lambda` grows with the sample count and no finite
sample budget can establish the shared-tau law in that channel.

