# Shared-Tau Ratio Theorem

Status: Counterexample candidate. This note promotes the M6 ratio target into
the next theorem/no-go test for the dynamic-visibility loophole. The main
result is now a sample-budget statement, not a sideband-existence statement.

## Normalized Observables

Status: Counterexample candidate. Work with deprojected and forcing-normalized
line amplitudes,

```math
L(\omega)=\frac{\hat O(\omega)}{\hat F(\omega)},
\qquad
S(u,v)=\frac{\hat O(u+v)}{\hat F(u)\hat F(v)} .
```

Status: Counterexample candidate. The one-state nonlinear-drive model predicts

```math
L_{\rm dyn}(z)=G(z)=c_Y+\frac{\beta}{1+\tau_\chi z},
\qquad
S_{\rm dyn}(u,v)=\frac{C_{\rm side}}{1+\tau_\chi(u+v)} .
```

Status: Counterexample candidate. The shared-tau ratio is therefore

```math
R_{\rm dyn}(u,v)
=
\frac{S_{\rm dyn}(u,v)}{L_{\rm dyn}(u)L_{\rm dyn}(v)}
=
\frac{C_{\rm side}}
{[1+\tau_\chi(u+v)]G(u)G(v)} .
```

## Exact Dynamic Diagnostic

Status: Proven. The dynamic model obeys the algebraic diagnostic

```math
R_{\rm dyn}(u,v)G(u)G(v)[1+\tau_\chi(u+v)]-C_{\rm side}=0 .
```

Status: Proven. This identity is stronger than sideband existence because it
requires the same `tau_chi` in the linear transfer factors `G(u),G(v)` and in
the generated sideband factor `1+\tau_chi(u+v)`.

## Static Nonlinear Comparator

Status: Conjectural. The finite shared static nonlinear comparator is modeled
as

```math
L_{\rm stat}(z)=P_N(z),
\qquad
S_{\rm stat}(u,v)=Q_M(u,v),
\qquad
R_{\rm stat}(u,v)=\frac{Q_M(u,v)}{P_N(u)P_N(v)} ,
```

where `P_N` is a finite univariate polynomial and `Q_M` is a finite bivariate
polynomial of total degree `M`.

Status: Proven. The complex interpolation budget is finite:

```math
\#L\hbox{-samples}\le N+1,
\qquad
\#S\hbox{-pairs}\le \frac{(M+1)(M+2)}{2}.
```

Status: Proven. Within that budget, a static nonlinear comparator can absorb a
finite set of measured line and sideband amplitudes without proving a dynamic
state.

Status: Proven. The explicit budget-breaking criterion is recorded in
[`sample-budget-theorem.md`](sample-budget-theorem.md). A shared-tau ratio is
not a novelty claim unless its linear samples, sideband-pair samples, and
projection nuisance budget satisfy that theorem's inequalities.

Status: Counterexample candidate. Beyond that budget, exact absorption of the
shared relaxation law requires a finite polynomial representation of rational
pole factors. That is impossible unless a collapse condition removes the pole
or the nuisance model is enlarged beyond finite shared coefficients.

## Collapse Boundary

Status: Proven. The shared-tau ratio test collapses when `C_side=0`, because no
generated line remains.

Status: Proven. The sideband pole collapses when `tau_chi=0`; then

```math
R_{\rm dyn}(u,v)=\frac{C_{\rm side}}{(c_Y+\beta)^2},
```

a constant ratio that a static nonlinear coefficient can absorb.

Status: Proven. The test also collapses when the comparator is allowed
arbitrary independent complex nuisance values at each linear frequency and
each generated sideband pair.

Status: Counterexample candidate. The test survives only when both the
projection and the static nonlinear comparator are finite shared-nuisance
models across more samples than their interpolation budgets.

## Verdict

Status: Counterexample candidate. The M7 theorem target is not that nonlinear
`chi` creates a sideband. The target is that one common `tau_chi` organizes
`L(u)`, `L(v)`, and `S(u,v)` through a shared pole relation that finite static
local polynomials cannot reproduce exactly once enough distinct samples are
available.

Status: Proven. The current orbital `n,2n,3n` dictionary and the current
two-tone dictionary are classified in [`orbital-budget-case.md`](orbital-budget-case.md)
and [`two-tone-budget-case.md`](two-tone-budget-case.md); neither beats a
realistic `N=1` linear comparator by itself.
