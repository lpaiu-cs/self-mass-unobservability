# Analytic Phase-Boundary Theorem

Status: Counterexample candidate. This note turns the finite nonlinear
robustness grid into explicit phase-boundary conditions. It does not introduce
a new model.

## Count And Bracket Boundary

Status: Proven. For a degree-`M` static nonlinear generated-line comparator
and `K_Lambda` finite shared projection nuisance parameters, the generated
sample-count boundary is

```math
H_{\rm count}=M+2+K_\Lambda .
```

Status: Proven. A single positive-harmonic eccentric dictionary can bracket
the internal phase-wrap location only if `rho>1`. The bracketing cutoff is

```math
H_{\rm bracket}(\rho)
=
\begin{cases}
\lfloor\rho\rfloor+1, & \rho>1,\ \rho\notin\mathbb Z,\\
\rho+1, & \rho\in\mathbb Z,\ \rho>1,\\
\hbox{none}, & \rho\le 1 .
\end{cases}
```

Status: Proven. The formal harmonic minimum is

```math
H_{\min}
=
\max\{H_{\rm count},H_{\rm bracket}(\rho)\},
```

when the bracket exists. If `rho<=1`, the single-system positive-harmonic
design has no bracketing phase-boundary solution.

## Usability Boundary

Status: Counterexample candidate. A formal harmonic cutoff does not count
unless the usable generated set

```math
U_\eta
=
\{k\le H:\ W_k/W_{\max}\ge\eta\}
```

satisfies both

```math
|U_\eta|\ge H_{\rm count},
\qquad
\exists k_-,k_+\in U_\eta:\ k_-<\rho<k_+ .
```

Status: Proven. If the first inequality fails, the phase is
`usable-count-underbudget`. If the second fails, the phase is
`usable-bracket-failure`.

## Component-Rank Boundary

Status: Counterexample candidate. In the fixed-projection amplitude submodel,
the three-column generated-component test uses

```math
\left[
A_k,\quad
A_k\rho^2/D_k,\quad
B_k/D_k
\right].
```

Status: Proven. For three harmonics, multiplying the determinant by
`D_1D_2D_3` gives the rank minor

```math
\Delta_\beta
=
\rho^2[
A_1A_2B_3(D_1-D_2)
-A_1A_3B_2(D_1-D_3)
+A_2A_3B_1(D_2-D_3)
].
```

Status: Proven. If `B_k=c A_k` over the sampled harmonics, then
`\Delta_\beta=0`; the generated component is absorbed into the existing
internal-response amplitude column. More generally, `\Delta_beta != 0` is a
simple sufficient condition that `Q_beta` is not just a static rescaling of
the linear harmonic content in the fixed-nuisance amplitude submodel.

Status: Counterexample candidate. With shared `rho`, `delta`, and optional
range `kappa` nuisance included, the full condition is the real Jacobian-rank
condition recorded in [`component-separability-theorem.md`](component-separability-theorem.md):
`Q_beta` must add rank and the full Jacobian must have rank equal to the
number of shared parameters.

## Phase Labels

Status: Proven. The phase-boundary audit uses the following ordered
classification:

1. `no-single-system-bracket`: `rho<=1`.
2. `harmonic-cutoff-under-minimum`: `H < H_min`.
3. `rank-degenerate`: `Q_beta` does not add rank after shared nuisance.
4. `nuisance-rank-degenerate`: `Q_beta` adds rank, but the full shared model
   is not full rank.
5. `usable-count-underbudget`: usable generated samples are below
   `H_count`.
6. `usable-bracket-failure`: usable generated samples do not bracket `rho`.
7. `analytic-positive`: count, bracket, amplitude, and rank conditions all
   hold.

Status: Proven. This boundary reproduces the finite robustness-map verdicts:
the deterministic audit maps all `162/162` default grid points to the same
positive/underdesigned/degenerate classification used by
[`nonlinear-robustness-map.md`](nonlinear-robustness-map.md).

## Interpretation

Status: Counterexample candidate. The nonlinear second-order positive branch
is controlled by four explicit boundaries:

- generated count: `M+2+K_Lambda`;
- bracketing: `rho>1` and usable harmonics on both sides of `rho`;
- amplitude: `W_k/W_max >= eta` or an SNR analogue;
- rank: `Q_beta` adds full shared-parameter Jacobian rank.

Status: Proven. Failure of any one boundary is enough to demote the design to
a no-go/collapse region. Passing all four is the current strongest
parameterized positive theorem in this worktree.
