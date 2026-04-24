# Conditional Collapse Lemma

## Statement

- Status: Conjectural. Fix a truncation order `\Delta \le \Delta_{\max}` and suppose the minimal free-fall sector admits a finite scalar invariant basis `Y^I` and a finite set of higher-rank local operators `{\cal O}_\alpha` at that order.
- Status: Conjectural. Assume the body is quasi-static, nearly spherical, nonspinning, parity-even, self-bound in equilibrium, and described by a local worldline EFT with no orbital-timescale internal state variable.
- Status: Conjectural. Assume the monopole coupling is analytic in the invariant coordinates `Y^I` near the reference background.

- Status: Conjectural. Then the free-fall worldline action can be written as

```math
S_A
=
\int d\tau \left[
-m_A(Y)
+ \sum_\alpha C_{A,\alpha}\,{\cal O}_\alpha
\right],
```

with

```math
m_A(Y)
=
m_A^{(0)}
+\sum_{n=1}^{N}
\frac{m_A^{(0)}}{n!}
s_{A,I_1\cdots I_n}
Y^{I_1}\cdots Y^{I_n}
+ O(\Delta_{\max}+1),
```

so the leading body-dependent free-fall deviation collapses to finitely many sensitivity coordinates `s_{A,I_1\cdots I_n}` plus finitely many higher-multipole Wilson coefficients `C_{A,\alpha}`.

- Status: Conjectural. The one-dimensional scalar-`s_A` EFT is only a corollary, obtained only if the finite response manifold above is rank one in the observable family under study.

## Why This Lemma Is Conditional

- Status: Proven. This lemma does not prove basis closure.
- Status: Proven. The burden moved out of this lemma is exactly the proof that the admissible operator space at fixed order has a finite normal-form basis.

## Proof Skeleton

1. Status: Conjectural. Write the most general local worldline action in the assumed finite basis `Y^I` and finite operator list `{\cal O}_\alpha`.
2. Status: Conjectural. Use the no-internal-state assumption to rule out additional explicit orbital-timescale coordinates and history functionals in the monopole sector.
3. Status: Conjectural. Expand the analytic monopole function `m_A(Y)` in the finitely many coordinates `Y^I`.
4. Status: Conjectural. Truncate at the chosen fixed order `\Delta_{\max}`. Only finitely many Taylor coefficients survive.
5. Status: Conjectural. Classify the remaining body dependence into higher-rank local operators with Wilson coefficients `C_{A,\alpha}`.

## Dependencies

- Status: Imported from prior work. [`../lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md)
- Status: Imported from prior work. [`../lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md)
- Status: Conjectural. [`../lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md)
- Status: Conjectural. [`power-counting.md`](power-counting.md)
- Status: Conjectural. [`../lemmas/05-finite-basis-closure.md`](../lemmas/05-finite-basis-closure.md)
