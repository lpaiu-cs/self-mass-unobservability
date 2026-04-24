# Theorem A Draft: Free-Fall Sensitivity Collapse

## Statement

Status: Conjectural. Let body `A` be a quasi-static, nearly spherical, nonspinning compact body with self-bound equilibrium. Assume:

1. Status: Conjectural. A local worldline EFT exists for the center-of-mass worldline `z_A^mu(tau)` after integrating out short-distance structure.
2. Status: Conjectural. The free-fall sector contains no orbital-timescale internal state variable beyond the worldline degrees of freedom retained explicitly in the EFT.
3. Status: Conjectural. The relevant couplings to external fields are analytic in a neighborhood of a reference background, expressed in a finite basis of scalarized local invariants `Y^I` at the fixed truncation order of interest.
4. Status: Imported from prior work. Literal internal self-unobservability is not used to alter the equilibrium equations, because that branch fails as a compact-body model.
5. Status: Imported from prior work. COM self-subtraction does not generate a new monopole force.

Status: Conjectural. Then the leading body-dependent deviation in the free-fall sector can be written as a finite-jet monopole term plus higher-multipole Wilson coefficients:

```math
S_A
=
\int d\tau \left[
-m_A^{(0)}
-m_A^{(0)}\sum_{n=1}^{N}\frac{1}{n!} s_{A,I_1 \cdots I_n} Y^{I_1}\cdots Y^{I_n}
+\sum_{\ell \ge 2} C_A^{L_\ell}\,{\cal O}_{L_\ell}
\right]
+ O(\epsilon^{N+1}),
```

- Status: Conjectural. The coefficients `s_{A,I_1 ... I_n}` are the monopole sensitivity coordinates of the truncated response manifold.
- Status: Conjectural. The coefficients `C_A^{L_ell}` are higher-multipole Wilson coefficients rather than extra scalar sensitivities.
- Status: Conjectural. No extra leading body-dependent free-fall observable remains unless at least one theorem assumption fails.

Status: Conjectural. The one-dimensional scalar-s EFT is only a corollary: it follows only if the relevant response manifold is rank one in the observable family under study.

## Proof Skeleton

1. Status: Conjectural. Choose a finite invariant basis `Y^I` at the desired truncation order from local external fields and their derivatives, reduced by symmetry and equations of motion.
2. Status: Conjectural. Write the most general local worldline scalar built from `u^mu`, the chosen invariants, and higher-multipole operators. This yields a generic action with a monopole function `m_A(Y)` plus multipole couplings.
3. Status: Imported from prior work. Use the internal-structure no-go lemma to exclude the branch where "self-unobservability" is implemented by directly modifying the equilibrium equations.
4. Status: Imported from prior work. Use the COM decoupling lemma to exclude spurious self-force terms that are pure internal bookkeeping artifacts.
5. Status: Conjectural. Because there is no orbital-timescale internal state variable, the monopole part of the action is algebraic in the instantaneous invariants `Y^I` rather than functional on a history or on an extra worldline field.
6. Status: Conjectural. Analyticity then supplies a finite Taylor jet for `m_A(Y)` at the chosen truncation order. The Taylor coefficients are precisely the sensitivity coordinates `s_{A,I}`, `s_{A,IJ}`, and higher jets.
7. Status: Conjectural. All remaining body dependence that is not in the monopole Taylor jet appears in higher-rank operators, so it is classified as higher-multipole Wilson data rather than as a new scalar sensitivity principle.

## Exact Unresolved Step

- Status: Conjectural. Step 6 is not yet a full proof because the repository still lacks a clean closure argument that the chosen invariant basis is complete and finite for the targeted observable order without hidden redundancies.
- Status: Conjectural. The smallest missing assumption is a precise fixed-order operator-counting rule; adding more than that would risk silently strengthening the theorem.

## Immediate Lemma Dependencies

- Status: Imported from prior work. [`lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md)
- Status: Imported from prior work. [`lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md)
- Status: Conjectural. [`lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md)

## Failure Triggers

- Status: Counterexample candidate. An explicit `chi_A` state invalidates Step 5.
- Status: Counterexample candidate. Nonanalytic activation invalidates Step 6.
- Status: Counterexample candidate. Hereditary couplings invalidate Step 5 by replacing a local function with a memory functional.
