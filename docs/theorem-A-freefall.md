# Theorem A Candidate: Fixed-Order Basis Closure In The Minimal Free-Fall Sector

## Statement

- Status: Conjectural. Work in the M2 sector only: nonspinning, nearly spherical, parity-even, local worldline EFT, fixed truncation order `\Delta \le \Delta_{\max}`, and finite external field content.
- Status: Conjectural. Assume the body is quasi-static, self-bound in equilibrium, and carries no orbital-timescale internal state variable in the free-fall sector.
- Status: Conjectural. Assume couplings to the external environment are analytic once the local operator space has been reduced to normal form at fixed order.
- Status: Imported from prior work. Exclude literal internal self-unobservability as an equilibrium mechanism.
- Status: Imported from prior work. Exclude COM self-subtraction as a source of a new monopole force.

- Status: Conjectural. Then the admissible local free-fall operator space at order `\Delta \le \Delta_{\max}` has a finite normal-form basis modulo total derivatives and lower-order equations of motion.
- Status: Conjectural. Consequently, by [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md), the leading body-dependent free-fall deviation collapses to finitely many sensitivity coordinates plus finitely many higher-multipole Wilson coefficients.

## Burden Split

- Status: Proven. This theorem candidate no longer assumes a finite invariant basis as an axiom.
- Status: Conjectural. The conditional EFT-collapse step lives in [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md).
- Status: Conjectural. The current theorem burden is `normal-form completeness modulo total derivatives and lower-order equations of motion`, organized by [`primitive-catalog.md`](primitive-catalog.md) and [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md).

## Proof Route

1. Status: Proven. Use the fixed-order counting rule in [`power-counting.md`](power-counting.md) to define a bounded-weight operator search space.
2. Status: Proven. Finite external field content and positive weights imply a finite set of decorated primitive building blocks at order `\Delta \le \Delta_{\max}`.
3. Status: Proven. Once a candidate primitive catalog is fixed, the set of admissible parity-even monomials and tensor contractions at that order is finite.
4. Status: Conjectural. Show `normal-form completeness modulo total derivatives and lower-order equations of motion` for the exact `Delta<=4` primitive catalog.
5. Status: Conjectural. Invoke the conditional collapse lemma to convert the finite normal-form basis into finitely many sensitivity coordinates and Wilson coefficients.

## Exact Remaining Unresolved Step

- Status: Conjectural. The unresolved step is no longer "is the candidate monomial set finite?".
- Status: Conjectural. The exact remaining burden is `normal-form completeness modulo total derivatives and lower-order equations of motion` for the exact `Delta<=4` primitive catalog.
- Status: Counterexample candidate. If this reduction fails, the smallest loophole should be recorded as an explicit obstruction operator or an explicit violated assumption rather than as a vague "infinite basis" slogan.

## Dependencies

- Status: Imported from prior work. [`../lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md)
- Status: Imported from prior work. [`../lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md)
- Status: Conjectural. [`../lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md)
- Status: Proven. [`../lemmas/05-finite-basis-closure.md`](../lemmas/05-finite-basis-closure.md)
- Status: Conjectural. [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md)
- Status: Conjectural. [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md)
- Status: Conjectural. [`power-counting.md`](power-counting.md)
- Status: Proven. [`primitive-catalog.md`](primitive-catalog.md)

## Failure Triggers

- Status: Counterexample candidate. An explicit `chi_A` state invalidates the no-state reduction step.
- Status: Counterexample candidate. Nonanalytic activation invalidates the analytic collapse step.
- Status: Counterexample candidate. Hereditary couplings invalidate locality before basis closure is even applied.
- Status: Counterexample candidate. Infinite or uncontrolled external field content invalidates the fixed-order closure theorem candidate.
