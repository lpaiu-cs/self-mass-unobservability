# Theorem A Candidate: Delta<=4 Contraction-Level Completeness For The Exact Current Primitive Set

## Statement

- Status: Conjectural. Work in the minimal free-fall sector only: nonspinning, nearly spherical, parity-even, local worldline EFT, fixed truncation order `\Delta \le 4`, and finite external field content.
- Status: Conjectural. Assume the body is quasi-static, self-bound in equilibrium, and carries no orbital-timescale internal state variable in the free-fall sector.
- Status: Conjectural. Assume couplings to the external environment are analytic once the local operator space has been reduced to normal form at fixed order.
- Status: Imported from prior work. Exclude literal internal self-unobservability as an equilibrium mechanism.
- Status: Imported from prior work. Exclude COM self-subtraction as a source of a new monopole force.

- Status: Conjectural. Then, for the exact current primitive set at `Delta_max = 4`, the admissible local free-fall scalar operator space closes on the corrected normal-form list

```math
\{E2,\ E3,\ E2^2,\ dotE2,\ gradE2,\ divE2,\ mixedGradE2\}
```

modulo total derivatives, lower-order equations of motion, and the explicitly stated algebraic identities.
- Status: Conjectural. Consequently, by [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md), the leading body-dependent free-fall deviation collapses to finitely many sensitivity coordinates plus finitely many higher-multipole Wilson coefficients.

## Burden Split

- Status: Proven. This theorem candidate does not assume finite scalar-basis closure as an axiom.
- Status: Conjectural. The conditional EFT-collapse step lives in [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md).
- Status: Conjectural. The current M4 path uses explicit contraction enumeration and explicit reduction rules rather than absence arguments, organized by [`primitive-catalog.md`](primitive-catalog.md), [`reduction-rules.md`](reduction-rules.md), [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md), and [`../lemmas/08-mixed-time-derivative-audit.md`](../lemmas/08-mixed-time-derivative-audit.md).

## Proof Route

1. Status: Proven. Use the fixed-order counting rule in [`power-counting.md`](power-counting.md) to define a bounded-weight operator search space.
2. Status: Proven. Finite external field content and positive weights imply a finite set of decorated primitive building blocks at order `\Delta \le 4`.
3. Status: Proven. [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) exhaustively enumerates all parity-even scalar contractions from the exact current primitive blocks through `Delta<=4`.
4. Status: Proven. [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) reduces all total-derivative, lower-order-EOM, and Cayley-Hamilton-reducible classes explicitly.
5. Status: Conjectural. The surviving classes are exactly the corrected seven-element normal-form list above.
6. Status: Conjectural. Invoke the conditional collapse lemma to convert that finite normal-form basis into finitely many sensitivity coordinates and Wilson coefficients.

## Exact Current Status

- Status: Proven. The previous five-element `Delta<=4` target list was incomplete.
- Status: Proven. The explicit omitted surviving operators are `divE2` and `mixedGradE2`.
- Status: Conjectural. Within the exact current primitive set, no smaller surviving obstruction has been found beyond those corrected gradient invariants.

## Dependencies

- Status: Imported from prior work. [`../lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md)
- Status: Imported from prior work. [`../lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md)
- Status: Conjectural. [`../lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md)
- Status: Proven. [`../lemmas/05-finite-basis-closure.md`](../lemmas/05-finite-basis-closure.md)
- Status: Conjectural. [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md)
- Status: Conjectural. [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md)
- Status: Proven. [`../lemmas/08-mixed-time-derivative-audit.md`](../lemmas/08-mixed-time-derivative-audit.md)
- Status: Conjectural. [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md)
- Status: Conjectural. [`power-counting.md`](power-counting.md)
- Status: Proven. [`primitive-catalog.md`](primitive-catalog.md)
- Status: Proven. [`reduction-rules.md`](reduction-rules.md)

## Failure Triggers

- Status: Counterexample candidate. An explicit `chi_A` state invalidates the no-state reduction step.
- Status: Counterexample candidate. Nonanalytic activation invalidates the analytic collapse step.
- Status: Counterexample candidate. Hereditary couplings invalidate locality before basis closure is even applied.
- Status: Counterexample candidate. Infinite or uncontrolled external field content invalidates the fixed-order closure theorem candidate.
- Status: Counterexample candidate. Any additional `Delta<=4` scalar contraction from the current primitive set that evades the explicit enumeration/reduction pipeline would be a direct obstruction to the current M4 path.
