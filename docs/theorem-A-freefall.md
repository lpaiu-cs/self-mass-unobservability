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

## Scope Separation

- Status: Proven. The current theorem candidate is only for the exact current primitive set fixed in [`primitive-catalog.md`](primitive-catalog.md).
- Status: Proven. It is not yet a theorem for a physically justified minimal free-fall sector.
- Status: Proven. [`primitive-set-adequacy.md`](primitive-set-adequacy.md) records the first explicit adequacy attack and shows that adding the magnetic tidal family `B_ij` produces a new survivor `B2`.

## Burden Split

- Status: Proven. This theorem candidate does not assume finite scalar-basis closure as an axiom.
- Status: Conjectural. The conditional EFT-collapse step lives in [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md).
- Status: Conjectural. The current M4 path uses explicit contraction enumeration and explicit reduction rules rather than absence arguments, organized by [`primitive-catalog.md`](primitive-catalog.md), [`reduction-rules.md`](reduction-rules.md), [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md), and [`../lemmas/08-mixed-time-derivative-audit.md`](../lemmas/08-mixed-time-derivative-audit.md).
- Status: Proven. The current M5 path adds [`../lemmas/09-survivor-independence-delta4.md`](../lemmas/09-survivor-independence-delta4.md) and separates survivor independence from primitive-set adequacy.

## Proof Route

1. Status: Proven. Use the fixed-order counting rule in [`power-counting.md`](power-counting.md) to define a bounded-weight operator search space.
2. Status: Proven. Finite external field content and positive weights imply a finite set of decorated primitive building blocks at order `\Delta \le 4`.
3. Status: Proven. [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) exhaustively enumerates all parity-even scalar contractions from the exact current primitive blocks through `Delta<=4`.
4. Status: Proven. [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) reduces all total-derivative, lower-order-EOM, and Cayley-Hamilton-reducible classes explicitly.
5. Status: Proven. [`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py) shows that the corrected seven survivors are linearly independent as operators for the exact current primitive set.
6. Status: Conjectural. Invoke the conditional collapse lemma to convert that finite normal-form basis into finitely many sensitivity coordinates and Wilson coefficients.
7. Status: Conjectural. Upgrade from the exact current primitive set to a physically justified minimal sector only after the adequacy burden in [`primitive-set-adequacy.md`](primitive-set-adequacy.md) is resolved.

## Exact Current Status

- Status: Proven. The previous five-element `Delta<=4` target list was incomplete.
- Status: Proven. The explicit omitted surviving operators are `divE2` and `mixedGradE2`.
- Status: Proven. Within the exact current primitive set, the corrected seven-scalar list is linearly independent.
- Status: Proven. The next explicit obstruction is no longer internal to that seven-scalar basis.
- Status: Proven. The active obstruction is primitive-set adequacy: a one-family magnetic extension produces the new survivor `B2`.

## Dependencies

- Status: Imported from prior work. [`../lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md)
- Status: Imported from prior work. [`../lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md)
- Status: Conjectural. [`../lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md)
- Status: Proven. [`../lemmas/05-finite-basis-closure.md`](../lemmas/05-finite-basis-closure.md)
- Status: Conjectural. [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md)
- Status: Conjectural. [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md)
- Status: Proven. [`../lemmas/08-mixed-time-derivative-audit.md`](../lemmas/08-mixed-time-derivative-audit.md)
- Status: Proven. [`../lemmas/09-survivor-independence-delta4.md`](../lemmas/09-survivor-independence-delta4.md)
- Status: Conjectural. [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md)
- Status: Conjectural. [`power-counting.md`](power-counting.md)
- Status: Proven. [`primitive-catalog.md`](primitive-catalog.md)
- Status: Proven. [`primitive-set-adequacy.md`](primitive-set-adequacy.md)
- Status: Proven. [`reduction-rules.md`](reduction-rules.md)

## Failure Triggers

- Status: Counterexample candidate. An explicit `chi_A` state invalidates the no-state reduction step.
- Status: Counterexample candidate. Nonanalytic activation invalidates the analytic collapse step.
- Status: Counterexample candidate. Hereditary couplings invalidate locality before basis closure is even applied.
- Status: Counterexample candidate. Infinite or uncontrolled external field content invalidates the fixed-order closure theorem candidate.
- Status: Counterexample candidate. Any additional `Delta<=4` scalar contraction from the current primitive set that evades the explicit enumeration/reduction pipeline would be a direct obstruction to the current M4 path.
- Status: Counterexample candidate. Any physically admissible primitive family outside the exact current primitive set that generates a new `Delta<=4` survivor obstructs promotion to a physically justified minimal-sector theorem.
