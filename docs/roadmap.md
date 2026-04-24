# Roadmap

Status: Conjectural. The repository roadmap is organized around theorem progress first and loophole progress second.

## M1: Free-Fall Sensitivity-Collapse Scaffold

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M1.1 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) | Precise theorem statement and proof skeleton exist. |
| M1.2 | Imported from prior work | [`lemmas/01-internal-structure-no-go.md`](../lemmas/01-internal-structure-no-go.md) | Imported obstruction to literal internal self-unobservability is recorded. |
| M1.3 | Imported from prior work | [`lemmas/02-com-decoupling.md`](../lemmas/02-com-decoupling.md) | Imported COM decoupling lemma is recorded. |
| M1.4 | Conjectural | [`lemmas/03-worldline-reduction.md`](../lemmas/03-worldline-reduction.md) | Generic local worldline action note exists. |
| M1.5 | Counterexample candidate | [`counterexamples/chi-state/README.md`](../counterexamples/chi-state/README.md) | At least one explicit loophole model is written down. |
| M1.6 | Proven | [`symbolic/worldline_expand.py`](../symbolic/worldline_expand.py) and [`symbolic/checks/test_symbolic.py`](../symbolic/checks/test_symbolic.py) | Symbolic scaffold runs without error. |

## M2: Non-Circular Fixed-Order Theorem Candidate

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M2.1 | Conjectural | [`docs/conditional-collapse-lemma.md`](conditional-collapse-lemma.md) | Conditional collapse is isolated from basis closure. |
| M2.2 | Conjectural | [`docs/power-counting.md`](power-counting.md) | Fixed-order counting rule is explicit and does not assume finite basis closure. |
| M2.3 | Proven | [`lemmas/05-finite-basis-closure.md`](../lemmas/05-finite-basis-closure.md) | Abstract candidate-set finiteness is separated from the physical completeness burden. |
| M2.4 | Proven | [`../symbolic/enumerate_basis.py`](../symbolic/enumerate_basis.py) | The chosen-order candidate monomial list enumerates without error. |
| M2.5 | Counterexample candidate | [`../counterexamples/chi-state/README.md`](../counterexamples/chi-state/README.md) | The smallest loophole model includes an adiabatic-collapse analysis. |

## M3: Delta<=4 Normal-Form Completeness Attack

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M3.1 | Proven | [`docs/primitive-catalog.md`](primitive-catalog.md) | Exact `Delta<=4` primitive and candidate operator content is fixed. |
| M3.2 | Proven | [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) | Total-derivative, lower-order EOM, and algebraic reductions run without error. |
| M3.3 | Proven | [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md) | The old five-element target is recognized as incomplete rather than silently assumed. |

## M4: Contraction-Level Exhaustiveness Audit

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M4.1 | Proven | [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) | All parity-even scalar contractions from the exact current primitive blocks through `Delta<=4` are enumerated. |
| M4.2 | Proven | [`docs/reduction-rules.md`](reduction-rules.md) | Every allowed reduction rule is listed explicitly. |
| M4.3 | Conjectural | [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md) | The gradient sector is audited and the surviving basis is identified explicitly. |
| M4.4 | Proven | [`../lemmas/08-mixed-time-derivative-audit.md`](../lemmas/08-mixed-time-derivative-audit.md) | Mixed `E-E-D_tau E` terms are audited explicitly. |
| M4.5 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) | The corrected seven-element `Delta<=4` normal-form path is stated explicitly. |

## M5: Independence And Primitive-Set Adequacy

| Item | Status | Deliverable | Exit condition |
| --- | --- | --- | --- |
| M5.1 | Proven | [`../lemmas/09-survivor-independence-delta4.md`](../lemmas/09-survivor-independence-delta4.md) | The corrected seven survivors are checked for linear independence rather than inferred from non-reducibility. |
| M5.2 | Proven | [`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py) | Exact polynomial coefficient rank confirms the current seven-scalar basis is linearly independent. |
| M5.3 | Proven | [`primitive-set-adequacy.md`](primitive-set-adequacy.md) | The exact-current-set theorem is separated sharply from any physically justified minimal-sector claim. |
| M5.4 | Proven | [`../symbolic/primitive_family_attack.py`](../symbolic/primitive_family_attack.py) | One physically motivated primitive-family extension is audited explicitly. |
| M5.5 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) and [`failure-ledger.md`](failure-ledger.md) | The main remaining risk is written uniformly as primitive-set adequacy, not survivor dependence. |

## Near-Term Sequence

1. Status: Conjectural. Decide whether the magnetic tidal family `B_ij` should be excluded by a new explicit assumption or admitted into the minimal-sector theorem target.
2. Status: Conjectural. If `B_ij` is admitted, rerun the normal-form and collapse pipeline on the enlarged primitive catalog rather than preserving the older electric-only theorem statement.
3. Status: Conjectural. Decide whether any extra gradient identity should be introduced explicitly, rather than silently collapsing `divE2` or `mixedGradE2`.
4. Status: Counterexample candidate. Promote the dominant loophole class only if it produces a genuine failure of the explicit enumeration or adequacy path.
5. Status: Conjectural. Revisit the scalar `s_A` corollary only after the higher-dimensional sensitivity manifold statement is sharp.

## Explicit Non-Goals For M1

- Status: Proven. No new LLR, MLRS, PEP, Nutimo, or TOA pipeline work belongs in M1.
- Status: Proven. No clock-sector theorem is attempted in M1.
- Status: Proven. No build or runtime chores count as theorem progress.

## Classification Rule

- `theorem progress`: theorem draft, lemma draft, assumption cleanup, symbolic reduction check
- `loophole progress`: explicit counterexample model, exact failed proof step, minimal violated assumption
