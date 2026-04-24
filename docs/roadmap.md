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
| M3.3 | Conjectural | [`../lemmas/06-normal-form-completeness-delta4.md`](../lemmas/06-normal-form-completeness-delta4.md) | A real Delta<=4 completeness path or explicit obstruction is recorded. |
| M3.4 | Conjectural | [`docs/theorem-A-freefall.md`](theorem-A-freefall.md) | The current burden is written uniformly as `normal-form completeness modulo total derivatives and lower-order equations of motion`. |

## Near-Term Sequence

1. Status: Conjectural. Close `normal-form completeness modulo total derivatives and lower-order equations of motion` for the exact `Delta<=4` primitive catalog.
2. Status: Conjectural. Either prove catalog exhaustiveness at `Delta<=4` or isolate the smallest additional admissible scalar obstruction.
3. Status: Counterexample candidate. Promote the dominant loophole class only if it produces a genuine failure of the `Delta<=4` normal-form target list.
4. Status: Conjectural. Revisit the scalar `s_A` corollary only after the higher-dimensional sensitivity manifold statement is sharp.

## Explicit Non-Goals For M1

- Status: Proven. No new LLR, MLRS, PEP, Nutimo, or TOA pipeline work belongs in M1.
- Status: Proven. No clock-sector theorem is attempted in M1.
- Status: Proven. No build or runtime chores count as theorem progress.

## Classification Rule

- `theorem progress`: theorem draft, lemma draft, assumption cleanup, symbolic reduction check
- `loophole progress`: explicit counterexample model, exact failed proof step, minimal violated assumption
