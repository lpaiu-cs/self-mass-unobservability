# Lemma 43: Rank-4 Exhaustiveness Check

## Scope

- Status: Note. This lemma compares the manual `Reven4+` survivor bookkeeping against the exhaustive generator introduced in [`../docs/high-rank-audit-methodology.md`](../docs/high-rank-audit-methodology.md).
- Status: Proven. The comparison uses the current primitive-family definition of a genuine parity-even fully symmetric trace-free primitive family `Q_{ijkl}` together with the baseline electric sector at `\Delta \le 4`.
- Status: Proven. Trace descendants reducible to lower-rank audited families and derivative-generated rank-4 descendants are still excluded from the primitive-family definition itself; the present issue is exhaustiveness inside that definition, not family misdefinition.

## Exhaustive Result

- Status: Proven. [`../symbolic/high_rank_family_enumerator.py`](../symbolic/high_rank_family_enumerator.py) generates `22` rank-4 candidate-survivor labels in the current normal-form candidate layer.
- Status: Proven. [`../symbolic/high_rank_diff_report.py`](../symbolic/high_rank_diff_report.py) shows that the current manual `Reven4+` bookkeeping matches only `15` of those `22` labels.
- Status: Proven. Therefore the current manual rank-4 survivor list is not exhaustive.

## First Omitted Contraction

- Status: Proven. The first omitted mixed cubic contraction class is `EEQ`, with signature `(\,E,\ E,\ Q\,)` and unique contraction pattern
  `('pure', ('E', 'E', 'Q'), (0, 2, 2))`.
- Status: Proven. This `EEQ` class exists under the current primitive-family definition because two rank-2 electric tensors can saturate a rank-4 STF tensor without using any trace descendant or derivative-generated block.
- Status: Proven. Therefore the feared `EEQ`-type mixed cubic omission is real under the current rank-4 primitive-family definition.

## Full Omitted Rank-4 Set

- Status: Proven. The exhaustive rank-4 audit finds the following omitted candidate labels:
  `\{EEQ,\ Q3,\ EEDtQ,\ EEEQ,\ EQQQ_1,\ EQQQ_2,\ GradEGradQ\}`.
- Status: Proven. `Q3` is a genuine omitted pure-self cubic invariant at weight `3`.
- Status: Proven. `EEDtQ` is the omitted mixed time-derivative normal-form candidate at weight `4`.
- Status: Proven. `EEEQ`, `EQQQ_1`, and `EQQQ_2` are omitted quartic mixed STF contractions at weight `4`.
- Status: Proven. `GradEGradQ` is an omitted mixed gradient invariant at weight `4`.

## Theorem-Layer Effect

- Status: Proven. The obstruction-class verdict is unchanged: the first self witness remains `Q2` at weight `2`.
- Status: Proven. The threshold type is unchanged: the rank-4 family is still self-only sharp because `W_{\mathrm{self}} = 2 w_Q` remains below the first mixed-witness weight.
- Status: Proven. What fails is not the rank-4 obstruction-class verdict but the claimed exhaustiveness of the manual rank-4 candidate-survivor bookkeeping.

## Immediate Consequence

- Status: Proven. The live theorem bottleneck cannot move on to `Reven6+` yet.
- Status: Proven. The first explicit omitted high-rank contraction `EEQ` is now the new bounded bottleneck that must be resolved before any further family-envelope march upward.
