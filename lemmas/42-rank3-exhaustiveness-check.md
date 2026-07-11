# Lemma 42: Rank-3 Exhaustiveness Check

## Scope

- Status: Note. This lemma compares the manual `Rodd+` survivor bookkeeping against the exhaustive generator introduced in [`../docs/high-rank-audit-methodology.md`](../docs/high-rank-audit-methodology.md).
- Status: Proven. The comparison is restricted to the current primitive-family definition used in the rank-3 audit: a genuine parity-even fully symmetric trace-free primitive family `T_{ijk}` plus the baseline electric sector at `\Delta \le 4`.

## Exhaustive Result

- Status: Proven. [`../symbolic/high_rank_family_enumerator.py`](../symbolic/high_rank_family_enumerator.py) generates exactly `14` rank-3 candidate-survivor labels.
- Status: Proven. [`../symbolic/high_rank_diff_report.py`](../symbolic/high_rank_diff_report.py) shows `14/14` matches between the exhaustive generated list and the current manual `Rodd+` survivor list.
- Status: Proven. No omitted rank-3 contraction class appears in the exhaustive comparison.

## Matched Rank-3 List

- Status: Proven. The exhaustive rank-3 candidate-survivor list is
  `\{T2,\ ETT,\ dotT2,\ ETDtT,\ E2T2,\ E2T2_{mixed,1},\ E2T2_{mixed,2},\ E2T2_{mixed,3},\ divT2,\ gradT2,\ mixedGradT2,\ T2^2,\ T4_{chain},\ T4_{tetra}\}`.
- Status: Proven. Therefore the current manual `Rodd+` survivor bookkeeping is complete at the level of candidate contraction classes, even though the raw rank-3 list is still not a corrected basis statement because of the already recorded null relations.

## Consequence

- Status: Proven. The rank-3 obstruction-class verdict itself is unchanged.
- Status: Proven. The rank-3 threshold formula `w_T \ge 3` is unchanged.
- Status: Proven. The rank-3 bookkeeping issue is closed for exhaustiveness, not only for internal consistency.
