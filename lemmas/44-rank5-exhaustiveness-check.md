# Lemma 44: Rank-5 Exhaustiveness Check

## Scope

- Status: Proven. This lemma compares the manual `Rodd5+` survivor bookkeeping against the exhaustive generator introduced in [`../docs/high-rank-audit-methodology.md`](../docs/high-rank-audit-methodology.md).
- Status: Proven. The comparison is restricted to the current primitive-family definition used in the rank-5 audit: a genuine parity-even fully symmetric trace-free primitive family `U_{ijklm}` plus the baseline electric sector at `\Delta \le 4`.

## Exhaustive Result

- Status: Proven. [`../symbolic/high_rank_family_enumerator.py`](../symbolic/high_rank_family_enumerator.py) generates exactly `16` rank-5 candidate-survivor labels.
- Status: Proven. [`../symbolic/high_rank_diff_report.py`](../symbolic/high_rank_diff_report.py) shows `16/16` matches between the exhaustive generated list and the current manual `Rodd5+` survivor list.
- Status: Proven. No omitted rank-5 contraction class appears in the exhaustive comparison.

## Matched Rank-5 List

- Status: Proven. The exhaustive rank-5 candidate-survivor list is
  `\{U2,\ EUU,\ dotU2,\ EUDtU,\ E2U2,\ E2U2_{mixed,1},\ E2U2_{mixed,2},\ E2U2_{mixed,3},\ divU2,\ gradU2,\ mixedGradU2,\ U2^2,\ U4_{balanced},\ U4_{bridge},\ U4_{chain},\ U4_{tetra}\}`.
- Status: Proven. Therefore the current manual `Rodd5+` survivor bookkeeping is complete at the level of candidate contraction classes, even though the raw rank-5 list is still not a corrected basis statement because of the already recorded null relations.

## Consequence

- Status: Proven. The rank-5 obstruction-class verdict itself is unchanged.
- Status: Proven. The rank-5 threshold formula `w_U \ge 3` is unchanged.
- Status: Proven. The rank-5 bookkeeping issue is closed for exhaustiveness, not only for internal consistency.
