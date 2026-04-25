# Lemma 30: Rank-3 Family Admission

## Scope

- Status: Proven. This lemma audits unsuppressed admission of a genuine local parity-even rank-3 family in the MVP free-fall sector at `\Delta \le 4`.
- Status: Proven. The chosen representative is the fully symmetric trace-free primitive family `T_{ijk}` defined in [`../docs/rank3-family-ordering.md`](../docs/rank3-family-ordering.md).
- Status: Proven. Derivative-generated rank-3 blocks already attached to audited families are excluded from the primitive-family definition and are not double-counted here.

## First Witnesses

- Status: Proven. The first self witness is `T2 = T_{ijk} T_{ijk}` at weight `2`.
- Status: Proven. The first mixed witness against the current enlarged audited-set baseline is `ETT = E_{ab} T_{acd} T_{bcd}` at weight `3`.
- Status: Proven. Therefore unsuppressed admission of the genuine rank-3 family creates a new low-order witness before any new composition rescue is attempted.

## Extended `\Delta \le 4` Survivor List

- Status: Proven. [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py) finds the new surviving labels
  `\{T2,\ ETT,\ dotT2,\ ETDtT,\ E2T2,\ E2T2_mixed_1,\ E2T2_mixed_2,\ E2T2_mixed_3,\ divT2,\ gradT2,\ mixedGradT2,\ T2^2,\ T4_chain,\ T4_tetra\}`
  beyond the enlarged audited-set baseline.
- Status: Proven. [`../symbolic/r3_survivor_rank_check.py`](../symbolic/r3_survivor_rank_check.py) shows that the raw R3-extended survivor list is not linearly independent; it has nullity `2`.
- Status: Proven. The first exact dependence relation already extracted is
  `-E2T2/2 + 2 E2T2_mixed_1 - E2T2_mixed_2 + E2T2_mixed_3 = 0`.
- Status: Proven. This linear-dependence correction does not rescue minimal-sector uniqueness, because the no-go already follows from the lower-weight witnesses `T2` and `ETT`.

## Theorem Layer Obstructed

- Status: Proven. The theorem layer obstructed by this audit is promotion of the enlarged audited-set result `{R2, R0a, R0b, R1}` to MVP-envelope sufficiency.
- Status: Proven. The direct obstruction witness is `T2`.
- Status: Proven. Therefore `Rodd+` is a genuine next obstruction class under the current MVP assumptions.

## Boundary

- Status: Proven. This lemma does not yet re-close audited-set composition after admitting the rank-3 threshold.
- Status: Proven. The next positive-step burden is a new enlarged audited-set composition audit that includes `Rodd+`.
