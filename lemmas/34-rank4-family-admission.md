# Lemma 34: Rank-4 Family Admission

## Scope

- Status: Note. This lemma audits unsuppressed admission of a genuine local parity-even rank-4 family in the MVP free-fall sector at `\Delta \le 4`.
- Status: Proven. The chosen representative is the fully symmetric trace-free primitive family `Q_{ijkl}` defined in [`../docs/rank4-family-ordering.md`](../docs/rank4-family-ordering.md).
- Status: Proven. Trace descendants reducible to lower-rank audited families and derivative-generated rank-4 blocks already attached to audited families are excluded from the primitive-family definition and are not double-counted here.

## First Witnesses

- Status: Proven. The first self witness is `Q2 = Q_{ijkl} Q_{ijkl}` at weight `2`.
- Status: Proven. The original manual audit already exhibited one weight-`3` mixed witness, `EQQ = E_{ab} Q_{acde} Q_{bcde}`.
- Status: Proven. The later exhaustive check in [`43-rank4-exhaustiveness-check.md`](43-rank4-exhaustiveness-check.md) shows that `EEQ` is an equally low weight-`3` mixed witness omitted from the manual bookkeeping.
- Status: Proven. Therefore unsuppressed admission of the genuine rank-4 family creates a new low-order witness before any new composition rescue is attempted.

## Extended `\Delta \le 4` Survivor List

- Status: Proven. [`../symbolic/r4_sector_delta4.py`](../symbolic/r4_sector_delta4.py) records the provisional manual survivor bookkeeping
  `\{Q2,\ EQQ,\ dotQ2,\ EQDtQ,\ E2Q2,\ E2Q2_mixed_1,\ E2Q2_mixed_2,\ E2Q2_mixed_3,\ divQ2,\ gradQ2,\ mixedGradQ2,\ Q2^2,\ Q4_bridge,\ Q4_chain,\ Q4_tetra\}`
  beyond the enlarged audited-set baseline.
- Status: Proven. [`43-rank4-exhaustiveness-check.md`](43-rank4-exhaustiveness-check.md) shows that this manual list is not exhaustive; the omitted labels are
  `\{EEQ,\ Q3,\ EEDtQ,\ EEEQ,\ EQQQ_1,\ EQQQ_2,\ GradEGradQ\}`.
- Status: Proven. [`../symbolic/r4_survivor_rank_check.py`](../symbolic/r4_survivor_rank_check.py) therefore applies only to the pre-exhaustiveness manual subset, not to a corrected exhaustive rank-4 basis.
- Status: Proven. On that manual subset, the raw R4 survivor bookkeeping is not linearly independent; it has raw count `22`, rank `19`, and nullity `3`.
- Status: Proven. The first exact dependence relation already extracted is
  `-E2Q2/2 + 2 E2Q2_mixed_1 - E2Q2_mixed_2 + E2Q2_mixed_3 = 0`.
- Status: Proven. Two additional exact quartic self relations are already verified:
  `3 Q2^2 / 10 + Q4_bridge - 8 Q4_chain / 5 = 0`
  and
  `-Q2^2 / 5 + 2 Q4_chain / 5 + Q4_tetra = 0`.
- Status: Proven. The raw R4 bookkeeping is therefore doubly non-final: it is neither exhaustive nor a corrected basis statement.
- Status: Note. This linear-dependence correction does not rescue minimal-sector uniqueness, because the no-go already follows from the lower-weight witnesses `Q2` and `EQQ`.

## Theorem Layer Obstructed

- Status: Proven. The theorem layer obstructed by this audit is promotion of the enlarged audited-set result `{R2, R0a, R0b, R1, Rodd+}` to MVP-envelope sufficiency.
- Status: Proven. The direct obstruction witness is `Q2`.
- Status: Proven. Therefore `Reven4+` is a genuine next obstruction class under the current MVP assumptions.

## Boundary

- Status: Note. This lemma does not yet re-close audited-set composition after admitting the rank-4 threshold.
- Status: Proven. The next positive-step burden is a new enlarged audited-set composition audit that includes `Reven4+`.
