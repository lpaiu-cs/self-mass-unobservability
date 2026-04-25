# Lemma 44: Rank-6 Family Admission

## Scope

- Status: Proven. This lemma audits unsuppressed admission of a genuine parity-even fully symmetric trace-free rank-6 primitive family `Z_{ijklmn}` in the MVP free-fall sector at `\Delta \le 4`.
- Status: Proven. The audit uses the exhaustive high-rank methodology from the start via [`../symbolic/r6_sector_delta4.py`](../symbolic/r6_sector_delta4.py), not manual survivor bookkeeping.
- Status: Proven. Trace descendants reducible to lower-rank audited families and derivative-generated rank-6 descendants are excluded from the primitive-family definition itself.

## Exhaustive Candidate Layer

- Status: Proven. [`../symbolic/r6_sector_delta4.py`](../symbolic/r6_sector_delta4.py) generates `22` new rank-6 candidate-survivor labels beyond the electric baseline.
- Status: Proven. The smallest new label is `Z2`.
- Status: Proven. The exhaustive first mixed layer contains exactly one audited weight-`3` mixed label: `EZZ`.

## First Witnesses

- Status: Proven. The first self witness is the quadratic norm `Z2` at weight `2`.
- Status: Proven. The first mixed witness is `EZZ` at weight `3`.
- Status: Proven. A cubic self invariant `Z3` also survives at weight `3`, but it does not precede `Z2`.

## Even-Rank Mixed-Class Test

- Status: Proven. No mixed quadratic scalar `EZ` exists because the electric rank `2` block cannot fully contract a genuine rank-6 STF block.
- Status: Proven. No `EEZ`-type mixed cubic survives under the present primitive-family definition.
- Status: Proven. The structural reason is that two rank-2 electric tensors cannot saturate a rank-6 STF block without introducing forbidden self-contractions or trace-descended structure.
- Status: Proven. The exhaustive rank-6 generator confirms that there is no additional even-rank mixed class at the same first mixed order or lower effective order beyond `EZZ`.

## Raw-List Bookkeeping

- Status: Conjectural. [`../symbolic/r6_survivor_rank_check.py`](../symbolic/r6_survivor_rank_check.py) finds a raw rank-6 candidate list with sample-stable rank `16` out of `22`.
- Status: Conjectural. The first revalidated null relation is
  `E2Z2 + 2 E2Z2_{mixed,1} - 2 E2Z2_{mixed,2} - 4 E2Z2_{mixed,3} = 0`.
- Status: Conjectural. Therefore the raw rank-6 list is bookkeeping only and is not presented as a corrected basis statement.

## Theorem-Layer Effect

- Status: Proven. Unsuppressed `Reven6+` admission is a genuine new obstruction class with smallest explicit witness `Z2`.
- Status: Proven. This obstructs promotion of the enlarged audited-set result to MVP-envelope sufficiency unless the rank-6 family is explicitly suppressed or excluded.
- Status: Proven. Because `Reven6+` is now audited and obstructive, the next theorem-layer step is a post-`Reven6+` enlarged audited-set composition re-close rather than an immediate move to the next family gate.
