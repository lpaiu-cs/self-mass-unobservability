# Lemma 42: STF Rank-`L` Admission Attempt

## Scope

- Status: Proven. This lemma abstracts the common admission structure behind the audited genuine parity-even fully symmetric trace-free primitive families of ranks `L = 3, 4, 5`.
- Status: Proven. The primitive family is denoted `Y_L`, with `L \ge 3`, parity even, no trace descendants counted as primitive, and no derivative-generated descendants counted as new primitive-family admissions.
- Status: Proven. The present statement is class-limited to the STF tower only; it is not a theorem for all higher-rank tensor families.

## Structural Statements That Do Hold

- Status: Proven. A linear scalar built from a single STF block `Y_L` does not exist for any `L \ge 1`, because any scalar contraction would require an internal trace.
- Status: Proven. Therefore the quadratic norm
  `Y2 = Y_{i_1 \dots i_L} Y_{i_1 \dots i_L}`
  is the first self witness candidate for every genuine STF family.
- Status: Proven. This quadratic self contraction exists for every `L \ge 1`, hence in particular for every `L \ge 3`.
- Status: Proven. A mixed quadratic scalar `E Y_L` exists only when `L = 2`, because a two-block full contraction requires equal ranks.
- Status: Proven. This is the exact structural place where the rank-2 STF family is exceptional.
- Status: Proven. For every `L \ge 3`, the mixed cubic contraction
  `EYY = E_{ab} Y_{aA} Y_{bA}`
  exists by contracting one electric index into each STF block and pairing the remaining `L-1` STF indices between the two `Y_L` copies.

## Lowest-Weight Admission Layer

- Status: Proven. For the unsuppressed audited admission weight `w_Y = 1`, the self witness `Y2` sits at weight `2`.
- Status: Proven. For the unsuppressed audited admission weight `w_Y = 1`, the cubic mixed witness `EYY` sits at weight `3`.
- Status: Proven. No parity-even scalar with a genuine STF family insertion exists at weight `1`.
- Status: Proven. No mixed quadratic scalar exists for `L \ge 3`.
- Status: Proven. Therefore, across the audited STF ranks, there is no parity-even scalar with a genuine `Y_L` insertion below weight `2`, and no mixed scalar below weight `3`.

## Exact Failing Step Of The Tower Theorem

- Status: Conjectural. The attempted stronger tower statement would say that the first mixed witness layer is exhausted by `EYY` for every `L \ge 3`.
- Status: Proven. That stronger step fails already at rank `L = 4`.
- Status: Proven. The explicit counterexample inside the audited STF tower is
  `EEQ`, a genuine parity-even mixed cubic contraction of one rank-4 STF block with two electric tensors.
- Status: Proven. Therefore the admission theorem does not close as a single universal STF-tower statement with unique mixed witness `EYY`.

## Result

- Status: Proven. The common audited STF structure is:
  `Y2` always exists,
  `EYY` always exists,
  and rank `L = 2` is exceptional because `EY` already exists.
- Status: Proven. The attempted universal higher-rank conclusion fails at the first audited even-rank member `L = 4` because `EEQ` survives.
