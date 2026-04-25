# Lemma 28: R1 Pairwise Composition Audit

## Scope

- Status: Proven. This lemma audits the three new pairwise combinations that arise once the genuine primitive vector class `R1` is added to the already audited threshold set at `\Delta \le 4`.
- Status: Proven. The imposed thresholds are:
  `R2`: `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness,
  `R0a`: `w_S \ge 5` or explicit exclusion of bare `S`,
  `R0b`: `w_D \ge 3` or explicit removal of the mixed derivative witnesses,
  `R1`: `w_V \ge 3` or explicit exclusion or absorption of the primitive vector family.
- Status: Proven. The vector family counted here is a genuine primitive parity-even rank-1 family `V_i`; derivative-generated vectors such as `\nabla_i S` or divergence descendants of audited STF families are not counted as a new primitive-family admission.

## Pair `(R1, R2)`

- Status: Proven. With `w_V \ge 3`, the first self vector witness `V2` is pushed to weight `6`, while the first mixed vector witness `EVV` is pushed to weight `7`, so no `R1`-containing scalar enters the `\Delta \le 4` window.
- Status: Proven. With `w_X \ge 4`, no `R2`-containing scalar enters the `\Delta \le 4` window either.
- Status: Proven. Therefore the pair `(R1, R2)` leaves exactly the baseline electric survivor list and no new cross-family witness.
- Status: Proven. Smallest cross-family surviving witness: none.

## Pair `(R1, R0a)`

- Status: Proven. With `w_V \ge 3`, every genuine primitive vector witness lies above `\Delta = 4`.
- Status: Proven. With `w_S \ge 5` or explicit exclusion of bare `S`, no bare-scalar operator remains inside the `\Delta \le 4` window.
- Status: Proven. Therefore the pair `(R1, R0a)` is sufficient at `\Delta \le 4`.
- Status: Proven. Smallest cross-family surviving witness: none.

## Pair `(R1, R0b)`

- Status: Proven. With `w_V \ge 3`, no `R1`-containing scalar survives at or below `\Delta = 4`.
- Status: Proven. With `w_D \ge 3`, the derivative-only family contributes only reducible singleton or lower-order-EOM candidates such as `DtS`, `Dt2S`, and `aGradS` at or below `\Delta = 4`.
- Status: Proven. Therefore the pair `(R1, R0b)` leaves no new survivor beyond the baseline electric sector.
- Status: Proven. Smallest cross-family surviving witness: none.

## Boundary

- Status: Proven. Absence of new R1-containing pairwise obstructions does not automatically prove triple or quadruple sufficiency.
- Status: Proven. The higher-order audited combinations are recorded separately in [`29-r1-augmented-composition-audit.md`](29-r1-augmented-composition-audit.md).
