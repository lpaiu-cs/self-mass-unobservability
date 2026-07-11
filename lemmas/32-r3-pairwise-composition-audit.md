# Lemma 32: Rodd+ Pairwise Composition Audit

## Scope

- Status: Note. This lemma audits the four new pairwise combinations that arise once the genuine primitive rank-3 STF class `Rodd+` is added to the already audited threshold set at `\Delta \le 4`.
- Status: Proven. The imposed thresholds are:
  `R2`: `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness,
  `R0a`: `w_S \ge 5` or explicit exclusion of bare `S`,
  `R0b`: `w_D \ge 3` or explicit removal of the mixed derivative witnesses,
  `R1`: `w_V \ge 3` or explicit exclusion or absorption of the primitive vector family,
  `Rodd+`: `w_T \ge 3` or explicit exclusion or absorption of the primitive rank-3 STF family.
- Status: Proven. The rank-3 family counted here is a genuine primitive parity-even fully symmetric trace-free family `T_{ijk}`; derivative-generated rank-3 blocks such as `\nabla E`, `\nabla B`, `\nabla V`, or scalar-derivative descendants are not counted as a new primitive-family admission.

## Pair `(Rodd+, R2)`

- Status: Proven. With `w_T \ge 3`, the first self rank-3 witness `T2` is pushed to weight `6`, while the first mixed rank-3 witness `ETT` is pushed to weight `7`, so no `Rodd+`-containing scalar enters the `\Delta \le 4` window.
- Status: Proven. With `w_X \ge 4`, no `R2`-containing scalar enters the `\Delta \le 4` window either.
- Status: Proven. Therefore the pair `(Rodd+, R2)` leaves exactly the baseline electric survivor list and no new cross-family witness.
- Status: Proven. Smallest cross-family surviving witness: none.

## Pair `(Rodd+, R0a)`

- Status: Proven. With `w_T \ge 3`, every genuine primitive rank-3 witness lies above `\Delta = 4`.
- Status: Proven. With `w_S \ge 5` or explicit exclusion of bare `S`, no bare-scalar operator remains inside the `\Delta \le 4` window.
- Status: Proven. Therefore the pair `(Rodd+, R0a)` is sufficient at `\Delta \le 4`.
- Status: Proven. Smallest cross-family surviving witness: none.

## Pair `(Rodd+, R0b)`

- Status: Proven. With `w_T \ge 3`, no `Rodd+`-containing scalar survives at or below `\Delta = 4`.
- Status: Proven. With `w_D \ge 3`, the derivative-only family contributes only reducible singleton or lower-order-EOM candidates such as `DtS`, `Dt2S`, and `aGradS` at or below `\Delta = 4`.
- Status: Proven. Therefore the pair `(Rodd+, R0b)` leaves no new survivor beyond the baseline electric sector.
- Status: Proven. Smallest cross-family surviving witness: none.

## Pair `(Rodd+, R1)`

- Status: Proven. With `w_T \ge 3`, every genuine primitive rank-3 witness lies above `\Delta = 4`.
- Status: Proven. With `w_V \ge 3`, every genuine primitive vector witness also lies above `\Delta = 4`.
- Status: Proven. Therefore the pair `(Rodd+, R1)` leaves exactly the baseline electric survivor list and no new cross-family witness.
- Status: Proven. Smallest cross-family surviving witness: none.

## Boundary

- Status: Proven. Absence of new Rodd+-containing pairwise obstructions does not automatically prove triple, quadruple, or five-family sufficiency.
- Status: Proven. The higher-order audited combinations are recorded separately in [`33-r3-augmented-composition-audit.md`](33-r3-augmented-composition-audit.md).
