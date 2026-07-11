# Lemma 22: Pairwise Composition Audit

## Scope

- Status: Note. This lemma audits the three pairwise combinations of the currently audited family classes at `\Delta \le 4`.
- Status: Proven. The imposed thresholds are:
  `R2`: `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness,
  `R0a`: `w_S \ge 5` or explicit exclusion of bare `S`,
  `R0b`: `w_D \ge 3` or explicit removal of the mixed derivative witnesses.

## Pair `(R2, R0a)`

- Status: Proven. With `w_X \ge 4` and `w_S \ge 5`, any scalar involving `X` or `S` lies above `\Delta = 4` or fails scalar contraction at that order.
- Status: Proven. The resulting survivor list is exactly the baseline electric survivor list.
- Status: Proven. Therefore the pair `(R2, R0a)` is sufficient at `\Delta \le 4`.
- Status: Proven. Smallest cross-family surviving witness: none.

## Pair `(R2, R0b)`

- Status: Proven. With `w_X \ge 4`, no `X`-containing scalar survives at or below `\Delta = 4`.
- Status: Proven. With `w_D \ge 3`, the only derivative-only candidates remaining at or below `\Delta = 4` are reducible singleton/EOM terms such as `DtS`, `Dt2S`, and `aGradS`.
- Status: Proven. Therefore the pair `(R2, R0b)` leaves no new survivor beyond the baseline electric sector.
- Status: Proven. Smallest cross-family surviving witness: none.

## Pair `(R0a, R0b)`

- Status: Proven. With `w_S \ge 5`, bare-scalar operators are pushed outside the `\Delta \le 4` window.
- Status: Proven. With `w_D \ge 3`, the derivative-only candidates remaining at or below `\Delta = 4` reduce away as in the single-family audit.
- Status: Proven. Therefore the pair `(R0a, R0b)` is sufficient at `\Delta \le 4`.
- Status: Proven. Smallest cross-family surviving witness: none.

## Boundary

- Status: Proven. Absence of pairwise obstructions does not automatically prove triple sufficiency.
- Status: Proven. The triple audit is recorded separately in [`23-three-family-composition-audit.md`](23-three-family-composition-audit.md).
