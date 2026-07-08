# Lemma 29: R1-Augmented Composition Audit

## Scope

- Status: Note. This lemma audits the R1-containing triple combinations `(R1, R2, R0a)`, `(R1, R2, R0b)`, `(R1, R0a, R0b)`, and the full enlarged audited set `(R2, R0a, R0b, R1)` at `\Delta \le 4`.
- Status: Proven. The imposed threshold set is the current audited one:
  `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness,
  `w_S \ge 5` or explicit exclusion of bare `S`,
  `w_D \ge 3` or explicit removal of the mixed derivative witnesses,
  and `w_V \ge 3` or explicit exclusion or absorption of the primitive vector family.

## Triple `(R1, R2, R0a)`

- Status: Proven. With `w_V \ge 3`, no genuine primitive vector scalar enters the `\Delta \le 4` window.
- Status: Proven. With `w_X \ge 4`, no `R2`-containing scalar enters the `\Delta \le 4` window.
- Status: Proven. With `w_S \ge 5` or explicit exclusion of bare `S`, no bare scalar enters the `\Delta \le 4` window.
- Status: Proven. Therefore the triple `(R1, R2, R0a)` leaves exactly the baseline electric survivor list.

## Triple `(R1, R2, R0b)`

- Status: Proven. With `w_V \ge 3`, no genuine primitive vector scalar enters the `\Delta \le 4` window.
- Status: Proven. With `w_X \ge 4`, no `R2`-containing scalar enters the `\Delta \le 4` window.
- Status: Proven. With `w_D \ge 3`, the derivative-only family contributes only reducible singleton or lower-order-EOM terms at `\Delta \le 4`.
- Status: Proven. Therefore the triple `(R1, R2, R0b)` leaves exactly the baseline electric survivor list.

## Triple `(R1, R0a, R0b)`

- Status: Proven. With `w_V \ge 3`, no genuine primitive vector scalar enters the `\Delta \le 4` window.
- Status: Proven. With `w_S \ge 5` or explicit exclusion of bare `S`, no bare scalar enters the `\Delta \le 4` window.
- Status: Proven. With `w_D \ge 3`, the derivative-only family contributes only reducible singleton or lower-order-EOM terms at `\Delta \le 4`.
- Status: Proven. Therefore the triple `(R1, R0a, R0b)` leaves exactly the baseline electric survivor list.

## Full Enlarged Audited Set `(R2, R0a, R0b, R1)`

- Status: Proven. The simultaneous four-family audit leaves no `R2`-, `R0a`-, `R0b`-, or genuine `R1`-containing survivor at or below `\Delta = 4` beyond the baseline electric sector.
- Status: Proven. The full enlarged audited set therefore leaves exactly the baseline electric survivor list and no new pairwise, triple, or quadruple cross-family survivor.

## Consequence

- Status: Proven. The current audited threshold set is jointly sufficient for minimal-sector uniqueness against the enlarged audited family set `{R2, R0a, R0b, R1}` at `\Delta \le 4`.
- Status: Proven. First explicit post-R1 cross-family obstruction: none found in the enlarged audited set.

## Boundary

- Status: Note. This is not a universal composition theorem for arbitrary family catalogs.
- Status: Proven. It is only a jointly sufficient audited-set result for the enlarged audited family set `{R2, R0a, R0b, R1}` at fixed order `\Delta \le 4`.
