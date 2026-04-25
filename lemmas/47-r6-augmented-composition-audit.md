# Lemma 47: Reven6+ Augmented Composition Audit

## Statement

- Status: Proven. Work at fixed order `\Delta \le 4` with the currently audited family classes `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}` only.
- Status: Proven. Impose the current audited thresholds:
  `w_X \ge 4` unless `EX = 0`,
  `w_S \ge 5` or exclude bare `S`,
  `w_D \ge 3` or remove the mixed derivative witnesses,
  `w_V \ge 3` or exclude/absorb the primitive vector family,
  `w_T \ge 3` or exclude/absorb the primitive rank-3 STF family,
  `w_Q \ge 3` or exclude/absorb the primitive rank-4 STF family,
  `w_U \ge 3` or exclude/absorb the primitive rank-5 STF family,
  and `w_Z \ge 3` or exclude/absorb the primitive rank-6 STF family.
- Status: Proven. Then every `Reven6+`-containing triple, quadruple, quintuple, six-family combination, seven-family combination, and the full enlarged audited set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}`, leaves exactly the baseline electric survivor list and no new cross-family witness.

## Triple Layer

- Status: Proven. The `Reven6+`-containing triples are exactly the `21` sets `(Reven6+, A, B)` with distinct `A, B` chosen from `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
- Status: Proven. Every such triple is audited explicitly by [`../symbolic/composition_attack_delta4.py`](../symbolic/composition_attack_delta4.py) and every such triple is sufficient.

## Quadruple Layer

- Status: Proven. The `Reven6+`-containing quadruples are exactly the `35` sets `(Reven6+, A, B, C)` with distinct `A, B, C` chosen from `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
- Status: Proven. Every such quadruple is audited explicitly and every such quadruple is sufficient.

## Quintuple Layer

- Status: Proven. The `Reven6+`-containing quintuples are exactly the `35` sets `(Reven6+, A, B, C, D)` with distinct `A, B, C, D` chosen from `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
- Status: Proven. Every such quintuple is audited explicitly and every such quintuple is sufficient.

## Six-Family, Seven-Family, And Full-Set Layer

- Status: Proven. The `Reven6+`-containing six-family combinations are exactly the `21` sets `(Reven6+, A, B, C, D, E)` with distinct `A, B, C, D, E` chosen from `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
- Status: Proven. Every such six-family set is audited explicitly and every such six-family set is sufficient.
- Status: Proven. The `Reven6+`-containing seven-family combinations are exactly the `7` sets `(Reven6+, A, B, C, D, E, F)` with distinct `A, B, C, D, E, F` chosen from `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}`.
- Status: Proven. Every such seven-family set is audited explicitly and every such seven-family set is sufficient.
- Status: Proven. The full enlarged audited set `(R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+)` is jointly sufficient at `\Delta \le 4`.
- Status: Proven. First explicit post-Reven6+ cross-family obstruction: none found.

## Why The Augmented Layer Closes

1. Status: Proven. The rank-6 threshold already removes every `Reven6+`-containing scalar from the `\Delta \le 4` window.
2. Status: Proven. The remaining audited families were already jointly sufficient before `Reven6+` was added, as recorded in [`40-r5-pairwise-composition-audit.md`](40-r5-pairwise-composition-audit.md) and [`41-r5-augmented-composition-audit.md`](41-r5-augmented-composition-audit.md).
3. Status: Proven. Therefore adjoining the thresholded rank-6 family does not create a new pairwise, higher-order, or full-set cross-family survivor at `\Delta \le 4`.

## Boundary

- Status: Proven. This lemma is only an audited-set composition theorem for the eight currently audited family classes.
- Status: Proven. It does not imply MVP-envelope completeness by itself.
- Status: Proven. After this closure, the live positive bottleneck moves back to family-envelope completeness and the next smallest unaudited family gate.
