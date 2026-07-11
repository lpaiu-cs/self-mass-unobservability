# Lemma 40: Rodd5+ Pairwise Composition Audit

## Statement

- Status: Proven. Work at fixed order `\Delta \le 4` with the currently audited family classes only.
- Status: Proven. Impose the current audited thresholds:
  `w_X \ge 4` unless `EX = 0`,
  `w_S \ge 5` or exclude bare `S`,
  `w_D \ge 3` or remove the mixed derivative witnesses,
  `w_V \ge 3` or exclude/absorb the primitive vector family,
  `w_T \ge 3` or exclude/absorb the primitive rank-3 STF family,
  `w_Q \ge 3` or exclude/absorb the primitive rank-4 STF family,
  and `w_U \ge 3` or exclude/absorb the primitive rank-5 STF family.
- Status: Proven. Then each new pairwise audited combination involving `Rodd5+`,
  `(Rodd5+, R2)`, `(Rodd5+, R0a)`, `(Rodd5+, R0b)`, `(Rodd5+, R1)`, `(Rodd5+, Rodd+)`, and `(Rodd5+, Reven4+)`,
  leaves exactly the baseline electric survivor list and no new cross-family witness.

## Audit Table

| Pair | Result | Smallest surviving cross-family witness |
| --- | --- | --- |
| `(Rodd5+, R2)` | Proven sufficient | None |
| `(Rodd5+, R0a)` | Proven sufficient | None |
| `(Rodd5+, R0b)` | Proven sufficient | None |
| `(Rodd5+, R1)` | Proven sufficient | None |
| `(Rodd5+, Rodd+)` | Proven sufficient | None |
| `(Rodd5+, Reven4+)` | Proven sufficient | None |

## Why The Pairwise Layer Closes

1. Status: Proven. The sharp rank-5 threshold `w_U \ge 3` pushes the first self witness `U2` to weight `6`, outside the `\Delta \le 4` window.
2. Status: Proven. The first mixed witness `EUU` is pushed to weight `7`, also outside the `\Delta \le 4` window.
3. Status: Proven. Every other explicit rank-5 scalar class in [`../symbolic/r5_sector_delta4.py`](../symbolic/r5_sector_delta4.py) has weight at least that large once the threshold is imposed.
4. Status: Proven. Therefore no `Rodd5+`-containing scalar survives inside the audited pairwise window at `\Delta \le 4`.

## Boundary

- Status: Note. This lemma is only a pairwise audited-set composition result.
- Status: Proven. It does not prove triple, quadruple, quintuple, six-family, or full-set sufficiency by itself.
- Status: Proven. It does not re-interpret the raw `23`-element rank-5 survivor list as a corrected basis statement.
