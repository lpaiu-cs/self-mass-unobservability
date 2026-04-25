# Lemma 46: Reven6+ Pairwise Composition Audit

## Statement

- Status: Proven. Work at fixed order `\Delta \le 4` with the currently audited family classes only.
- Status: Proven. Impose the current audited thresholds:
  `w_X \ge 4` unless `EX = 0`,
  `w_S \ge 5` or exclude bare `S`,
  `w_D \ge 3` or remove the mixed derivative witnesses,
  `w_V \ge 3` or exclude/absorb the primitive vector family,
  `w_T \ge 3` or exclude/absorb the primitive rank-3 STF family,
  `w_Q \ge 3` or exclude/absorb the primitive rank-4 STF family,
  `w_U \ge 3` or exclude/absorb the primitive rank-5 STF family,
  and `w_Z \ge 3` or exclude/absorb the primitive rank-6 STF family.
- Status: Proven. Then each new pairwise audited combination involving `Reven6+`,
  `(Reven6+, R2)`, `(Reven6+, R0a)`, `(Reven6+, R0b)`, `(Reven6+, R1)`, `(Reven6+, Rodd+)`, `(Reven6+, Reven4+)`, and `(Reven6+, Rodd5+)`,
  leaves exactly the baseline electric survivor list and no new cross-family witness.

## Audit Table

| Pair | Result | Smallest surviving cross-family witness |
| --- | --- | --- |
| `(Reven6+, R2)` | Proven sufficient | None |
| `(Reven6+, R0a)` | Proven sufficient | None |
| `(Reven6+, R0b)` | Proven sufficient | None |
| `(Reven6+, R1)` | Proven sufficient | None |
| `(Reven6+, Rodd+)` | Proven sufficient | None |
| `(Reven6+, Reven4+)` | Proven sufficient | None |
| `(Reven6+, Rodd5+)` | Proven sufficient | None |

## Why The Pairwise Layer Closes

1. Status: Proven. The sharp rank-6 threshold `w_Z \ge 3` pushes the first self witness `Z2` to weight `6`, outside the `\Delta \le 4` window.
2. Status: Proven. The first mixed witness `EZZ` is pushed to weight `7`, also outside the `\Delta \le 4` window.
3. Status: Proven. The exhaustive rank-6 audit in [`44-rank6-family-admission.md`](44-rank6-family-admission.md) finds no audited `EEZ`-type analogue at the first mixed order that could evade that threshold.
4. Status: Proven. Therefore no `Reven6+`-containing scalar survives inside the audited pairwise window at `\Delta \le 4`.

## Boundary

- Status: Proven. This lemma is only a pairwise audited-set composition result.
- Status: Proven. It does not prove triple, quadruple, quintuple, six-family, seven-family, or full-set sufficiency by itself.
- Status: Proven. It does not re-interpret the raw `22`-label rank-6 survivor bookkeeping as a corrected basis statement.
