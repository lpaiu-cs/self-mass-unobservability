# Threshold Formulas

## Purpose

- Status: Note. This note lifts the audited witness-threshold results from case labels into formula-level statements for the currently audited family classes.
- Status: Proven. The formulas are class-limited to the audited families only; they are not universal statements for every imaginable primitive family.

## Notation

- Status: Proven. Let `w_{\mathcal F}` denote the effective counting weight of the primitive family block admitted on top of the chosen baseline catalog.
- Status: Proven. Define

```math
W_{\mathrm{self}}(\mathcal F; w_{\mathcal F})
```

as the smallest surviving self-witness weight for family class `\mathcal F`,

```math
W_{\mathrm{mix}}(\mathcal F; w_{\mathcal F})
```

as the smallest surviving mixed-witness weight for the same family class against the relevant baseline catalog, and

```math
W_{\min}(\mathcal F; w_{\mathcal F})
:=
\min\!\bigl(
W_{\mathrm{self}}(\mathcal F; w_{\mathcal F}),
W_{\mathrm{mix}}(\mathcal F; w_{\mathcal F})
\bigr).
```

- Status: Proven. Minimal-sector uniqueness up to `\Delta_{\max}` requires

```math
W_{\min}(\mathcal F; w_{\mathcal F}) > \Delta_{\max}
```

for every admitted family class under the currently allowed rules.

## Threshold-Type Language

- Status: Proven. A class is self-only when the current necessary threshold is governed by `W_{\mathrm{self}}` and the mixed branch is weaker.
- Status: Proven. A class is mixed-aware when the current necessary threshold is stricter than the self-only lower bound because `W_{\mathrm{mix}}` controls `W_{\min}`.
- Status: Proven. A class is tied-sharp when the self and mixed branches yield the same sharp threshold at the fixed audited `\Delta_{\max}` even though the witnesses are distinct.

## Audited Family Formulas

| Family class | `W_{\mathrm{self}}` | `W_{\mathrm{mix}}` | `W_{\min}` | Threshold type at `\Delta_{\max}=4` |
| --- | --- | --- | --- | --- |
| Rank-2 STF | `2 w_X` | `w_X + 1` | `\min(2 w_X,\ w_X + 1)` | Mixed-aware |
| Rank-0 bare scalar | `w_S` | `w_S + 2` | `w_S` | Self-only |
| Rank-0 derivative-only scalar | `2 w_D` | `w_D + 2` | `\min(2 w_D,\ w_D + 2)` | Tied-sharp |
| Rank-1 vector | `2 w_V` | `2 w_V + 1` | `2 w_V` | Self-only |
| Rank-3 STF | `2 w_T` | `2 w_T + 1` | `2 w_T` | Self-only |
| Rank-4 STF | `2 w_Q` | `2 w_Q + 1` | `2 w_Q` | Self-only |
| Rank-5 STF | `2 w_U` | `2 w_U + 1` | `2 w_U` | Self-only |

## Immediate Consequences At `\Delta_{\max}=4`

- Status: Proven. Rank-2 STF class:
  the self-only lower bound is `w_X \ge 3`, but the mixed-aware necessary threshold is `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness.
- Status: Proven. Rank-0 bare scalar class:
  `W_{\min} = w_S`, so the current threshold is self-only and sharp: `w_S \ge 5`, or explicit exclusion of bare `S`.
- Status: Proven. Rank-0 derivative-only scalar class:
  `W_{\mathrm{self}} = 2 w_D` and `W_{\mathrm{mix}} = w_D + 2` give the same sharp budget at `\Delta_{\max}=4`, namely `w_D \ge 3`.
- Status: Proven. Rank-1 vector class:
  `W_{\min} = 2 w_V`, so the current threshold is self-only and sharp: `w_V \ge 3`, or explicit exclusion or absorption of the primitive vector family.
- Status: Proven. Rank-3 STF class:
  `W_{\min} = 2 w_T`, so the current threshold is self-only and sharp: `w_T \ge 3`, or explicit exclusion or absorption of the primitive rank-3 family.
- Status: Proven. Rank-4 STF class:
  `W_{\min} = 2 w_Q`, so the current threshold is self-only and sharp: `w_Q \ge 3`, or explicit exclusion or absorption of the primitive rank-4 family.
- Status: Proven. Rank-5 STF class:
  `W_{\min} = 2 w_U`, so the current threshold is self-only and sharp: `w_U \ge 3`, or explicit exclusion or absorption of the primitive rank-5 family.

## Boundary

- Status: Proven. These formulas classify the audited family classes only.
- Status: Proven. They do not by themselves prove uniqueness once the thresholds are met.
- Status: Proven. The positive theorem target remains finite-family fixed-order collapse.
