# Sharp Threshold Status

## Purpose

- Status: Proven. This note records whether the current suppression budget for each audited family class is already sharp, merely a lower bound, or still unresolved.
- Status: Proven. The classification below is class-limited to the currently audited family classes at `\Delta_{\max} = 4`.

| Family class | First self witness | Self weight | First mixed witness | Mixed weight | `W_min` | Threshold class | Current consistency statement |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Rank-2 STF family | `X2` | `2` | `EX` | `2` | `2` | Mixed-aware | Self-only lower bound `w_X \ge 3`; mixed-aware necessary threshold `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness. |
| Rank-0 scalar family with bare source | `S` | `1` | `SE2` | `3` | `1` | Self-only | `w_S \ge 5`, or explicit exclusion of bare `S`. |
| Rank-0 derivative-only scalar family | `dotS2` | `4` | `DtS_E2` | `4` | `4` | Tied-sharp | `w_D \ge 3`; the self and mixed branches yield the same sharp budget at `\Delta_{\max}=4`. |
| Rank-1 vector family | `V2` | `2` | `EVV` | `3` | `2` | Self-only | `w_V \ge 3`, or explicit exclusion or absorption of the primitive vector family. |
| Rank-3 STF family | `T2` | `2` | `ETT` | `3` | `2` | Self-only | `w_T \ge 3`, or explicit exclusion or absorption of the primitive rank-3 family. |
| Rank-4 STF family | `Q2` | `2` | `EQQ` | `3` | `2` | Self-only | `w_Q \ge 3`, or explicit exclusion or absorption of the primitive rank-4 family. |
| Rank-5 STF family | `U2` | `2` | `EUU` | `3` | `2` | Self-only | `w_U \ge 3`, or explicit exclusion or absorption of the primitive rank-5 family. |

## Current Readout

- Status: Proven. No audited family class is currently unresolved.
- Status: Proven. The rank-2 STF class is no longer described by a self-only threshold; its consistent budget is mixed-aware.
- Status: Proven. The bare scalar class remains self-only.
- Status: Proven. The derivative-only scalar class is tied-sharp at the current fixed-order threshold.
- Status: Proven. The genuine vector class is self-only at the current fixed-order threshold.
- Status: Proven. The genuine rank-3 STF class is self-only at the current fixed-order threshold.
- Status: Proven. The genuine rank-4 STF class is self-only at the current fixed-order threshold.
- Status: Proven. The genuine rank-5 STF class is self-only at the current fixed-order threshold.
- Status: Proven. No audited family class currently exhibits a mixed witness below its first self witness at `\Delta \le 4`.
- Status: Proven. The sharp-threshold classification bottleneck is closed for the currently audited family classes.
- Status: Proven. The repo-level live bottleneck has moved from threshold sharpness to the omitted rank-4 contraction `EEQ` and the resulting high-rank exhaustiveness patch before any move to `Reven6+`.
- Status: Conjectural. What remains open is not witness sharpness inside the currently audited classes, but family-envelope completeness beyond the currently audited and composition-closed family set.
