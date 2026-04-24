# Mixed-Witness Thresholds

## Definitions

- Status: Proven. Fix a baseline admitted primitive catalog `\mathcal C_0` and adjoin one new primitive family class `\mathcal F`.
- Status: Proven. A self witness for `\mathcal F` is a surviving normal-form scalar operator built only from primitives belonging to `\mathcal F`.
- Status: Proven. A mixed witness for `\mathcal F` is a surviving normal-form scalar operator containing at least one primitive from `\mathcal F` and at least one primitive from the baseline catalog `\mathcal C_0`.
- Status: Proven. Define

```math
W_{\mathrm{self}}(\mathcal F)
\ :=\
\min \{ \Delta[\mathcal O] : \mathcal O \text{ is a self witness for } \mathcal F \},
```

```math
W_{\mathrm{mix}}(\mathcal F)
\ :=\
\min \{ \Delta[\mathcal O] : \mathcal O \text{ is a mixed witness for } \mathcal F \},
```

```math
W_{\min}(\mathcal F) := \min\!\bigl(W_{\mathrm{self}}(\mathcal F),\ W_{\mathrm{mix}}(\mathcal F)\bigr).
```

- Status: Proven. The M8 suppression-budget theorem only required one audited witness weight `W_{\mathrm{wit}}(\mathcal F)` satisfying `W_{\mathrm{wit}}(\mathcal F) \le \Delta_{\max}`.
- Status: Proven. The M9 sharpness question is whether the advertised class threshold is controlled by the first self witness, the first mixed witness, or both inside the currently audited family classes.

## Sharpness Language

- Status: Proven. A threshold is sharp within an audited family class if the current advertised obstruction weight equals `W_{\min}(\mathcal F)`.
- Status: Proven. A threshold is self-dominated when `W_{\mathrm{self}}(\mathcal F) < W_{\mathrm{mix}}(\mathcal F)`.
- Status: Proven. A threshold is mixed-dominated when `W_{\mathrm{mix}}(\mathcal F) < W_{\mathrm{self}}(\mathcal F)`.
- Status: Proven. A threshold is tied when `W_{\mathrm{self}}(\mathcal F) = W_{\mathrm{mix}}(\mathcal F)`.
- Status: Proven. A necessary threshold from M8 is only a lower bound if the currently advertised witness sits above `W_{\min}(\mathcal F)`.
- Status: Proven. A family class remains unresolved only if the currently audited notes do not yet determine the order relation between `W_{\mathrm{self}}` and `W_{\mathrm{mix}}`.

## M9 Question

- Status: Proven. The current M8 theorem is a necessary-threshold statement, not a sufficient-threshold statement.
- Status: Proven. M9 does not ask whether witnesses exist; that question is already settled for the audited family classes.
- Status: Proven. M9 asks whether the currently advertised thresholds are sharp within those audited classes once self and mixed witnesses are separated.

## Current Audited Outcome At `\Delta_{\max} = 4`

- Status: Proven. For the audited STF rank-2 family class, `W_{\mathrm{self}} = 2` and `W_{\mathrm{mix}} = 2`, so the current threshold is sharp within the audited class but not purely self-witness-driven.
- Status: Proven. For the audited bare scalar family class, `W_{\mathrm{self}} = 1` and `W_{\mathrm{mix}} = 3`, so the current threshold is sharp and self-dominated within the audited class.
- Status: Proven. For the audited derivative-only scalar family class, `W_{\mathrm{self}} = 4` and `W_{\mathrm{mix}} = 4`, so the current threshold is sharp within the audited class but not purely self-witness-driven.
- Status: Proven. No audited family class currently exhibits `W_{\mathrm{mix}} < W_{\mathrm{self}}` at `\Delta \le 4`.
- Status: Conjectural. The remaining negative-branch burden is to determine whether the same sharp class-limited pattern persists beyond the currently audited families, not to re-establish witness existence.
