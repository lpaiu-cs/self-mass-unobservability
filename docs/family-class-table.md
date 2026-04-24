# Family-Class Table

## Purpose

- Status: Proven. This table compresses the current audited family classes into a single witness-threshold view.
- Status: Proven. The thresholds listed below are necessary budgets for minimal-sector uniqueness up to `\Delta_{\max} = 4`.
- Status: Proven. The self-versus-mixed split below is class-limited to the currently audited family classes.

| Family class | Class status | First self witness | Self weight | First mixed witness | Mixed weight | `W_min` | Sharpness status | Required suppression budget | Current confidence level | Theorem layer obstructed | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Rank-2 STF family | Proven | `X2` | `2` | `EX` | `2` | `2` | Sharp within audited class; self and mixed tie | Effective family weight `w_X \ge 3`, or explicit exclusion / background restriction | Proven | Promotion of electric-only exact-current-set theorem to a physically justified minimal-sector theorem | [`../lemmas/15-rank2-family-admission.md`](../lemmas/15-rank2-family-admission.md), [`../lemmas/18-rank2-mixed-witness-threshold.md`](../lemmas/18-rank2-mixed-witness-threshold.md) |
| Rank-0 scalar family with bare source | Proven | `S` | `1` | `SE2` | `3` | `1` | Sharp within audited class; self dominates mixed | Bare-scalar weight `w_S \ge 5`, or explicit rule forbidding bare `S` | Proven | Promotion of corrected `E/B` exact-current-set theorem to a physically justified minimal-sector theorem | [`../lemmas/16-rank0-family-admission.md`](../lemmas/16-rank0-family-admission.md), [`../lemmas/19-rank0-mixed-witness-threshold.md`](../lemmas/19-rank0-mixed-witness-threshold.md) |
| Rank-0 derivative-only scalar family | Proven | `dotS2` | `4` | `DtS_E2` | `4` | `4` | Sharp within audited class; self and mixed tie | Raise derivative-family block weight above `2`, for example to `3`, or explicitly remove the mixed derivative witnesses | Proven | Rescue of minimal-sector uniqueness after removing bare `S` | [`../lemmas/14-derivative-only-scalar-audit.md`](../lemmas/14-derivative-only-scalar-audit.md), [`../lemmas/19-rank0-mixed-witness-threshold.md`](../lemmas/19-rank0-mixed-witness-threshold.md) |

## Reading Rule

- Status: Proven. The table records necessary thresholds only.
- Status: Proven. Crossing the listed threshold does not automatically prove uniqueness; it only removes the current lowest audited witness from the `\Delta \le 4` window.
- Status: Proven. The new self/mixed columns explain whether the lowest audited witness is self-dominated or tied with a mixed witness.
