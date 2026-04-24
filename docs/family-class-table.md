# Family-Class Table

## Purpose

- Status: Proven. This table compresses the current audited family classes into a single witness-threshold view.
- Status: Proven. The thresholds listed below are necessary budgets for minimal-sector uniqueness up to `\Delta_{\max} = 4`.

| Family class | Status | Canonical audited witness | Current witness weight | Necessary budget for uniqueness up to `\Delta \le 4` | If unmet | Evidence |
| --- | --- | --- | --- | --- | --- | --- |
| Rank-2 STF family | Proven | `X2` | `2` | Effective family weight `w_X \ge 3`, or explicit exclusion / background restriction | Minimal-sector uniqueness fails | [`../lemmas/15-rank2-family-admission.md`](../lemmas/15-rank2-family-admission.md) |
| Rank-0 scalar family with bare source | Proven | `S` | `1` | Bare-scalar weight `w_S \ge 5`, or explicit rule forbidding bare `S` | Minimal-sector uniqueness fails | [`../lemmas/16-rank0-family-admission.md`](../lemmas/16-rank0-family-admission.md) |
| Rank-0 derivative-only scalar family | Proven | `dotS2` | `4` | Raise derivative-family block weight above `2`, for example to `3`, or explicitly remove the mixed derivative witnesses | Minimal-sector uniqueness fails | [`../lemmas/14-derivative-only-scalar-audit.md`](../lemmas/14-derivative-only-scalar-audit.md) |

## Reading Rule

- Status: Proven. The table records necessary thresholds only.
- Status: Proven. Crossing the listed threshold does not automatically prove uniqueness; it only removes the current lowest audited witness from the `\Delta \le 4` window.
