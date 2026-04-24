# Lemma 17: Witness-Threshold Classification

## Statement

- Status: Proven. For the currently audited family classes, the smallest explicit witness determines a necessary suppression threshold for minimal-sector uniqueness up to fixed order `\Delta_{\max}`.
- Status: Proven. If the class witness remains at weight `\le \Delta_{\max}`, then minimal-sector uniqueness fails against that class.
- Status: Proven. Therefore any attempted uniqueness salvage must at least push that witness above `\Delta_{\max}` or remove it by an explicit symmetry, ordering, or background rule.

## Class-Limited Thresholds

- Status: Proven. Rank-2 STF class:
  the self witness gives `W_{\mathrm{self}} = 2 w_X`, while the mixed witness gives `W_{\mathrm{mix}} = w_X + 1`.
- Status: Proven. Therefore `W_{\min} = \min(2 w_X,\ w_X + 1)`.
- Status: Proven. The self-only lower bound is `2 w_X > \Delta_{\max}`.
- Status: Proven. At `\Delta_{\max} = 4`, that self-only lower bound is `w_X \ge 3`.
- Status: Proven. The current mixed-aware necessary threshold is stricter:
  `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness.

- Status: Proven. Bare rank-0 scalar class:
  `W_{\mathrm{self}} = w_S` from the linear witness `S`, while `W_{\mathrm{mix}} = w_S + 2`.
- Status: Proven. Therefore `W_{\min} = w_S`.
- Status: Proven. Hence uniqueness up to `\Delta_{\max}` requires `w_S > \Delta_{\max}`.
- Status: Proven. At `\Delta_{\max} = 4`, this gives the self-only sharp threshold `w_S \ge 5`, unless bare `S` is forbidden outright.

- Status: Proven. Derivative-only rank-0 scalar class:
  the first self branch is `W_{\mathrm{self}} = 2 w_D`, while the first mixed branch is `W_{\mathrm{mix}} = w_D + 2`.
- Status: Proven. Therefore `W_{\min} = \min(2 w_D,\ w_D + 2)`.
- Status: Proven. At `\Delta_{\max} = 4`, both branches give the same sharp threshold `w_D \ge 3`.
- Status: Proven. Therefore the derivative-only class is tied-sharp at the current fixed-order threshold.

## Interpretation

- Status: Proven. The threshold map explains why `B2`, `S`, and `dotS2` belong to one structural pattern.
- Status: Proven. They differ in whether the current threshold is self-only, mixed-aware, or tied-sharp.
- Status: Proven. For the rank-2 class, the old `w_X \ge 3` line is only a self-only lower bound, not the current mixed-aware necessary threshold.
- Status: Proven. This is a classification lemma for the audited classes, not a universal theorem for all primitive families.

## Boundary

- Status: Proven. These thresholds are necessary conditions for uniqueness, not guaranteed sufficient conditions.
- Status: Proven. The positive finite-family collapse branch remains alive even when uniqueness thresholds are not met.
