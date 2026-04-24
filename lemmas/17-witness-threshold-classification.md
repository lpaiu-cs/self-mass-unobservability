# Lemma 17: Witness-Threshold Classification

## Statement

- Status: Proven. For the currently audited family classes, the smallest explicit witness determines a necessary suppression threshold for minimal-sector uniqueness up to fixed order `\Delta_{\max}`.
- Status: Proven. If the class witness remains at weight `\le \Delta_{\max}`, then minimal-sector uniqueness fails against that class.
- Status: Proven. Therefore any attempted uniqueness salvage must at least push that witness above `\Delta_{\max}` or remove it by an explicit symmetry, ordering, or background rule.

## Class-Limited Thresholds

- Status: Proven. Rank-2 STF class:
  `W_{\mathrm{wit}} = 2 w_X` from the quadratic witness `X2`.
- Status: Proven. Hence uniqueness up to `\Delta_{\max}` requires `2 w_X > \Delta_{\max}`.
- Status: Proven. At `\Delta_{\max} = 4`, this gives the necessary threshold `w_X \ge 3`.

- Status: Proven. Bare rank-0 scalar class:
  `W_{\mathrm{wit}} = w_S` from the linear witness `S`.
- Status: Proven. Hence uniqueness up to `\Delta_{\max}` requires `w_S > \Delta_{\max}`.
- Status: Proven. At `\Delta_{\max} = 4`, this gives the necessary threshold `w_S \ge 5`, unless bare `S` is forbidden outright.

- Status: Proven. Derivative-only rank-0 scalar class:
  the current lowest audited witnesses are mixed one-family terms such as `DtS_E2` and `divEGradS`, together with pure derivative terms such as `dotS2` and `gradS2`.
- Status: Proven. Their current minimal audited weight is `4`.
- Status: Proven. Therefore uniqueness up to `\Delta_{\max} = 4` requires an additional derivative-family suppression beyond the present audit, for example raising derivative-family insertions from weight `2` to weight `3`.

## Interpretation

- Status: Proven. The threshold map explains why `B2`, `S`, and `dotS2` belong to one structural pattern.
- Status: Proven. They differ in witness weight and therefore in the size of the suppression budget required to evade them.
- Status: Proven. This is a classification lemma for the audited classes, not a universal theorem for all primitive families.

## Boundary

- Status: Proven. These thresholds are necessary conditions for uniqueness, not guaranteed sufficient conditions.
- Status: Proven. The positive finite-family collapse branch remains alive even when uniqueness thresholds are not met.
