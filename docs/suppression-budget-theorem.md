# Suppression-Budget Theorem

## Scope

- Status: Proven. This note is a class-limited theorem for the currently audited family classes only.
- Status: Proven. It does not claim a universal statement for every imaginable primitive family.
- Status: Proven. The theorem gives necessary suppression budgets for minimal-sector uniqueness up to fixed order `\Delta_{\max}`, not automatic sufficient conditions for every subclass.

## Sharp Necessary Condition

- Status: Proven. Let `\mathcal F` be one of the audited family classes and let `W_{\mathrm{wit}}(\mathcal F)` denote the smallest audited witness weight for that class under the currently allowed reduction rules.
- Status: Proven. If `W_{\mathrm{wit}}(\mathcal F) \le \Delta_{\max}`, then admitting `\mathcal F` unsuppressed obstructs minimal-sector uniqueness up to `\Delta_{\max}`.
- Status: Proven. Therefore minimal-sector uniqueness up to `\Delta_{\max}` requires an explicit suppression budget, ordering rule, symmetry rule, or background restriction strong enough to push the audited witness weight above `\Delta_{\max}` or remove that witness entirely.

## Audited Thresholds At `\Delta_{\max} = 4`

- Status: Proven. Rank-2 STF family class:
  the quadratic witness is `X2`, so uniqueness up to `\Delta \le 4` requires effective family weight `w_X \ge 3`, or an explicit rule that excludes the family on the allowed backgrounds.
- Status: Proven. Unsuppressed bare scalar class:
  the witness is `S`, so uniqueness up to `\Delta \le 4` requires `w_S \ge 5`, or a symmetry/background rule that forbids bare `S`.
- Status: Proven. Derivative-only scalar class:
  the first audited witnesses are `\{DtS_B2,\ dotS2,\ DtS_E2,\ divEGradS,\ gradS2\}` at weight `4`, so uniqueness up to `\Delta \le 4` requires derivative-family blocks to be pushed above their current weight-`2` admission, for example to effective weight `3`, or an explicit rule removing those mixed derivative witnesses.

## M9 Sharpness Status Within The Audited Classes

- Status: Proven. For the audited rank-2 STF class, the current threshold is sharp within the class because `W_{\mathrm{self}} = W_{\mathrm{mix}} = 2`.
- Status: Proven. For the audited bare scalar class, the current threshold is sharp within the class because `W_{\mathrm{self}} = 1 < 3 = W_{\mathrm{mix}}`.
- Status: Proven. For the audited derivative-only scalar class, the current threshold is sharp within the class because `W_{\mathrm{self}} = W_{\mathrm{mix}} = 4`.
- Status: Proven. Therefore the current class-limited threshold map is no longer only about witness existence; it also records whether the first obstruction is self-dominated or tied with a mixed witness.

## What This Theorem Buys

- Status: Proven. The repeated family-attack pattern is now compressed into a witness-threshold statement.
- Status: Proven. The theorem explains why `B2`, `S`, and `dotS2` belong to the same structural pattern rather than being isolated curiosities.
- Status: Proven. The theorem does not rescue minimal-sector uniqueness.
- Status: Proven. It only tells us the minimum budget any attempted rescue must at least pay for the audited family classes.

## What Still Remains For Sharpness

- Status: Proven. The theorem above remains a necessary-threshold theorem even after the M9 sharpness split.
- Status: Proven. Sharpness within the audited classes does not make the listed budgets sufficient for uniqueness.
- Status: Conjectural. It remains open whether broader subclasses or future family audits could produce a smaller mixed witness than the currently audited class-limited `W_{\min}` values.
- Status: Conjectural. It remains open whether there is an abstract theorem computing `W_{\min}(\mathcal F)` beyond the currently audited family classes.

## Boundary

- Status: Proven. This theorem is about minimal-sector uniqueness only.
- Status: Proven. Failure to meet the threshold reinforces the negative uniqueness branch.
- Status: Proven. Meeting the threshold does not by itself prove the stronger theorem, because additional witnesses may still survive in some subclasses.
- Status: Proven. The positive finite-family collapse branch remains separate and is not falsified merely by the appearance of new witnesses below threshold.
