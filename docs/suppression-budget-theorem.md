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
  the self-only lower bound is `w_X \ge 3`, but the mixed-aware necessary threshold is `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness.
- Status: Proven. Unsuppressed bare scalar class:
  the witness is `S`, so uniqueness up to `\Delta \le 4` requires `w_S \ge 5`, or a symmetry/background rule that forbids bare `S`.
- Status: Proven. Derivative-only scalar class:
  the first audited witnesses are `\{DtS_B2,\ dotS2,\ DtS_E2,\ divEGradS,\ gradS2\}` at weight `4`, so uniqueness up to `\Delta \le 4` requires derivative-family blocks to be pushed above their current weight-`2` admission, for example to effective weight `3`, or an explicit rule removing those mixed derivative witnesses.
- Status: Proven. Genuine rank-1 vector class:
  the first audited witness is `V2`, so MVP-envelope sufficiency up to `\Delta \le 4` requires `w_V \ge 3`, or an explicit rule excluding or absorbing the primitive vector family.

## M9 Sharpness Status Within The Audited Classes

- Status: Proven. For the audited rank-2 STF class, the current threshold is mixed-aware: `W_{\mathrm{self}} = W_{\mathrm{mix}} = 2` in the unsuppressed audit, but the true `\Delta_{\max}=4` budget is governed by the mixed branch.
- Status: Proven. For the audited bare scalar class, the current threshold is self-only and sharp because `W_{\mathrm{self}} = 1 < 3 = W_{\mathrm{mix}}`.
- Status: Proven. For the audited derivative-only scalar class, the current threshold is tied-sharp because `W_{\mathrm{self}} = W_{\mathrm{mix}} = 4` in the current audit and both branches give `w_D \ge 3`.
- Status: Proven. For the audited rank-1 vector class, the current threshold is self-only and sharp because `W_{\mathrm{self}} = 2 < 3 = W_{\mathrm{mix}}`.
- Status: Proven. Therefore the current class-limited threshold map is no longer only about witness existence; it also records whether the active threshold is self-only, mixed-aware, or tied-sharp.

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

## Audited-Set Joint Sufficiency At `\Delta_{\max} = 4`

- Status: Proven. The current pairwise and triple composition audits for the audited family classes `R2`, `R0a`, and `R0b` find no new surviving operator beyond the baseline electric sector once the current thresholds are imposed together.
- Status: Proven. Therefore the current thresholds are jointly sufficient for the currently audited family classes at `\Delta \le 4`.
- Status: Proven. This audited-set joint sufficiency does not upgrade the theorem into a universal composition statement for arbitrary family catalogs.

## Boundary

- Status: Proven. This theorem is about minimal-sector uniqueness only.
- Status: Proven. Failure to meet the threshold reinforces the negative uniqueness branch.
- Status: Proven. Meeting the threshold does not by itself prove the stronger theorem, because additional witnesses may still survive in some subclasses.
- Status: Proven. The positive finite-family collapse branch remains separate and is not falsified merely by the appearance of new witnesses below threshold.
