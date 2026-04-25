# Lemma 39: Rank-5 Threshold Formula

## Scope

- Status: Proven. This note classifies the witness threshold for the genuine local parity-even fully symmetric trace-free rank-5 family `U_{ijklm}` audited in [`38-rank5-family-admission.md`](38-rank5-family-admission.md).
- Status: Proven. The formula is class-limited to the audited `Rodd5+` family class.

## Formulae

- Status: Proven. `W_{\mathrm{self}}(Rodd5+; w_U) = 2 w_U`.
- Status: Proven. `W_{\mathrm{mix}}(Rodd5+; w_U) = 2 w_U + 1`.
- Status: Proven. Therefore `W_{\min}(Rodd5+; w_U) = 2 w_U`.

## `\Delta_{\max} = 4` Consequence

- Status: Proven. At fixed order `\Delta_{\max} = 4`, minimal-sector uniqueness against the audited rank-5 family class requires `2 w_U > 4`.
- Status: Proven. Equivalently, the current class-limited necessary threshold is `w_U \ge 3`, or an explicit rule excluding or absorbing the primitive rank-5 family.

## Threshold Type

- Status: Proven. The audited `Rodd5+` threshold is self-only sharp.
- Status: Proven. In the unsuppressed case `w_U = 1`, one has `W_{\mathrm{self}} = 2 < 3 = W_{\mathrm{mix}}`.

## Boundary

- Status: Proven. This threshold is necessary for the audited rank-5 family class only; it is not sufficient for harmless enlarged-set composition.
- Status: Proven. Meeting `w_U \ge 3` does not by itself prove enlarged-set harmlessness; the separate composition closure for the enlarged audited set is recorded in [`40-r5-pairwise-composition-audit.md`](40-r5-pairwise-composition-audit.md) and [`41-r5-augmented-composition-audit.md`](41-r5-augmented-composition-audit.md).
