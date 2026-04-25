# Lemma 35: Rank-4 Threshold Formula

## Scope

- Status: Proven. This lemma derives the `\Delta_{\max} = 4` witness-threshold formulas for the audited genuine rank-4 family class represented by a fully symmetric trace-free primitive family `Q_{ijkl}`.
- Status: Proven. The baseline for the mixed witness is the current enlarged audited-set baseline, which reduces back to the electric sector once the already audited thresholds are imposed.

## Formulae

- Status: Proven. If the primitive rank-4 family is assigned effective block weight `w_Q`, then the first self witness is `Q2`, so

```math
W_{\mathrm{self}}(Reven4+; w_Q) = 2 w_Q.
```

- Status: Proven. The first mixed witness is `EQQ`, so

```math
W_{\mathrm{mix}}(Reven4+; w_Q) = 2 w_Q + 1.
```

- Status: Proven. Therefore

```math
W_{\min}(Reven4+; w_Q) = \min(2 w_Q,\ 2 w_Q + 1) = 2 w_Q.
```

## `\Delta_{\max} = 4` Consequence

- Status: Proven. At unsuppressed weight `w_Q = 1`, the current case is
  `(W_{\mathrm{self}}, W_{\mathrm{mix}}, W_{\min}) = (2, 3, 2)`.
- Status: Proven. To keep the audited rank-4 witnesses outside `\Delta \le 4`, minimal-sector uniqueness requires
  `2 w_Q > 4`, hence `w_Q \ge 3`.
- Status: Proven. The corresponding class-limited necessary threshold statement is:
  `w_Q \ge 3`, or an explicit rule excluding or absorbing the primitive rank-4 family.

## Threshold Type

- Status: Proven. The rank-4 family threshold is self-only, because `W_{\mathrm{self}} < W_{\mathrm{mix}}` in the audited unsuppressed class.
- Status: Proven. It is also sharp within the audited class, because the advertised threshold is governed exactly by `W_{\min}`.

## Boundary

- Status: Proven. This threshold theorem is class-limited to the audited genuine STF rank-4 family.
- Status: Proven. Meeting the threshold does not by itself prove enlarged-set harmlessness; the separate composition closure for the enlarged audited set is recorded in [`36-r4-pairwise-composition-audit.md`](36-r4-pairwise-composition-audit.md) and [`37-r4-augmented-composition-audit.md`](37-r4-augmented-composition-audit.md).
