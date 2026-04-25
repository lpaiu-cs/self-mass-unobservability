# Lemma 31: Rank-3 Threshold Formula

## Scope

- Status: Proven. This lemma derives the `\Delta_{\max} = 4` witness-threshold formulas for the audited genuine rank-3 family class represented by a fully symmetric trace-free primitive family `T_{ijk}`.
- Status: Proven. The baseline for the mixed witness is the current enlarged audited-set baseline, which reduces back to the electric sector once the already audited thresholds are imposed.

## Formulae

- Status: Proven. If the primitive rank-3 family is assigned effective block weight `w_T`, then the first self witness is `T2`, so

```math
W_{\mathrm{self}}(Rodd+; w_T) = 2 w_T.
```

- Status: Proven. The first mixed witness is `ETT`, so

```math
W_{\mathrm{mix}}(Rodd+; w_T) = 2 w_T + 1.
```

- Status: Proven. Therefore

```math
W_{\min}(Rodd+; w_T) = \min(2 w_T,\ 2 w_T + 1) = 2 w_T.
```

## `\Delta_{\max} = 4` Consequence

- Status: Proven. At unsuppressed weight `w_T = 1`, the current case is
  `(W_{\mathrm{self}}, W_{\mathrm{mix}}, W_{\min}) = (2, 3, 2)`.
- Status: Proven. To keep the audited rank-3 witnesses outside `\Delta \le 4`, minimal-sector uniqueness requires
  `2 w_T > 4`, hence `w_T \ge 3`.
- Status: Proven. The corresponding class-limited necessary threshold statement is:
  `w_T \ge 3`, or an explicit rule excluding or absorbing the primitive rank-3 family.

## Threshold Type

- Status: Proven. The rank-3 family threshold is self-only, because `W_{\mathrm{self}} < W_{\mathrm{mix}}` in the audited unsuppressed class.
- Status: Proven. It is also sharp within the audited class, because the advertised threshold is governed exactly by `W_{\min}`.

## Boundary

- Status: Proven. This threshold theorem is class-limited to the audited genuine STF rank-3 family.
- Status: Proven. Meeting the threshold does not yet prove harmless composition; a new enlarged audited-set composition re-close is required next.
