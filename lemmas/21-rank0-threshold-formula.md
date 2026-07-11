# Lemma 21: Rank-0 Threshold Formula

## Setup

- Status: Proven. Consider the audited scalar-like rank-0 family classes in the free-fall MVP sector at fixed order `\Delta \le 4`.
- Status: Proven. The two audited subclasses remain:
  `R0a` bare scalar admission with primitive `S`,
  and `R0b` derivative-only admission with primitives `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S` but no bare `S`.

## Subclass `R0a`: Bare Scalar

- Status: Proven. Let `w_S` denote the counting weight assigned to the primitive `S`.
- Status: Proven. The first self witness is

```math
S,
\qquad
W_{\mathrm{self}}(\mathrm{R0a}; w_S) = w_S.
```

- Status: Proven. The first mixed witness on the audited baseline is `SE2`, so

```math
W_{\mathrm{mix}}(\mathrm{R0a}; w_S) = w_S + 2.
```

- Status: Proven. Therefore

```math
W_{\min}(\mathrm{R0a}; w_S) = w_S.
```

- Status: Proven. At `\Delta_{\max}=4`, the current threshold is self-only and sharp:
  `w_S > 4`, hence `w_S \ge 5`, or explicit exclusion of bare `S`.

## Subclass `R0b`: Derivative-Only Scalar

- Status: Proven. Let `w_D := \Delta[D_\tau S] = \Delta[\nabla_i S]`.
- Status: Proven. In the audited derivative-only counting, `\Delta[D_\tau^2 S] = w_D + 1`.
- Status: Proven. The first self witnesses are `dotS2` and `gradS2`, so

```math
W_{\mathrm{self}}(\mathrm{R0b}; w_D) = 2 w_D.
```

- Status: Proven. The first mixed witnesses are `DtS_E2`, `divEGradS`, and, on the corrected `E/B` baseline, also `DtS_B2`, so

```math
W_{\mathrm{mix}}(\mathrm{R0b}; w_D) = w_D + 2.
```

- Status: Proven. Therefore

```math
W_{\min}(\mathrm{R0b}; w_D)
=
\min(2 w_D,\ w_D + 2).
```

- Status: Proven. In the current audited case `w_D = 2`,

```math
W_{\mathrm{self}} = W_{\mathrm{mix}} = W_{\min} = 4.
```

- Status: Proven. At `\Delta_{\max}=4`, both branches give the same sharp budget:
  `2 w_D > 4` and `w_D + 2 > 4` both reduce to `w_D \ge 3`.
- Status: Proven. Therefore the current derivative-only threshold is tied-sharp rather than purely self-only or strictly mixed-dominated.

## Boundary

- Status: Note. These formulas are class-limited to the audited scalar subclasses only.
- Status: Proven. The lemma does not imply sufficiency of the thresholds for uniqueness.
- Status: Proven. The positive finite-family fixed-order collapse target remains separate.
