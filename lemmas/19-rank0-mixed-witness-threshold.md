# Lemma 19: Rank-0 Mixed-Witness Threshold

## Setup

- Status: Proven. Consider a scalar-like rank-0 primitive family admitted into the free-fall MVP sector at fixed order `\Delta \le 4`.
- Status: Proven. The currently audited subclasses are:
  `R0a` unsuppressed bare scalar admission with primitive `S`,
  and `R0b` derivative-only or shift-symmetric admission with primitives `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S` but no bare `S`.
- Status: Proven. No hidden scalar-background rule such as `S = 0`, `S = const.`, or `\nabla_i S = 0` is assumed.

## Subclass `R0a`: Unsuppressed Bare Scalar

### First Self Witness

- Status: Proven. The operator `S` itself is a parity-even scalar of weight `1`.
- Status: Proven. Under the stated reduction rules, `S` is not removed by total derivatives, lower-order EOM, or any listed algebraic identity.
- Status: Proven. Therefore `W_{\mathrm{self}}(\mathrm{R0a}) = 1`.

### First Mixed Witness

- Status: Proven. The first mixed scalar on the electric baseline is `SE2`, of weight `3`.
- Status: Proven. If the magnetic family is also admitted, `SB2` ties at the same mixed weight and does not appear earlier.
- Status: Proven. No mixed witness can appear at weight `1` or `2`, because the baseline parity-even scalar sector has no weight-`0` scalar and no one-block contraction with `S` alone produces a mixed operator.
- Status: Proven. Therefore `W_{\mathrm{mix}}(\mathrm{R0a}) = 3`.

### Consequence

- Status: Proven. The subclass `R0a` satisfies

```math
W_{\mathrm{self}}(\mathrm{R0a}) = 1
< 3 = W_{\mathrm{mix}}(\mathrm{R0a}),
\qquad
W_{\min}(\mathrm{R0a}) = 1.
```

- Status: Proven. The current threshold `w_S \ge 5` or explicit exclusion of bare `S` is sharp within the audited bare-scalar class.
- Status: Proven. The self witness `S` is the true sharp canonical obstruction in this subclass.
- Status: Proven. The obstructed theorem layer is the promotion of the corrected `E/B` exact-current-set theorem candidate to a physically justified minimal-sector theorem.

## Subclass `R0b`: Derivative-Only Or Shift-Symmetric Scalar

### First Self Witness

- Status: Proven. The self operators `dotS2` and `gradS2` both survive at weight `4`.
- Status: Proven. The current audited canonical self representative is `dotS2`.
- Status: Proven. There is no surviving self witness below weight `4`, because `DtS` and `Dt2S` are total derivatives and `GradS` alone requires an acceleration contraction that vanishes by the lower-order EOM.
- Status: Proven. Therefore `W_{\mathrm{self}}(\mathrm{R0b}) = 4`.

### First Mixed Witness

- Status: Proven. On the electric baseline, the mixed survivor `DtS_E2` appears at weight `4`. (Corrected 2026-07-12: the former companion `divEGradS` contains the vanishing vacuum trace `\nabla_i E^{ij}` and is not generated under the STF-3 gradient kinematics.)
- Status: Proven. On the corrected `E/B` baseline, `DtS_B2` joins the same first mixed weight and does not appear earlier.
- Status: Proven. No mixed witness appears below weight `4`, because every weight-`2` or weight-`3` candidate either reduces by total derivatives, vanishes by the lower-order EOM, or fails to produce a parity-even scalar.
- Status: Proven. Therefore `W_{\mathrm{mix}}(\mathrm{R0b}) = 4`.

### Consequence

- Status: Proven. The subclass `R0b` satisfies

```math
W_{\mathrm{self}}(\mathrm{R0b}) = 4
= W_{\mathrm{mix}}(\mathrm{R0b}),
\qquad
W_{\min}(\mathrm{R0b}) = 4.
```

- Status: Proven. The current derivative-only threshold is sharp within the audited derivative-only scalar class.
- Status: Proven. The self witness `dotS2` remains a canonical representative, but it is not uniquely fundamental because mixed witnesses such as `DtS_E2` appear at the same effective order.
- Status: Proven. The obstructed theorem layer is the rescue of minimal-sector uniqueness after forbidding bare `S`.

## Boundary

- Status: Note. These conclusions are class-limited to the audited bare-scalar and derivative-only scalar subclasses.
- Status: Note. This lemma does not claim that scalar-family admission destroys finite-family fixed-order collapse.
