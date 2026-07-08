# Lemma 18: Rank-2 Mixed-Witness Threshold

## Setup

- Status: Proven. Let `X_ij` be an unsuppressed parity-even STF rank-2 primitive family admitted at weight `1` on top of the electric baseline catalog.
- Status: Proven. The baseline catalog already contains the electric tidal STF block `E_ij`.
- Status: Proven. No cross-family orthogonality rule such as `Tr(EX) = 0` is assumed.
- Status: Proven. No additional transversality, vacuum/Bianchi identity, or background restriction is assumed beyond [`../docs/assumptions-ledger.md`](../docs/assumptions-ledger.md) and [`../docs/reduction-rules.md`](../docs/reduction-rules.md).

## First Self Witness

- Status: Proven. The quadratic invariant

```math
X2 := Tr(X^2)
```

is a parity-even scalar of weight `2`.
- Status: Proven. `X2` is not a total derivative because it contains no `D_\tau`.
- Status: Proven. `X2` is not removed by the lower-order free-fall EOM because it carries no acceleration factor.
- Status: Proven. `X2` is not removed by the currently allowed algebraic identities, which do not collapse a quadratic norm.
- Status: Proven. Therefore `W_{\mathrm{self}}(\mathrm{R2}) \le 2`.

## First Mixed Witness

- Status: Proven. The mixed quadratic invariant

```math
EX := Tr(EX)
```

is also a parity-even scalar of weight `2`.
- Status: Proven. `EX` is not a total derivative because it contains no `D_\tau`.
- Status: Proven. `EX` is not removed by the lower-order free-fall EOM because it carries no acceleration factor.
- Status: Proven. `EX` is not removed by the currently stated algebraic identities, because no cross-family quadratic identity has been assumed.
- Status: Proven. Therefore `W_{\mathrm{mix}}(\mathrm{R2}) \le 2`.

## Minimality

- Status: Proven. No weight-`1` parity-even scalar can be built from `X_ij`, because a single STF rank-2 tensor cannot produce a scalar without another rank-2 object.
- Status: Proven. Therefore no self or mixed witness appears below weight `2` in the audited rank-2 class.
- Status: Proven. Hence

```math
W_{\mathrm{self}}(\mathrm{R2}) = 2,\qquad
W_{\mathrm{mix}}(\mathrm{R2}) = 2,\qquad
W_{\min}(\mathrm{R2}) = 2.
```

## Consequence

- Status: Proven. The current unsuppressed audit is tied at weight `2`, but the threshold consequence at `\Delta_{\max} = 4` is mixed-aware rather than self-only.
- Status: Proven. In particular, the old self-only line `w_X \ge 3` is only a lower bound; the mixed-aware necessary threshold is derived separately in [`20-rank2-threshold-formula.md`](20-rank2-threshold-formula.md).
- Status: Proven. The self witness `X2` remains a canonical obstruction, but it is not uniquely fundamental because the mixed witness `EX` appears at the same effective order.
- Status: Proven. The obstructed theorem layer is the promotion of the electric-only exact-current-set theorem candidate to a physically justified minimal-sector theorem.

## Boundary

- Status: Note. This lemma is class-limited to the audited parity-even STF rank-2 family and does not claim a universal statement for every rank-2 tensor family.
- Status: Note. This lemma does not say that finite-family fixed-order collapse fails.
