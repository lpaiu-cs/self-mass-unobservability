# Lemma 15: Rank-2 STF Family Admission

## Family Class

- Status: Proven. Let `X_ij` be an unsuppressed symmetric trace-free rank-2 primitive family admitted at weight `1` into the parity-even free-fall scalar sector.
- Status: Proven. No extra `X`-specific identity, ordering rule, or background restriction is assumed beyond the currently stated reduction rules in [`../docs/reduction-rules.md`](../docs/reduction-rules.md).
- Status: Proven. This class includes parity-even STF families directly and covers the already audited magnetic-family pattern at the level of the parity-even scalar sector through the explicit witness `B2`.

## Generic Witness

- Status: Proven. The quadratic invariant

```math
X2 := Tr(X^2)
```

is a parity-even scalar of weight `2`.
- Status: Proven. `X2` is not a worldline total derivative because it contains no `D_\tau`.
- Status: Proven. `X2` is not removed by lower-order free-fall EOM because it carries no explicit acceleration factor.
- Status: Proven. `X2` is not removed by the stated Cayley-Hamilton-type identities, which first reduce quartic trace structures rather than the quadratic norm.
- Status: Proven. Therefore `X2` is a surviving low-order operator in the admitted family-conditioned normal form unless an extra explicit assumption kills it.

## Consequence

- Status: Proven. Any unsuppressed STF rank-2 family of this class obstructs minimal-sector uniqueness.
- Status: Proven. The magnetic-family audit provides the explicit audited instance `X = B`, for which the witness is `B2`.
- Status: Conjectural. Additional mixed operators such as `EX` or `Tr(EX^2)` may appear for some subclasses, but the quadratic witness `X2` already suffices for the no-go.

## Boundary

- Status: Proven. This lemma does not say that fixed-order finiteness fails.
- Status: Proven. It only says that the stronger minimal-sector theorem is not stable under unsuppressed admission of the stated rank-2 STF family class.
