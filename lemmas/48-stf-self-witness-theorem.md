# Lemma 48: STF Self-Witness Theorem

## Scope

- Status: Note. This lemma applies only to genuinely new parity-even fully symmetric trace-free primitive families `Y_L` with rank `L \ge 3`.
- Status: Proven. Derived tensors such as gradients of lower-rank families and trace descendants are excluded from the primitive-family definition and are not part of this lemma.

## Claim

- Status: Proven. No parity-even scalar linear in a single `Y_L` block exists.
- Status: Proven. No mixed quadratic scalar built from one electric block `E_{ij}` and one `Y_L` block exists when `L \ge 3`.
- Status: Proven. The quadratic norm `Y2` always exists.

## Proof

### Linear scalar

- Status: Proven. A scalar linear in `Y_L` would have to contract all `L` indices of a single STF block using only invariant tensors already allowed in the current primitive-family definition.
- Status: Proven. In the parity-even nonspinning free-fall sector, any parity-even scalar contraction with no additional primitive block reduces to contractions with Kronecker deltas.
- Status: Proven. Every delta contraction of two indices of the same STF block is a trace.
- Status: Proven. Because `Y_L` is trace-free, every such contraction vanishes.
- Status: Proven. Therefore no nonzero scalar linear in a single genuine STF block exists for any `L \ge 1`.

### Mixed quadratic `E Y_L`

- Status: Proven. Consider a candidate scalar built from one `E_{ij}` and one `Y_{k_1 \cdots k_L}`.
- Status: Proven. The electric block is itself STF, so contracting its two indices with each other gives `Tr(E) = 0`.
- Status: Proven. Therefore each electric index must contract into an index of `Y_L`.
- Status: Proven. After those two contractions, `L - 2` indices of `Y_L` remain unpaired.
- Status: Proven. To obtain a scalar with no other primitive block, those remaining indices must contract among themselves through deltas.
- Status: Proven. Every such self-contraction is a trace of `Y_L`, hence zero because `Y_L` is STF.
- Status: Proven. Therefore no nonzero mixed quadratic scalar `E Y_L` exists for `L \ge 3`.

### Quadratic norm

- Status: Proven. Two copies of the same STF block can always be fully cross-contracted:

```math
Y2 := Y_{i_1 \cdots i_L} Y^{i_1 \cdots i_L}.
```

- Status: Proven. This contraction does not use any internal trace of a single copy of `Y_L`.
- Status: Proven. Hence it is nonzero in the generic algebraic class and defines a parity-even scalar witness.

## Rank-2 Exception

- Status: Proven. The previous mixed-quadratic argument fails at `L = 2`, because after contracting the two indices of `E_{ij}` into `X_{kl}`, no leftover indices remain.
- Status: Proven. Therefore `E_{ij} X^{ij}` exists for the rank-2 STF class.
- Status: Note. This is exactly why rank `2` must remain outside the present `L \ge 3` theorem.

## Consequence

- Status: Proven. For every genuine parity-even STF primitive family with `L \ge 3`, the first unavoidable family witness is at least quadratic in the family block.
