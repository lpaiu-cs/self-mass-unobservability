> SUPPORTING THEOREM NOTE
>
> Status: Proven. The authoritative front-door theorem statement is [`theorem-package.md`](theorem-package.md).

# STF Self-Witness Threshold Theorem

## Scope

- Status: Proven. This note is class-limited to genuinely new parity-even fully symmetric trace-free primitive families `Y_L` of rank `L \ge 3` in the MVP free-fall sector at fixed order `\Delta \le 4`.
- Status: Proven. It does not claim anything about mixed-symmetry, antisymmetric, trace-descended, or derivative-generated higher-rank structures.
- Status: Proven. It is a threshold theorem only.
- Status: Proven. It is not a theorem about the full first mixed-witness pattern.

## Theorem Statement

- Status: Proven. For every genuinely new parity-even STF primitive family `Y_L` with `L \ge 3`, no parity-even scalar linear in a single `Y_L` block exists.
- Status: Proven. For every genuinely new parity-even STF primitive family `Y_L` with `L \ge 3`, no mixed quadratic scalar built from one electric block `E_{ij}` and one `Y_L` block exists.
- Status: Proven. For every genuinely new parity-even STF primitive family `Y_L` with `L \ge 3`, the quadratic norm

```math
Y2 := Y_{i_1 \cdots i_L} Y^{i_1 \cdots i_L}
```

exists and is parity-even.
- Status: Proven. Therefore the first unavoidable self witness in this class occurs at weight

```math
W_{\mathrm{self}}(L \ge 3; w_Y) = 2 w_Y.
```

- Status: Proven. Because the fixed-order counting rule assigns positive intrinsic weight `w_Y \ge 1` and retains only operators with `\Delta \le 4`, minimal-sector uniqueness can survive this class only if

```math
2 w_Y > 4,
```

hence

```math
w_Y \ge 3.
```

## What This Theorem Does Not Say

- Status: Proven. This theorem does not say that the first mixed witness is always `EYY`.
- Status: Proven. This theorem does not say that all higher-rank tensor families fall into one universal mixed-pattern class.
- Status: Proven. This theorem does not rescue minimal-sector uniqueness.
- Status: Proven. It only identifies a universal necessary threshold layer for the genuine parity-even STF branch with `L \ge 3`.

## Audited Instances

- Status: Proven. The audited ranks `L = 3, 4, 5, 6` all satisfy the same self-witness threshold structure: first self witness `Y2`, no mixed quadratic `EY_L`, and threshold `w_Y \ge 3` at `\Delta \le 4`.
- Status: Proven. Rank `4` remains a mixed-pattern exception because `EEQ` survives at the first mixed cubic layer.
- Status: Proven. That exception does not affect the threshold theorem, because it does not create a witness below the quadratic norm layer.

## Supporting Check

- Status: Proven. [`../symbolic/stf_self_witness_check.py`](../symbolic/stf_self_witness_check.py) records the audited ranks `L = 3, 4, 5, 6` as explicit checks of the class theorem.
