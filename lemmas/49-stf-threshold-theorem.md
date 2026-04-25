# Lemma 49: STF Threshold Theorem

## Scope

- Status: Proven. This lemma is class-limited to genuinely new parity-even fully symmetric trace-free primitive families `Y_L` with rank `L \ge 3`.
- Status: Proven. It uses the fixed-order counting rule of [`../docs/power-counting.md`](../docs/power-counting.md).

## Formula Layer

- Status: Proven. By Lemma 48, no linear scalar witness exists for a single `Y_L` block.
- Status: Proven. By Lemma 48, no mixed quadratic scalar `E Y_L` exists for `L \ge 3`.
- Status: Proven. The quadratic norm `Y2` exists and has weight

```math
W_{\mathrm{self}}(L \ge 3; w_Y) = 2 w_Y.
```

- Status: Proven. Therefore the universal threshold layer for this class is

```math
W_{\min}(L \ge 3; w_Y) = 2 w_Y.
```

- Status: Proven. This theorem is self-only rather than mixed-aware, because no lower mixed quadratic witness exists in the class.

## `Delta <= 4` Consequence

- Status: Proven. The current power-counting rule assigns each primitive species a positive intrinsic weight `w_Y \ge 1`.
- Status: Proven. Minimal-sector uniqueness at `\Delta \le 4` can survive a genuine STF family only if its first unavoidable witness lies above the retained window:

```math
2 w_Y > 4.
```

- Status: Proven. Since the current counting rule is integer-valued on primitive species, this becomes

```math
w_Y \ge 3.
```

## Audited Checks

- Status: Proven. The audited ranks `L = 3, 4, 5, 6` all satisfy the same threshold consequence `w_Y \ge 3`.
- Status: Proven. Rank `4` remains a mixed-pattern exception because `EEQ` survives at the cubic layer.
- Status: Proven. That exception does not change the threshold formula, because it still lies above the quadratic self witness in the `w_Y = 1` ordering and above the protected window once `w_Y \ge 3`.

## Verdict

- Status: Proven. All genuinely new parity-even STF primitive families of rank `L \ge 3` satisfy the same universal self-witness threshold `w_Y \ge 3` at `\Delta \le 4`.
