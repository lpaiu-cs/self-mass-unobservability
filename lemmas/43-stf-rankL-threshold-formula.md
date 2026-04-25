# Lemma 43: STF Rank-`L` Threshold-Formula Attempt

## Target Formula

- Status: Conjectural. The attempted STF-tower threshold theorem would state that for every genuine parity-even STF primitive family `Y_L` with `L \ge 3`,

```math
W_{\mathrm{self}}(L \ge 3; w_Y) = 2 w_Y,\qquad
W_{\mathrm{mix}}(L \ge 3; w_Y) = 2 w_Y + 1,\qquad
W_{\min}(L \ge 3; w_Y) = 2 w_Y.
```

- Status: Conjectural. That would imply a single self-only sharp threshold
  `w_Y \ge 3`
  at `\Delta \le 4`.

## What Survives

- Status: Proven. The self formula

```math
W_{\mathrm{self}}(L \ge 3; w_Y) = 2 w_Y
```

  matches all audited STF ranks `L = 3, 4, 5`.
- Status: Proven. The odd-rank audited cases `L = 3` and `L = 5` also match

```math
W_{\mathrm{mix}}(w_Y) = 2 w_Y + 1,
```

  with canonical mixed witness `EYY`.

## Exact Failure

- Status: Proven. The uniform mixed formula fails at the first audited even-rank member `L = 4`.
- Status: Proven. The reason is the additional mixed cubic contraction `EEQ`, whose effective weight is

```math
W_{EEQ}(w_Q) = w_Q + 2.
```

- Status: Proven. The audited rank-4 mixed floor is therefore

```math
W_{\mathrm{mix}}(L = 4; w_Q) = \min(w_Q + 2,\ 2 w_Q + 1),
```

  not the claimed universal value `2 w_Q + 1`.
- Status: Proven. Hence the single STF-tower mixed formula for all `L \ge 3` is false.

## `\Delta \le 4` Consequence

- Status: Proven. Despite that failure, the audited `\Delta \le 4` threshold still comes out as
  `w_Y \ge 3`
  for the audited ranks `L = 3, 4, 5`.
- Status: Proven. For `L = 3` and `L = 5`, this follows from the odd-rank formula
  `W_{\min} = 2 w_Y`.
- Status: Proven. For `L = 4`, it follows from the corrected floor
  `W_{\min}(4; w_Q) = \min(2 w_Q,\ w_Q + 2)`,
  whose `\Delta \le 4` exclusion condition is still `w_Q \ge 3`.

## Verdict

- Status: Proven. The attempted theorem "all genuine parity-even STF rank-`L` primitive families with `L \ge 3` fall into one self-only sharp theorem class with mixed formula `2 w_Y + 1`" fails at `L = 4`.
- Status: Proven. The precise failing reason is the extra mixed cubic `EEQ`.
- Status: Proven. The later audited rank-6 result does not remove that failure; it only shows that the audited even-rank branch is split into a rank-4 exception plus a rank-6 return to the `EYY`-type mixed layer.
