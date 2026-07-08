# Lemma 27: Rank-1 Threshold Formula

## Formula Layer

- Status: Proven. Let `w_V` be the primitive weight assigned to a genuine parity-even vector family `V_i`.
- Status: Proven. The first self witness is `V2`, so

```math
W_{\mathrm{self}}(R1; w_V) = 2 w_V.
```

- Status: Proven. The first mixed witness is `EVV`, so

```math
W_{\mathrm{mix}}(R1; w_V) = 1 + 2 w_V.
```

- Status: Proven. Therefore

```math
W_{\min}(R1; w_V) = 2 w_V.
```

## `\Delta_{\max} = 4` Threshold

- Status: Proven. Minimal-sector uniqueness up to `\Delta \le 4` requires

```math
2 w_V > 4,
```

so the current class-limited necessary threshold is

```math
w_V \ge 3.
```

- Status: Proven. The mixed branch gives `1 + 2 w_V > 4`, which also first clears the MVP window at `w_V \ge 3`, but the unsuppressed witness ordering remains self-first because `2 < 3`.

## Classification

- Status: Proven. The current unsuppressed vector-family obstruction is self-dominated because `W_{\mathrm{self}} = 2 < 3 = W_{\mathrm{mix}}`.
- Status: Proven. The current `R1` threshold is therefore self-only and sharp inside the audited vector-family class.
- Status: Proven. No mixed-aware rescue statement is needed for the present `R1` audit because no mixed witness appears below `V2`.

## Boundary

- Status: Note. This formula applies only to a genuinely new primitive vector family.
- Status: Proven. It does not apply to derivative-generated vectors already attached to `R0b` or `R2`.
