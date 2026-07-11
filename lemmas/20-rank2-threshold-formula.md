# Lemma 20: Rank-2 Threshold Formula

## Setup

- Status: Proven. Let `X_ij` be a parity-even STF rank-2 family admitted on top of the electric baseline catalog.
- Status: Proven. Let `w_X` denote the counting weight assigned to the primitive `X_ij`.
- Status: Proven. The electric baseline block `E_ij` keeps weight `1`.
- Status: Proven. No explicit cross-family orthogonality rule such as `Tr(EX) = 0` is active in [`../docs/assumptions-ledger.md`](../docs/assumptions-ledger.md).

## Formula Derivation

- Status: Proven. The first self witness is

```math
X2 = Tr(X^2),
\qquad
W_{\mathrm{self}}(\mathrm{R2}; w_X) = 2 w_X.
```

- Status: Proven. The first mixed witness is

```math
EX = Tr(EX),
\qquad
W_{\mathrm{mix}}(\mathrm{R2}; w_X) = w_X + 1.
```

- Status: Proven. Therefore

```math
W_{\min}(\mathrm{R2}; w_X)
=
\min(2 w_X,\ w_X + 1).
```

- Status: Proven. The current unsuppressed audited case is `w_X = 1`, for which

```math
W_{\mathrm{self}} = W_{\mathrm{mix}} = W_{\min} = 2.
```

## Threshold Consequence At `\Delta_{\max}=4`

- Status: Proven. The self-only lower bound is obtained from `2 w_X > 4`, namely `w_X \ge 3`.
- Status: Note. This self-only bound is not the true current necessary threshold, because it does not remove the mixed witness `EX`.
- Status: Proven. The mixed-aware necessary threshold comes from

```math
\min(2 w_X,\ w_X + 1) > 4,
```

which for integer `w_X \ge 1` reduces to

```math
w_X + 1 > 4
\quad\Longrightarrow\quad
w_X \ge 4.
```

- Status: Proven. Therefore the current advertised rank-2 threshold is mixed-aware and sharp:
  `w_X \ge 4` unless an explicit `EX = 0`-type rule removes the mixed quadratic witness.

## Boundary

- Status: Proven. If one explicitly adds a cross-family rule killing `EX`, then the old self-only lower bound `w_X \ge 3` can reappear, but that rule is not active in the current repo.
- Status: Note. This lemma is class-limited to the audited parity-even STF rank-2 family.
- Status: Proven. The lemma obstructs the stronger minimal-sector uniqueness claim, not the positive finite-family collapse branch.
