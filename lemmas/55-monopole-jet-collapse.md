# Lemma 55: Monopole Jet Collapse

## Statement

- Status: Proven. Let `Y^I` be any finite basis of the reduced scalar normal-form space `{\cal N}({\cal F}, \Delta \le 4)`.
- Status: Proven. Let each coordinate `Y^I` carry a positive operator weight `\Delta_I \ge 1`.
- Status: Proven. Assume locality `A3` and analyticity `A5`, so that the monopole response is a local analytic function `m_A(Y)` near the reference background.
- Status: Proven. Then the fixed-order monopole sector collapses to finitely many Taylor coefficients at `\Delta \le 4`.

## Proof

1. Status: Proven. Because the basis `Y^I` is finite, there are only finitely many coordinate directions in the monopole argument list.
2. Status: Proven. Because `m_A(Y)` is analytic by `A5`, it admits a Taylor expansion in the coordinates `Y^I`.
3. Status: Proven. The retained multi-indices are exactly those satisfying

```math
\sum_I n_I \Delta_I \le 4.
```

4. Status: Proven. Since the set of coordinates is finite and each `\Delta_I \ge 1`, only finitely many multi-indices satisfy that bound.
5. Status: Proven. Therefore the truncated jet

```math
m_A(Y)
=
m_A^{(0)}
+
\sum_{\sum_I n_I \Delta_I \le 4}
\frac{m_{A,\mathbf{n}}}{\mathbf{n}!}
\prod_I (Y^I)^{n_I}
+
O(\Delta > 4)
```

contains only finitely many coefficients `m_{A,\mathbf{n}}`.

## Sensitivity Interpretation

- Status: Proven. The monopole Taylor coefficients `m_{A,\mathbf{n}}` are the sensitivity coordinates of the finite response manifold at fixed order.
- Status: Proven. They are not Wilson coefficients.
- Status: Proven. They depend on analyticity `A5`, not on minimal-sector uniqueness.

## Exact Failure Point If The Lemma Is Dropped

- Status: Counterexample candidate. If analyticity `A5` fails, the exact failing step is the existence of the finite Taylor jet, not the finiteness of the operator space itself.
- Status: Counterexample candidate. The current explicit loophole model for this step remains the nonanalytic activation branch in [`../counterexamples/nonanalytic-activation/README.md`](../counterexamples/nonanalytic-activation/README.md).
