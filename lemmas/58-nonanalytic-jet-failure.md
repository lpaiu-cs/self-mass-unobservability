# Lemma 58: Nonanalytic Jet Failure

## Statement

- Status: Proven. Let `Y` be any corrected scalar normal-form coordinate of positive operator weight in the finite reduced basis.
- Status: Proven. Consider the one-coordinate local monopole model

```math
m_A(Y) = m_0 + \alpha \,\phi_{\mathrm{flat}}(Y),
\qquad
\phi_{\mathrm{flat}}(Y) =
\begin{cases}
0, & Y \le 0,\\
e^{-1/Y^2}, & Y > 0.
\end{cases}
```

- Status: Proven. This model preserves locality, finite admitted family content, fixed-order operator finiteness, and finite scalar normal-form closure.
- Status: Proven. It violates only analyticity `A5`.
- Status: Proven. Therefore it is the smallest explicit local counterexample to the analytic Taylor-jet step of [`55-monopole-jet-collapse.md`](55-monopole-jet-collapse.md).

## Proof

1. Status: Proven. The model is local because `m_A` depends only on the instantaneous scalar coordinate `Y(\tau)` and not on any hereditary kernel or orbital-timescale internal state variable.
2. Status: Proven. The model introduces no new primitive family and no new local operator class. It changes only the functional dependence of the monopole coefficient on the already admitted coordinate `Y`.
3. Status: Proven. The function `\phi_{\mathrm{flat}}` is `C^\infty` at `Y = 0`, and every derivative at `Y = 0` vanishes.
4. Status: Proven. Hence every Taylor coefficient of `m_A(Y)` about `Y = 0` beyond the constant term is zero, so the formal Taylor jet equals the constant `m_0`.
5. Status: Proven. But for every `Y > 0`, one has `\phi_{\mathrm{flat}}(Y) > 0`, so `m_A(Y) > m_0`.
6. Status: Proven. Therefore no finite analytic Taylor jet about `Y = 0` captures the local monopole response in any punctured neighborhood of the reference point.
7. Status: Proven. The exact failing step in [`55-monopole-jet-collapse.md`](55-monopole-jet-collapse.md) is the use of analyticity `A5` to assert the existence of a Taylor expansion that represents the response.

## Comparison Models

- Status: Proven. The cusp model `m_A(Y)=m_0+\alpha |Y|` is an even simpler nonsmooth counterexample, but it breaks differentiability as well as analyticity.
- Status: Proven. The threshold model `m_A(Y)=m_0+\alpha \Theta(Y-Y_c)\sqrt{Y-Y_c}` is another explicit local counterexample, but it introduces extra activation-point data `Y_c`.
- Status: Proven. The smooth flat model above is the sharpest smallest counterexample because it preserves smoothness and locality while breaking only analyticity.
