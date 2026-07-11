# Nonanalytic Model Table

| Model class | Locality kept? | Analytic near reference? | Finite jet valid? | Generalized replacement data | Theorem layer broken |
| --- | --- | --- | --- | --- | --- |
| Analytic control `m_A(Y)=m_0+\alpha Y+\beta Y^2` | Yes | Yes | Yes | Ordinary finite Taylor coefficients | None |
| Smooth flat one-coordinate activation `m_A(Y)=m_0+\alpha \phi_{\mathrm{flat}}(Y)` with `\phi_{\mathrm{flat}}(Y)=0` for `Y\le 0` and `e^{-1/Y^2}` for `Y>0` | Yes | No | No | Finite branch-profile data only if the model class is fixed explicitly; otherwise a non-Taylor function germ is needed | Analytic monopole jet collapse (`A5`, Lemma 55) |
| Threshold square-root activation `m_A(Y)=m_0+\alpha \Theta(Y-Y_c)\sqrt{Y-Y_c}` | Yes | No | No | Threshold location `Y_c`, branch exponent `1/2`, amplitude `\alpha`, and branch label | Analytic monopole jet collapse (`A5`, Lemma 55) |

## Readout

- Status: Proven. All three models keep locality and finite family/operator closure separate from analyticity.
- Status: Proven. Only the analytic control collapses to an ordinary finite sensitivity jet.
- Status: Proven. The smooth flat model is the smallest sharp counterexample because it breaks analyticity without introducing nonlocality or a discontinuous activation rule.
- Status: Proven. The threshold model is a stronger but less minimal nonanalytic escape route.
