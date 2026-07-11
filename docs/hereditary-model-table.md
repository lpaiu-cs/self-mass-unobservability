# Hereditary Model Table

| Model class | Locality kept? | Instantaneous analyticity kept? | Finite primitive-family envelope kept? | Finite local jet in `Y^I(\tau)`? | Finite-state Markovianizable? | Generalized replacement data | Theorem layer broken |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Local analytic control `m_A(Y)=m_0+\alpha Y+\beta Y^2` | Yes | Yes | Yes | Yes | Not needed | Ordinary finite Taylor coefficients | None |
| Exponential memory control `m_A[Y](\tau)=m_0+\alpha Y(\tau)+\lambda_{\exp}\int_{-\infty}^{\tau} d\tau' e^{-(\tau-\tau')/T_h}Y(\tau')/T_h` | No | Yes | Yes | No | Yes | Finite auxiliary-state data `(m_0,\alpha,\lambda_{\exp},T_h,\chi_A)` | Local monopole reduction in `Y^I(\tau)` alone; collapses back into an `A4`-type finite-state extension |
| Power-law hereditary kernel `m_A[Y](\tau)=m_0+\alpha Y(\tau)+\lambda_{\gamma}\int_{-\infty}^{\tau} d\tau' K_{\gamma}(\tau-\tau')Y(\tau')` with `0<\gamma<1` | No | Yes | Yes | No | No | Kernel branch data `(m_0,\alpha,\lambda_{\gamma},\gamma,T_h)` or a continuous spectral density | Local monopole reduction to `m_A(Y(\tau))` fails, so Lemmas 55 and 56 no longer apply |

## Readout

- Status: Proven. The local analytic control remains inside the proved theorem domain.
- Status: Proven. The exponential kernel is not a sharp `A3` escape because its memory can be encoded by finitely many local states.
- Status: Proven. The power-law kernel is the smallest genuinely hereditary counterexample because it keeps one-coordinate analytic instantaneous dependence while requiring nonlocal kernel data that are not finitely state-extendable.
