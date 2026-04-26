# nonanalytic activation loophole

## Smallest Explicit Local Model

Status: Proven. The smallest explicit local nonanalytic counterexample to the analytic monopole-jet step uses a single corrected scalar normal-form coordinate `Y` and the smooth flat response

```math
m_A(Y)
=
m_0 + \alpha\,\phi_{\mathrm{flat}}(Y),
\qquad
\phi_{\mathrm{flat}}(Y)=
\begin{cases}
0, & Y \le 0,\\
e^{-1/Y^2}, & Y > 0.
\end{cases}
```

## Why It Escapes Theorem A

- Status: Proven. The response is local and keeps the finite admitted family/operator envelope intact.
- Status: Proven. Every Taylor coefficient about `Y=0` vanishes beyond the constant term, so the formal Taylor jet is just `m_0`.
- Status: Proven. Nevertheless `m_A(Y) > m_0` for every `Y>0`, so no finite analytic Taylor jet captures the response in a punctured neighborhood of the reference point.
- Status: Proven. The violated theorem assumption is A5.

## Rougher Comparison Model

Status: Proven. A stronger but less minimal threshold activation model is

```math
m_A(Y)
=
m_0 \left[
1 + a \, \Theta(Y - Y_c) \sqrt{Y - Y_c}
\right].
```

- Status: Proven. This also violates `A5` by creating a branch point at `Y=Y_c`.
- Status: Proven. It is not the smallest counterexample because it introduces explicit activation-point data `Y_c` and a rougher singularity.
