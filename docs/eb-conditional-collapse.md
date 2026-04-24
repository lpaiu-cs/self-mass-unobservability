# `E/B` Conditional Collapse

## Statement

- Status: Conjectural. Fix `\Delta \le 4` and suppose the admissible local free-fall scalar operator space is reduced to the corrected finite `E/B` basis

```math
Y^{I}_{EB}
=
\{E2,\ B2,\ E3,\ EB2,\ E2^2,\ B2^2,\ dotE2,\ dotB2,\ EBDtB,\ E2B2,\ EB\_sq,\ TrE2B2,\ gradE2,\ divE2,\ mixedGradE2,\ gradB2,\ divB2,\ mixedGradB2\}.
```

- Status: Conjectural. Assume the body is quasi-static, nearly spherical, nonspinning, parity-even, self-bound in equilibrium, locally described by a worldline EFT, and carries no orbital-timescale internal state variable.
- Status: Conjectural. Assume the monopole coupling is analytic in the corrected `E/B` coordinates `Y^I_{EB}` near the reference background.

- Status: Conjectural. Then the free-fall action can be written as

```math
S_A
=
\int d\tau \left[
-m_A(Y_{EB})
+ \sum_\alpha C_{A,\alpha}\,{\cal O}_\alpha
\right],
```

with a finite Taylor jet

```math
m_A(Y_{EB})
=
m_A^{(0)}
+\sum_{n=1}^{N}
\frac{m_A^{(0)}}{n!}
s^{(EB)}_{A,I_1\cdots I_n}
Y^{I_1}_{EB}\cdots Y^{I_n}_{EB}
+ O(\Delta_{\max}+1),
```

so the leading body-dependent free-fall deviation still collapses to finitely many `E/B`-expanded sensitivity coordinates plus finitely many Wilson coefficients.

## What This Does And Does Not Claim

- Status: Proven. This note is conditional on the corrected `E/B` basis; it does not prove magnetic-family ordering.
- Status: Proven. This note does not rescue the electric-only minimal-sector theorem.
- Status: Conjectural. This note shows the broader finite-dimensional collapse program remains structurally plausible after explicit magnetic-family admission.

## Proof Skeleton

1. Status: Proven. [`../lemmas/11-eb-survivor-independence-delta4.md`](../lemmas/11-eb-survivor-independence-delta4.md) supplies a corrected finite `18`-element `E/B` basis.
2. Status: Conjectural. Use the no-internal-state assumption to exclude extra explicit orbital-timescale coordinates and hereditary functionals in the monopole sector.
3. Status: Conjectural. Expand the analytic monopole function `m_A(Y_{EB})` in the finitely many corrected `E/B` coordinates.
4. Status: Conjectural. Truncate at fixed order `\Delta \le 4`; only finitely many Taylor coefficients survive.
5. Status: Conjectural. Classify the remaining body dependence into higher-rank local operators with Wilson coefficients.

## Dependencies

- Status: Conjectural. [`conditional-collapse-lemma.md`](conditional-collapse-lemma.md)
- Status: Proven. [`../lemmas/11-eb-survivor-independence-delta4.md`](../lemmas/11-eb-survivor-independence-delta4.md)
- Status: Proven. [`magnetic-family-ordering.md`](magnetic-family-ordering.md)
- Status: Proven. [`reduction-rules.md`](reduction-rules.md)
