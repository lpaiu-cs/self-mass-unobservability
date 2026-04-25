# Nonlinear Second-Order Internal Mode

Status: Counterexample candidate. This note starts M12: a minimal nonlinear
extension of the second-order internal state. The target is not sideband
existence. The target is whether generated lines and linear lines share the
same quadratic internal denominator after static nonlinear comparator budgets
are counted.

## Minimal Model

Status: Counterexample candidate. The nonlinear-drive branch is

```math
\mu_\chi \ddot\chi
+\gamma_\chi \dot\chi
+\omega_\chi^2(\chi-\alpha F)
=\beta_{F2}F^2 .
```

For derivative variable `z`,

```math
D_2(z)=\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2 .
```

Status: Counterexample candidate. With linear readout

```math
q=c_YF+c_\chi\chi ,
```

the linear transfer is

```math
G_L(z)
=
c_Y+
\frac{c_\chi\alpha\omega_\chi^2}{D_2(z)} .
```

For a generated sideband at `z=u+v`, the nonlinear-drive contribution is

```math
S_\beta(z)
=
\frac{c_\chi\beta_{F2}}{D_2(z)}\,\widehat{F^2}(z).
```

Status: Proven. A single generated line still does not count as a uniquely
dynamic observable: a static nonlinear coefficient multiplying `F^2` can match
one complex sideband sample.

## Nonlinear Readout Branch

Status: Counterexample candidate. The minimal nonlinear-readout branch is

```math
q
=c_YF+c_\chi\chi+\lambda_{F\chi}F\chi+\lambda_{\chi2}\chi^2 .
```

For two source frequencies `u` and `v`, ignoring the common product
`\hat F(u)\hat F(v)`, the generated line contains

```math
S_{\rm read}(u,v)
=
\lambda_{F\chi}\,[H_2(u)+H_2(v)]
+\lambda_{\chi2}H_2(u)H_2(v),
```

where

```math
H_2(z)=\frac{\alpha\omega_\chi^2}{D_2(z)}.
```

Status: Counterexample candidate. This branch can be stronger than pure
nonlinear drive because the generated line can inherit resonance from the
source frequencies, not only from the generated frequency. The new readout
coefficients `lambda_Fchi` and `lambda_chi2` are shared nuisance parameters
and must be counted.

## Shared-Denominator Diagnostic

Status: Proven. For the nonlinear-drive branch, after deprojection and after
removing the static direct piece `c_Y`, the dynamic model obeys

```math
[G_L(z_i)-c_Y]S_\beta(z_j)
-
[G_L(z_j)-c_Y]S_\beta(z_i)
=0 .
```

This is not a sideband-existence claim. It is a shared-denominator law: the
same `D_2(z)` organizes the linear response and the nonlinear generated
response.

Status: Proven. The diagnostic collapses if `c_\chi\alpha\omega_\chi^2=0`,
`c_\chi\beta_{F2}=0`, or if every frequency is given an independent static
nonlinear or projection nuisance.

Status: Counterexample candidate. If `c_Y` is not independently calibrated,
it is a shared nuisance parameter. It may be fitted jointly, but it cannot be
silently ignored in a sample-budget claim.

## Static Comparator Boundary

Status: Proven. A finite static nonlinear polynomial comparator

```math
P_M(z)=\sum_{m=0}^M b_m z^m
```

with `K` shared projection nuisances can absorb at most `M+1+K` complex
generated-sideband samples.

Status: Proven. A generated-sideband line-shape claim requires

```math
N_{\rm gen}\ge M+2+K
```

and, for a resonance claim, usable samples on both sides of the relevant
second-order denominator.

Status: Counterexample candidate. A joint nonlinear second-order claim must
either beat the finite generated-sideband budget and the finite linear
line-shape budget, or use a calibrated cross-ratio/diagnostic in which all
shared nuisance parameters are explicitly counted.

## Current Verdict

Status: Counterexample candidate. Minimal nonlinear second-order dynamics is a
stronger loophole candidate than one-state nonlinear sidebands because the
generated lines carry either the same quadratic denominator `D_2(z)` or
products/sums of the same `H_2` factors that also fit the linear response.

Status: Proven. This does not yet establish a unique observable. Underbudget
finite sideband data are still static nonlinear mimicry, not theorem-breaking
evidence.

Status: Counterexample candidate. The next useful test is a resonance-assisted
sideband sample design: count generated lines, linear lines, readout nuisance
coefficients, and projection nuisance, then ask whether the shared `D_2`
residual survives beyond that budget.
