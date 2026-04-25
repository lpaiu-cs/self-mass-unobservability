# Dynamic Internal Visibility as a Loophole to Sensitivity Collapse

Status: Conjectural. This is a prose-first dynamic-loophole manuscript draft.
It is separate from `manuscript.md`, which records the older manuscript path.
Every scientific claim is labeled to preserve the repository contract.

## Abstract

Status: Conjectural. We study the theorem assumption that compact-body
free-fall deviations contain no orbital-timescale internal state variable.
A first-order relaxation state produces a phase-lagged shared pole, but under
the current orbital and two-tone forcing dictionaries it remains underbudget
against finite shared static comparators. We therefore treat the one-state
branch as a control/no-go result rather than as the final novelty claim.

Status: Counterexample candidate. We then analyze a second-order internal
mode. Its quadratic denominator survives finite shared free-fall-style
projection, and exact-in-e eccentric forcing supplies parameterized
resonance-aware detectability regions. With minimal nonlinear drive or
readout, generated harmonics can carry the same denominator, exceed finite
static nonlinear comparator budgets, and remain locally separable from the
linear response in a non-fine-tuned dimensionless parameter region.

Status: Proven. The result is not an instrument forecast and does not claim
identifiability under arbitrary per-frequency projection nuisance. It is a
theorem-level counterexample candidate to the no-internal-state assumption
A4.

## Front-Matter Claim Ledger

Status: Proven. This manuscript uses the following claim boundary.

| Claim | Status | Evidence route | Non-claim |
| --- | --- | --- | --- |
| One-state relaxation is a control/no-go branch for the current dictionaries. | Proven | Shared-pole response, nonlinear sideband mimicry, and sample-budget audits. | It is not the final positive dynamic signature. |
| Static nonlinear comparators can mimic sideband existence inside their interpolation budget. | Proven | Nonlinear comparator and sample-budget audits. | Sideband existence alone is not novelty. |
| Second-order internal visibility preserves a shared quadratic denominator under finite shared projection. | Counterexample candidate | Projection and resonance comparator audits. | It does not survive arbitrary per-frequency projection nuisance. |
| Exact-in-e eccentric forcing can supply budget-breaking samples for the second-order branch. | Counterexample candidate | Resonance-aware and amplitude-weighted detectability maps. | Formal harmonics do not count unless usable by amplitude or SNR. |
| Nonlinear second-order generated lines can share the same denominator as the linear lines. | Counterexample candidate | Nonlinear second-order mode and generated-line detectability audits. | A single generated line is still static-nonlinear-degenerate. |
| Generated components are locally separable in part of the parameter space. | Counterexample candidate | Joint harmonic Jacobian rank audit. | Local separability is not enough if generated samples are underbudget. |
| The positive region has analytic phase boundaries. | Proven | Count, bracket, amplitude-usability, and rank boundaries reproduce `162/162` grid verdicts. | The grid fraction is not a physical population measure. |

## Figure And Table Plan

Status: Conjectural. The manuscript should use four main figures and three
front-loaded tables.

| Item | Status | Purpose |
| --- | --- | --- |
| Figure 1 | Conjectural | Assumption map: static sensitivity collapse assumptions and the A4 dynamic violation. |
| Figure 2 | Proven | One-state branch flow: phase lag, sideband mimicry, and sample-budget no-go. |
| Figure 3 | Counterexample candidate | Second-order denominator geometry: resonance bracketing and finite shared projection. |
| Figure 4 | Counterexample candidate | Nonlinear second-order phase map: robust-positive, underdesigned, and degenerate regions. |
| Table 1 | Proven | Comparator budgets for linear, nonlinear, and projection nuisance models. |
| Table 2 | Counterexample candidate | Default robustness-map counts by projection channel. |
| Table 3 | Proven | Analytic phase-boundary conditions and collapse verdicts. |

## 1. Introduction

Status: Imported from prior work. The static sensitivity-collapse program
asks whether leading body-dependent free-fall deviations reduce to a finite
sensitivity-manifold EFT under explicit assumptions. In that program,
minimal-sector uniqueness is fragile once magnetic and scalar primitive
families are admitted, while broader fixed-order finite-family collapse
remains the positive theorem target.

Status: Counterexample candidate. The present manuscript attacks a different
assumption. It asks what happens when A4, the absence of orbital-timescale
internal state variables in the free-fall sector, is dropped.

Status: Proven. This change matters because more static primitive-family
audits do not test A4. A static family audit can enlarge or close a finite
basis at fixed order, but it does not decide whether a body carries memory or
internal dynamics on the orbital timescale.

Status: Counterexample candidate. The dynamic target is therefore not a new
static coefficient. The target is an observable relation that requires a
shared internal transfer structure across multiple frequency samples or
generated lines.

Status: Proven. A finite line or a finite ratio is not enough by itself. If a
finite shared static comparator has enough coefficients to interpolate the
observations, the purported dynamic observable remains a no-go for that
design.

Status: Counterexample candidate. The final positive branch in this
manuscript is nonlinear second-order internal visibility. It supplies a
shared quadratic denominator, a generated-line sector, component
separability, a finite robustness map, and analytic phase boundaries.

## 2. Comparator Discipline

Status: Proven. The comparator discipline is the central methodological
constraint. A dynamic signal counts only after the relevant finite shared
comparator and finite shared projection nuisance budgets have been exceeded.

Status: Proven. For a linear shared polynomial comparator of degree `N`, the
linear sample budget is exceeded only after more than `N+1` complex frequency
samples once projection nuisance has been counted.

Status: Proven. For a static nonlinear generated-line comparator of degree
`M`, the generated sample budget is

```math
H_{\rm count}=M+2+K_\Lambda,
```

where `K_Lambda` is the number of shared projection nuisance parameters.

Status: Proven. A projection function that is arbitrary at every frequency
collapses every finite-sample dynamic claim. The projection must be known,
calibrated, or finite-dimensional with shared nuisance parameters.

Status: Proven. This comparator rule is why one-frequency quadrature,
sideband existence, and a single generated line are not accepted as final
novelty claims.

## 3. One-State Relaxation As Control/No-Go

Status: Counterexample candidate. The smallest A4-violating model is a
one-state relaxation coordinate:

```math
\tau_\chi \dot\chi+\chi=\alpha F(Y),
\qquad
m_A=m_A^{(0)}[1+c_YF(Y)+c_\chi\chi].
```

Status: Proven. Under monochromatic forcing `F(t)=F_0 cos(Omega t)`, the
steady response has transfer

```math
H_1(i\Omega)=\frac{\alpha}{1+i\Omega\tau_\chi}.
```

Status: Proven. This creates an in-phase component, a quadrature component,
and a phase lag. It collapses to the static sensitivity description when
`tau_chi=0`, `Omega=0`, or `alpha c_chi=0`.

Status: Proven. In the adiabatic limit `Omega tau_chi << 1`, the response is
absorbed order by order into a local derivative EFT.

Status: Proven. One-frequency quadrature is not decisive if the comparator
allows a local derivative coefficient. The stronger one-state target is a
shared `tau_chi` relation across multiple frequencies.

Status: Proven. Minimal nonlinear one-state extensions can generate
sidebands, but sideband existence itself is static-nonlinear-degenerate.
Static terms such as `a_2 F^2` can generate `2Omega`,
`Omega_1+Omega_2`, and `|Omega_1-Omega_2|` inside their finite interpolation
budget.

Status: Proven. The current orbital dictionary `n,2n -> 3n` and the current
two-tone dictionary remain underbudget against realistic comparator classes.
Thus the one-state branch is a control/no-go branch for this manuscript.

## 4. Second-Order Internal Mode

Status: Counterexample candidate. The next minimal dynamic model is a
second-order internal mode:

```math
\mu_\chi\ddot\chi+\gamma_\chi\dot\chi
+\omega_\chi^2(\chi-\alpha F)=0.
```

Status: Counterexample candidate. Its transfer denominator is

```math
D_2(z)=\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2,
\qquad
H_2(z)=\frac{\alpha\omega_\chi^2}{D_2(z)}.
```

Status: Proven. The first-order relaxation model is recovered as a collapse
limit when `mu_chi=0` after parameter matching. In an adiabatic sampled band,
the second-order response also collapses into local derivative coefficients.

Status: Counterexample candidate. Away from those limits, the second-order
mode supplies resonance line shape, phase wrapping, and a shared quadratic
denominator rather than a single relaxation pole.

Status: Proven. A resonance peak or phase sign flip alone is still not a
novelty claim. The samples must exceed the finite shared comparator budget
and must bracket the denominator feature.

## 5. Free-Fall-Style Projection

Status: Counterexample candidate. The branch uses two concrete
free-fall-style readouts:

```math
\delta\hat a=\Gamma\hat q,
\qquad
\ddot R+\kappa^2R=\Gamma q_A,
\qquad
\Lambda_R(z)=\frac{\Gamma}{\kappa^2+z^2}.
```

Status: Counterexample candidate. If `Gamma` and `kappa` are known,
calibrated, or finite shared nuisance parameters, the projection does not
become an arbitrary complex number at each sampled frequency.

Status: Counterexample candidate. Under those finite-nuisance conditions,
the internal denominator can survive in the observable channel.

Status: Proven. The projection claim collapses if the projection is allowed
to be an arbitrary per-frequency complex nuisance, if `Gamma=0`, if the
sampled point hits a projection singularity, or if the projection cancels the
internal denominator.

## 6. Exact-In-E Eccentric Forcing

Status: Counterexample candidate. The concrete forcing dictionary is

```math
F(Y(t))=(a/r(t))^p.
```

Status: Counterexample candidate. Its exact-in-e harmonic coefficients are

```math
A_k(p,e)=
\frac{1}{\pi}\int_0^{2\pi}
(1-e\cos E)^{1-p}
\cos[k(E-e\sin E)]\,dE.
```

Status: Proven. Single-system positive harmonics can bracket the internal
phase-wrap location only when

```math
\rho=\Omega_{\rm wrap}/n>1.
```

Status: Proven. The bracketing cutoff is

```math
H_{\rm bracket}(\rho)=
\begin{cases}
\lfloor\rho\rfloor+1,& \rho>1,\ \rho\notin\mathbb Z,\\
\rho+1,& \rho\in\mathbb Z,\ \rho>1,\\
\hbox{none},& \rho\le1.
\end{cases}
```

Status: Proven. The formal harmonic minimum is

```math
H_{\min}=\max\{N+2+K_\Lambda,H_{\rm bracket}(\rho)\}
```

for the linear second-order branch, with the analogous generated-line budget
for the nonlinear branch.

Status: Counterexample candidate. Formal harmonics do not count unless they
are usable. The relative cutoff and SNR versions replace mere nonzero
coefficients with `W_k/W_max >= eta` or `SNR_k >= SNR_min`.

Status: Counterexample candidate. In the parameterized physical detectability
map, moderate eccentricity designs such as `p=2`, `e=0.3`, `rho=3/2`,
`delta=0.2`, and `H=6` can cross the sample budget under finite shared
projection nuisance, while `e=0.1` remains underbudget in the default audits.

## 7. Nonlinear Second-Order Generated Lines

Status: Counterexample candidate. The minimal nonlinear second-order drive is

```math
\mu_\chi\ddot\chi+\gamma_\chi\dot\chi
+\omega_\chi^2(\chi-\alpha F)=\beta_{F2}F^2.
```

Status: Counterexample candidate. A nonlinear readout alternative is

```math
m/m_0=1+c_YF+c_\chi\chi+\lambda_{F\chi}F\chi+\lambda_{\chi2}\chi^2.
```

Status: Counterexample candidate. The important point is not merely that
generated harmonics exist. The important point is that generated components
can carry the same denominator `D_2(z)` that organizes the linear response.

Status: Counterexample candidate. For `F=(a/r)^p`, the nonlinear generated
coefficient is represented by

```math
B_k(p,e)=A_k(2p,e).
```

Status: Proven. A generated-line claim must exceed

```math
|U|\ge M+2+K_\Lambda
```

for a degree-`M` static nonlinear comparator and finite shared projection
nuisance count `K_Lambda`.

Status: Counterexample candidate. In the default generated-line audit with
`p=2`, `rho=3/2`, `delta=0.2`, `H_gen=6`, and `M=1`, `e=0.3` and `e=0.6`
are generated-budget-breaking under the relative cutoff, while `e=0.1`
remains generated-underbudget.

Status: Counterexample candidate. In the SNR example with `SNR_min=5`, the
default moderate-eccentricity minimum generated-amplitude scales are
`Q_beta=0.429654` for the acceleration channel and `Q_beta=0.118468` for the
range channel.

## 8. Component Separability

Status: Counterexample candidate. Exact-in-e linear and nonlinear generated
contributions can occupy the same harmonic. The observable model must
therefore separate the generated component from the linear component.

Status: Counterexample candidate. The joint harmonic model is

```math
\hat O_k=
\Lambda_k\left[
Q_YA_k(p,e)
+Q_LA_k(p,e)H_2(k)
+Q_\beta B_k(p,e)D_2(k)^{-1}
\right].
```

Status: Proven. A generated component is locally separable only if the
Jacobian column for `Q_beta` adds rank relative to the model without
`Q_beta`.

Status: Proven. A full local separability design also requires the real
Jacobian rank to equal the number of shared parameters in the joint model.

Status: Counterexample candidate. At `e=0.3`, the acceleration channel is
component-separable-and-budget-ready with minimum cutoff `H=4`, while fixed
and free `kappa` range channels are component-separable-and-budget-ready with
minimum cutoff `H=5`.

Status: Counterexample candidate. At `e=0.1`, acceleration and fixed-`kappa`
range can be locally separable but generated-budget-underdesigned, while
free-`kappa` range is generated-component-degenerate in the default audit.

## 9. Robustness And Analytic Phase Boundaries

Status: Counterexample candidate. The nonlinear robustness map tests whether
the positive branch is isolated. The default grid is

```text
p in {1.5, 2, 3}
e in {0.1, 0.3, 0.6}
rho in {4/3, 3/2, 2}
delta in {0.1, 0.2}
projection in {acceleration, range-fixed-kappa, range-free-kappa}
H = 6
eta = 10^-3
M = 1
```

Status: Counterexample candidate. The grid has `162` points. The audit finds
`107` robust-positive points, `40` component-separable-but-budget-underdesigned
points, and `15` generated-component-degenerate points.

Status: Counterexample candidate. By projection channel, the robust-positive
counts are `37/54` for acceleration, `35/54` for range with fixed `kappa`,
and `35/54` for range with free `kappa`.

Status: Proven. The positive fraction is not a physical prior or population
measure. It is a stress test showing that the positive branch is not confined
to a single tuned point in the chosen dimensionless grid.

Status: Proven. The analytic phase-boundary theorem explains the finite grid
using four boundaries:

1. Status: Proven. Count: `H_count=M+2+K_Lambda`.
2. Status: Proven. Bracket: `rho>1` and usable harmonics on both sides of
   `rho`.
3. Status: Counterexample candidate. Amplitude usability: `W_k/W_max >= eta`
   or an SNR analogue.
4. Status: Counterexample candidate. Rank: `Q_beta` adds full shared-parameter
   Jacobian rank.

Status: Proven. The boundary audit reproduces all `162/162` default
robustness-map verdicts.

## 10. Discussion

Status: Counterexample candidate. The novelty of this branch is not another
static primitive family. It is an explicit A4-violating internal-state model
whose observable content can survive finite shared projection and finite
static comparator budgets in part of parameter space.

Status: Proven. The one-state branch is valuable because it prevents
overclaiming. It shows that phase lag, sidebands, and even shared-pole ratios
are insufficient under finite sample budgets unless the design beats the
comparator.

Status: Counterexample candidate. The nonlinear second-order branch is the
current strongest dynamic counterexample candidate because it combines a
shared quadratic denominator, generated lines, budget breaking, local
component separability, finite-grid robustness, and analytic phase boundaries.

Status: Proven. This manuscript does not claim real-system detectability.
The design variables `Q_0`, `Q_beta`, `sigma_0`, `Gamma`, `kappa`, cadence,
and harmonic extraction error remain unfixed until a system-anchored support
note is written.

Status: Proven. This manuscript does not claim identifiability under
arbitrary per-frequency projection nuisance. Finite shared nuisance is an
essential assumption.

Status: Conjectural. The natural next support artifact is a system-anchored
detectability note for one concrete free-fall-style channel. That note would
support the manuscript but should not replace the theorem-level claim.

## Appendix Plan

Status: Conjectural. Appendix A should collect the one-state response,
adiabatic collapse, and sample-budget no-go.

Status: Conjectural. Appendix B should collect the exact-in-e coefficient
definition and the harmonic usability criteria.

Status: Conjectural. Appendix C should collect the projection formulas for
acceleration-like and range-like channels.

Status: Conjectural. Appendix D should collect the Jacobian rank and
rank-minor calculations used in the component-separability and analytic
phase-boundary audits.
