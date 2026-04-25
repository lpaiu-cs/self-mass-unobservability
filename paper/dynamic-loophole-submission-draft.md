# Dynamic Internal Visibility as a Loophole to Sensitivity Collapse

## Abstract

Sensitivity-collapse arguments usually assume that the free-fall response of
a compact body can be represented by local, instantaneous sensitivity
coordinates and higher-multipole Wilson coefficients. This paper isolates one
assumption behind that conclusion: the absence of an orbital-timescale
internal state variable. We first analyze the smallest relaxation model. It
produces phase lag and a shared pole, but under the current orbital and
two-tone forcing dictionaries its finite line data remain underbudget against
finite shared static comparators. The one-state branch is therefore a control
and no-go branch, not the main positive claim.

We then study a second-order internal mode. Its quadratic denominator can
survive finite shared free-fall-style projection, and exact-in-e eccentric
forcing can supply resonance-aware harmonic samples. With minimal nonlinear
drive or readout, generated harmonics can carry the same denominator as the
linear response. In moderate-eccentricity, finite-nuisance designs these
generated lines can exceed static nonlinear comparator budgets, remain
locally separable from the linear response, and fall inside robustness regions
explained by analytic count, bracket, amplitude, and rank boundaries.

The result is a parameterized design theorem and a counterexample candidate
to the no-internal-state assumption. It is not an instrument forecast, not a
real-data detection claim, and not an identifiability claim under arbitrary
per-frequency projection nuisance.

## 1. Introduction

The static theorem program asks whether body-dependent deviations in the
free-fall sector collapse to a finite-dimensional sensitivity-manifold EFT.
That program has a natural static direction: audit the admissible primitive
families and determine when fixed-order closure survives. It also has a
different and sharper failure mode. If a compact body carries an
orbital-timescale internal state in the free-fall sector, the response need
not be an instantaneous function of local external invariants.

This paper studies that dynamic failure mode. The target is not another
static Wilson coefficient and not a larger static primitive basis. The target
is a shared internal transfer law that organizes multiple observable
frequency samples or generated lines in a way that a finite shared static
comparator cannot absorb.

The main lesson is deliberately conservative. A dynamic model is not novel
just because it produces a phase lag or a sideband. A finite static
comparator can interpolate finite data until its coefficient budget is
exceeded, and a projection model with arbitrary complex freedom at every
frequency can absorb any finite transfer anomaly. The paper therefore tracks
sample budgets, projection nuisance dimensions, generated-line amplitudes,
and local identifiability as part of the theorem statement.

Within that discipline, the one-state relaxation model closes as a
control/no-go branch under the current forcing dictionaries. The nonlinear
second-order internal mode remains the live counterexample candidate. It
supplies a shared quadratic denominator, generated-line observables,
component separability, robustness regions, and analytic phase boundaries.

## 2. Why The Static Line Motivates The A4 Attack

The static sensitivity-collapse theorem assumes, among other conditions, a
local worldline EFT and no orbital-timescale internal state variable in the
free-fall sector. Static primitive-family audits test the finite-basis side
of the theorem. They do not test whether a body can carry dynamical memory
that appears directly in the observable response.

The dynamic branch attacks this no-internal-state assumption directly. The
claim is not that all static collapse statements fail. The claim is that an
explicit internal state can violate the assumption set in a way that produces
observable transfer relations outside a purely instantaneous sensitivity
description, provided the resulting samples exceed the appropriate finite
comparator budgets.

This separation also explains the paper split. A static theorem or
classification paper should explain which finite-family collapse statements
survive and which minimal-sector uniqueness statements fail. The present
dynamic paper asks how the collapse can fail when an internal state is
allowed.

## 3. One-State Relaxation As The Control/No-Go Branch

The smallest dynamic model is a first-order relaxation coordinate,

```math
\tau_\chi \dot\chi+\chi=\alpha F(Y),
\qquad
m_A=m_A^{(0)}[1+c_YF(Y)+c_\chi\chi].
```

For a monochromatic drive, the response contains the transfer

```math
H_1(i\Omega)=\frac{\alpha}{1+i\Omega\tau_\chi}.
```

This produces in-phase and quadrature components and a phase lag. It also
has sharp collapse limits: the response becomes static when `tau_chi=0`,
`Omega=0`, or `alpha c_chi=0`, and it is absorbed order by order into a
local derivative EFT in the adiabatic regime `Omega tau_chi << 1`.

The one-state branch becomes scientifically useful because it prevents
overclaiming. A single quadrature component can be mimicked by a local
derivative coefficient. A single generated sideband can be mimicked by a
static nonlinear term. Even a shared one-pole relation has to be tested
against the finite interpolation budget of the chosen comparator.

Under the current orbital dictionary, which supplies `n` and `2n` and then
generates `3n`, the line and generated-sideband counts are underbudget. The
same conclusion holds for the current two-tone dictionary once realistic
first-derivative or first-degree static nonlinear comparators are admitted.
Thus the first-order model is not the paper's positive branch. It is the
control case that fixes the budget discipline.

## 4. Second-Order Internal Mode

The next dynamic model is a second-order internal mode,

```math
\mu_\chi\ddot\chi+\gamma_\chi\dot\chi
+\omega_\chi^2(\chi-\alpha F)=0.
```

Its transfer function contains the quadratic denominator

```math
D_2(z)=\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2,
\qquad
H_2(z)=\frac{\alpha\omega_\chi^2}{D_2(z)}.
```

This model has the expected collapse limits. When `mu_chi=0`, it reduces to a
first-order relaxation form after parameter matching. In an adiabatic sampled
band, the response again collapses into local derivative coefficients. Away
from those limits, however, the quadratic denominator gives a resonance-aware
line shape and phase wrapping that are richer than the one-state pole.

A resonance peak or a phase sign change is still not sufficient by itself.
The sampled frequencies must both exceed the comparator budget and bracket
the denominator feature. This turns the problem into a resonance-aware sample
design problem rather than a qualitative "peak exists" argument.

## 5. Projection Survival

The paper uses two free-fall-style projection channels:

```math
\delta\hat a=\Gamma\hat q,
\qquad
\ddot R+\kappa^2R=\Gamma q_A,
\qquad
\Lambda_R(z)=\frac{\Gamma}{\kappa^2+z^2}.
```

These channels are useful because their nuisance parameters can be finite and
shared. If `Gamma` and `kappa` are known, calibrated, or fitted as shared
parameters, the projection is not an arbitrary complex nuisance at each
frequency. In that setting the internal quadratic denominator can remain
visible in the projected observable.

The caveat is essential. If projection is allowed to vary freely at each
frequency, finite-sample identifiability collapses. The paper therefore makes
no claim under arbitrary per-frequency projection nuisance.

## 6. Exact-In-E Forcing And Detectability Region

The concrete orbital forcing template is

```math
F(Y(t))=(a/r(t))^p.
```

The exact-in-e harmonic coefficients are

```math
A_k(p,e)=
\frac{1}{\pi}\int_0^{2\pi}
(1-e\cos E)^{1-p}\cos[k(E-e\sin E)]\,dE.
```

For a single eccentric system with positive harmonics `Omega_k=kn`, the
phase-wrap ratio is `rho=Omega_wrap/n`. Bracketing is possible only for
`rho>1`. The formal bracketing cutoff is

```math
H_{\rm bracket}(\rho)=
\begin{cases}
\lfloor\rho\rfloor+1,& \rho>1,\ \rho\notin\mathbb Z,\\
\rho+1,& \rho\in\mathbb Z,\ \rho>1,\\
\hbox{none},& \rho\le1.
\end{cases}
```

For the linear second-order branch, the formal harmonic minimum is

```math
H_{\min}=\max\{N+2+K_\Lambda,H_{\rm bracket}(\rho)\}.
```

Formal harmonic count is not enough. Harmonics must clear a relative
amplitude cutoff or a parameterized SNR threshold. With `p=2`, `rho=3/2`,
`delta=0.2`, and `H=6`, the default audits find that low eccentricity
(`e=0.1`) remains underbudget while moderate eccentricity (`e=0.3`) can
clear the acceleration and range sample counts under the stated
dimensionless or parameterized SNR criteria. These are design examples, not
instrument forecasts.

## 7. Nonlinear Generated-Line Theorem

The minimal nonlinear second-order extension can be placed either in the
internal drive,

```math
\mu_\chi\ddot\chi+\gamma_\chi\dot\chi
+\omega_\chi^2(\chi-\alpha F)=\beta_{F2}F^2,
```

or in the readout,

```math
m/m_0=1+c_YF+c_\chi\chi+\lambda_{F\chi}F\chi+\lambda_{\chi2}\chi^2.
```

The generated lines are not important merely because they are new harmonics.
Static nonlinear comparators can also create harmonics. They become useful
only when they carry a shared denominator relation and when enough generated
samples are usable to beat the static nonlinear comparator budget.

For `F=(a/r)^p`, the generated coefficient is represented as

```math
B_k(p,e)=A_k(2p,e).
```

For a degree-`M` static nonlinear generated-line comparator and `K_Lambda`
shared projection nuisance parameters, the generated-line count must satisfy

```math
|U|\ge M+2+K_\Lambda.
```

In the default generated-line audit with `p=2`, `rho=3/2`, `delta=0.2`,
`H_gen=6`, and `M=1`, the moderate-eccentricity cases `e=0.3` and `e=0.6`
are generated-budget-breaking under the relative cutoff, while `e=0.1`
remains generated-underbudget. In the parameterized SNR examples at `e=0.3`,
the minimum generated-amplitude scales are `Q_beta=0.429654` for the
acceleration channel and `Q_beta=0.118468` for the range channel.

## 8. Component Separability, Robustness, And Analytic Boundaries

Generated nonlinear components can overlap the exact-in-e linear harmonics.
The observable model therefore has to separate them. The joint harmonic model
is

```math
\hat O_k=
\Lambda_k\left[
Q_YA_k(p,e)
+Q_LA_k(p,e)H_2(k)
+Q_\beta B_k(p,e)D_2(k)^{-1}
\right].
```

The generated component is locally separable only if the Jacobian column for
`Q_beta` adds rank relative to the model without `Q_beta`. A full design also
requires the shared-parameter Jacobian to have full rank. This rank condition
separates three cases: generated-component degeneracy, local separability
that is still budget-underdesigned, and component-separable budget readiness.

In the default audit, `e=0.3` is component-separable-and-budget-ready with
minimum cutoff `H=4` in the acceleration channel and `H=5` in the range
channels. The low-eccentricity case remains underdesigned or degenerate
depending on the projection nuisance model.

The robustness map then asks whether the positive branch is isolated. The
default finite grid varies `p`, `e`, `rho`, `delta`, projection channel, and
nuisance mode. It has `162` design points. The audit finds `107`
robust-positive points, `40` component-separable-but-budget-underdesigned
points, and `15` generated-component-degenerate points. This fraction is not
a physical population measure; it only shows that the positive branch is not
confined to a single tuned point in the chosen dimensionless grid.

The analytic phase-boundary theorem explains the grid with four boundaries:

```math
H_{\rm count}=M+2+K_\Lambda,
```

resonance bracketing with usable harmonics on both sides of `rho`,
amplitude usability by a relative cutoff or SNR analogue, and Jacobian rank
addition by `Q_beta`. The boundary audit reproduces all `162/162` default
robustness-map verdicts.

## 9. Discussion: What Is Shown And What Is Not Shown

The paper shows that the no-internal-state assumption has a concrete dynamic
loophole. The one-state model is not strong enough under the current
dictionaries, but it establishes the correct comparator discipline. The
nonlinear second-order branch then supplies a conditional positive result:
under exact-in-e eccentric forcing and finite shared projection nuisance,
there are parameterized regions with budget-breaking, locally separable
generated-line observables and analytic phase boundaries.

The paper does not show that an effect has been detected. It does not attach
the design variables to a real instrument or data set. It does not claim
identifiability under arbitrary per-frequency projection nuisance. It does
not treat sideband existence alone as a dynamic discriminator. It does not
interpret the robustness-grid fraction as a population probability.

The dynamic paper is therefore distinct from the static theorem and
classification paper. The static paper asks which instantaneous finite-family
collapse statements survive. The dynamic paper asks how such a theorem can
fail when A4 is dropped and an internal state becomes visible through shared
frequency-domain structure.

## 10. Conclusion

The current dynamic-loophole package closes as follows. First-order
relaxation is a control/no-go branch for the current forcing and projection
dictionaries. Nonlinear second-order internal visibility is the live
counterexample candidate. Exact-in-e eccentric forcing plus finite shared
projection nuisance can yield generated-line observables that are
budget-breaking, locally separable, robust over a finite dimensionless grid,
and governed by analytic phase boundaries.

This is the strongest honest claim supported by the repository. It is a
parameterized design theorem and a counterexample candidate to A4, not an
instrument forecast or a real-data detection claim.
