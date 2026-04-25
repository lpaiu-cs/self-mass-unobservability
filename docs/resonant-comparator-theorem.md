# Resonant Comparator Theorem

Status: Counterexample candidate. This note states the second-order analogue
of the shared-tau sample-budget theorem.

## Target

Status: Counterexample candidate. The target observable is the calibrated
line shape

```math
H_2(z)=
\frac{\alpha\omega_\chi^2}
{\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2}.
```

Unlike the one-state pole, this transfer can carry a resonance peak and phase
wrapping.

## Static Comparator Budget

Status: Conjectural. The conservative static line-shape comparator is a finite
shared local derivative polynomial

```math
P_N(z)=\sum_{j=0}^{N}a_j z^j
```

plus `K_Lambda` shared projection nuisance parameters.

Status: Proven. A finite sample claim is underbudget unless

```math
N_\Omega > N+1+K_\Lambda .
```

Thus the minimum complex frequency-sample count is

```math
N_\Omega^{\min}=N+2+K_\Lambda .
```

Status: Proven. For the second-order range channel, `K_Lambda=2` for
`Gamma,kappa^2` if they are not separately calibrated; for the
acceleration-like channel, `K_Lambda=1` for `Gamma`.

## Resonance-Aware Design Rule

Status: Counterexample candidate. A resonance-aware design should also bracket
the phase-wrap pole

```math
\Omega^2=\omega_\chi^2/\mu_\chi,
```

and the amplitude peak exists only when

```math
\gamma_\chi^2<2\mu_\chi\omega_\chi^2.
```

Status: Proven. Bracketing or observing a sign flip is not by itself a novelty
claim: an underbudget finite polynomial comparator can interpolate finite
phase samples.

Status: Counterexample candidate. The useful design condition is therefore
both:

```math
N_\Omega \ge N+2+K_\Lambda
```

and frequency placement on both sides of the phase-wrap denominator.

## Example Budgets

| Comparator/projection class | `K_Lambda` | Minimum samples | Status |
| --- | ---: | ---: | --- |
| `N=1`, calibrated projection | `0` | `3` | Proven |
| `N=1`, acceleration `Gamma` nuisance | `1` | `4` | Proven |
| `N=1`, range `Gamma,kappa^2` nuisance | `2` | `5` | Proven |
| `N=2`, range `Gamma,kappa^2` nuisance | `2` | `6` | Proven |

## Collapse Boundary

Status: Proven. If `mu_chi=0`, the model collapses to the first-order
relaxation transfer.

Status: Proven. If the sampled band is far below the internal scale, the
transfer has a local derivative expansion and collapses order-by-order into
static coefficients.

Status: Proven. If projection nuisance is arbitrary at every frequency, the
line-shape theorem cannot be established by a finite sample design.

Status: Counterexample candidate. The second-order program is stronger than
the one-state program only if it samples the resonance/phase-wrap structure
with enough points to exceed the finite shared static comparator budget.

