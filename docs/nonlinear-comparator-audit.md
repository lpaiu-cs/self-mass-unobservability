# Nonlinear Comparator Audit

Status: Counterexample candidate. This note asks whether the minimal nonlinear
`chi` sidebands are genuinely dynamic observables or only static nonlinear
response written in another basis.

## Models Being Compared

Status: Counterexample candidate. The dynamic nonlinear branch is

```math
\tau_\chi \dot\chi + \chi = \alpha F + \beta_{F2} F^2,
\qquad
q_{\rm dyn}=c_Y F + c_\chi\chi .
```

Status: Counterexample candidate. For a generated sideband at complex
frequency variable `z`, the nonlinear-drive piece has the transfer

```math
Q_{\rm dyn}(z)=\frac{C_{\rm dyn}}{1+\tau_\chi z},
\qquad
C_{\rm dyn}=c_\chi\beta_{F2}\,\widehat{F^2}(z).
```

Status: Conjectural. The static nonlinear comparator is any finite local
derivative ansatz of the schematic form

```math
q_{\rm stat}
=a_1F+a_2F^2+a_3F\dot F+a_4\dot F^2+\cdots
```

with a finite shared coefficient set across all sampled sidebands.

## Immediate Degeneracy

Status: Proven. Sideband existence alone is not enough. A static nonlinear
term `a_2F^2` can produce `2\Omega`, `\Omega_1+\Omega_2`,
`|\Omega_1-\Omega_2|`, and the orbital `3n` line whenever the required input
products are present.

Status: Proven. At a single generated sideband `z_*`, a free complex nonlinear
coefficient can match the dynamic amplitude exactly:

```math
a_{\rm eff}(z_*)=\frac{C_{\rm dyn}}{1+\tau_\chi z_*}.
```

Status: Proven. Therefore the M6 target is not "sideband exists." The target
is a shared `tau_chi` law across sidebands and the already-measured linear
lines.

## Shared-Pole Obstruction

Status: Proven. Let the static nonlinear comparator be a finite shared
polynomial in the sideband frequency variable,

```math
P_N(z)=\sum_{m=0}^N b_m z^m .
```

Status: Proven. It can interpolate the dynamic pole at `N+1` distinct complex
samples, but a new distinct sample `z_*` has residual

```math
\frac{C_{\rm dyn}(-\tau_\chi)^{N+1}
\prod_{j=0}^{N}(z_*-z_j)}
{(1+\tau_\chi z_*)\prod_{j=0}^{N}(1+\tau_\chi z_j)} .
```

Status: Proven. The residual vanishes only if `C_dyn=0`, `tau_chi=0`,
`z_*` duplicates a fitted node, or a pole/zero singular case is being sampled.
Otherwise the finite polynomial comparator cannot represent the relaxation
pole exactly beyond its interpolation budget.

## Ratio Targets

Status: Counterexample candidate. Define the linear transfer

```math
G(\omega)=c_Y+\frac{\beta}{1+i\omega\tau_\chi}.
```

Status: Counterexample candidate. For two-tone forcing, the deprojected
sum-sideband ratio target is

```math
\frac{\hat O(\Omega_1+\Omega_2)}
{\hat O(\Omega_1)\hat O(\Omega_2)}
=
\frac{C_{\rm dyn}}
{[1+i(\Omega_1+\Omega_2)\tau_\chi]G(\Omega_1)G(\Omega_2)} .
```

Status: Counterexample candidate. For the orbital forcing template with input
harmonics `n` and `2n`, the corresponding generated `3n` target is

```math
\frac{\hat O(3n)}{\hat O(n)\hat O(2n)}
=
\frac{C_{\rm dyn}}
{(1+3in\tau_\chi)G(n)G(2n)} .
```

Status: Counterexample candidate. These are stronger than sideband existence
because they require the same `tau_chi` to organize the fundamental lines and
the generated line.

## Verdict

Status: Proven. Minimal nonlinear dynamic `chi` sidebands defeat a purely
linear projection degeneracy because a linear time-invariant projection cannot
create frequencies absent from its input.

Status: Proven. The same sidebands do not defeat an unrestricted static
nonlinear comparator by existence alone.

Status: Counterexample candidate. The surviving M6 novelty candidate is the
shared pole law: sidebands plus linear lines must be fit by one common
`tau_chi`, while a finite shared static derivative polynomial fails once the
number of distinct generated samples exceeds its interpolation budget.

Status: Counterexample candidate. The promoted ratio version of this target is
recorded in [`shared-tau-ratio-theorem.md`](shared-tau-ratio-theorem.md).
