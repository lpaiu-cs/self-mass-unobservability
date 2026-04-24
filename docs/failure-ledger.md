# Failure Ledger

Status: Imported from prior work. This ledger replaces static primitive-family bookkeeping as the main failure surface for this worktree.

## Current Minimal Model

Status: Counterexample candidate. The tested model is

```text
tau_chi * d chi / dt + chi = alpha F(t)
m / m0 = 1 + c_Y F(t) + c_chi chi(t)
```

with `F(t) = F0 cos(Omega t)` for the first derivation.

## Exact Collapse Conditions

| Status | Condition | Exact failing step | Classification |
| --- | --- | --- | --- |
| Proven | `tau_chi = 0` | The equation becomes `chi = alpha F`, so `m/m0 = 1 + (c_Y + c_chi alpha)F`. | Static sensitivity redefinition. |
| Proven | `Omega = 0` | The quadrature coefficient `alpha F0 Omega tau_chi / (1 + Omega^2 tau_chi^2)` vanishes. | Static or secular shift only. |
| Proven | `alpha c_chi = 0` | The internal state either is not driven or is not read out in the mass response. | No dynamic observable. |
| Proven | `Omega tau_chi << 1` with a local derivative EFT truncated at order `N` | `chi = alpha sum_{n=0}^N (-tau_chi d/dt)^n F + O((Omega tau_chi)^{N+1})`. | Order-by-order derivative collapse. |
| Proven | Single known drive frequency with an unconstrained local `dot F` Wilson coefficient | The leading quadrature can be fit by one derivative coefficient without proving an internal state. | Static-basis-degenerate signal. |
| Proven | Single known drive frequency with unconstrained `{F, dot F}` coefficients | Both cosine and sine quadratures are fit by `a0 F + a1 dot F`, but the fitted `a1` depends on `Omega`. | Single-frequency derivative degeneracy. |
| Proven | Linear first-order `chi` plus linear two-frequency forcing and linear readout | Superposition leaves only the input frequencies; no `Omega_1 +/- Omega_2` terms appear. | No sideband loophole in the linear MVP. |

## Non-Collapse Candidate

Status: Counterexample candidate. At finite `Omega tau_chi`, the exact transfer

```text
H(Omega) = alpha / (1 + i Omega tau_chi)
```

has a relaxation pole and cannot be represented exactly by a finite instantaneous static sensitivity basis.

Status: Counterexample candidate. If the allowed comparator is a finite local derivative EFT, the loophole is not a single quadrature point; it is the frequency-dependent transfer relation across drives or the need for an infinite derivative tower to reproduce the pole exactly.

## Failed Sideband Attempt

Status: Proven. The linear two-frequency MVP does not produce sidebands.

Exact failing step:

```text
F(t) = F1 cos(Omega1 t) + F2 cos(Omega2 t)
```

implies

```text
chi(t) = chi_1(t) + chi_2(t)
```

by linearity, so the readout `c_Y F + c_chi chi` contains only `Omega1` and `Omega2`.

Minimal missing assumption for sidebands:

Status: Conjectural. Add a nonlinear drive or readout, such as `m/m0` containing `c_chi2 chi^2`, or a nonlinear `F(Y)` evaluated on a multi-frequency invariant.

## Next Escalations

Status: Counterexample candidate. If the one-state model is judged insufficient against a derivative-EFT comparator, the next smallest escalations are:

1. hereditary kernel with non-rational memory response,
2. second-order internal mode with resonant phase behavior,
3. nonlinear readout or nonlinear drive to generate sidebands,
4. nonanalytic threshold or hysteresis to break local analytic expansion.
