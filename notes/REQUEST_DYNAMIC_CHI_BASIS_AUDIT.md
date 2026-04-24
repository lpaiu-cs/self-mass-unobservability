# Dynamic Chi Basis Audit

## Goal

Status: Counterexample candidate. The MVP found a finite-`Omega tau_chi`
phase-lag response in the one-state model. This audit asks whether that signal
survives against progressively broader comparator bases.

## Model

Status: Proven. The model is

```text
tau_chi dot chi + chi = alpha F(t)
delta(m/m0) = c_Y F(t) + c_chi chi(t)
```

with the monochromatic drive

```text
F(t) = F0 cos(Omega t).
```

Status: Proven. The readout transfer is

```text
G(Omega)
= c_Y + c_chi alpha/(1 + i Omega tau_chi).
```

Equivalently, the periodic readout contains

```text
F0 [c_Y + c_chi alpha/(1 + Omega^2 tau_chi^2)] cos(Omega t)
+ F0 [c_chi alpha Omega tau_chi/(1 + Omega^2 tau_chi^2)] sin(Omega t).
```

## Comparator Bases

Status: Proven. Against an instantaneous static `F` coefficient, the in-phase
component is degenerate and the quadrature component remains outside the basis.

Status: Proven. Against a one-frequency local derivative basis

```text
a0 F + a1 dot F,
```

one can choose

```text
a0 = c_Y + c_chi alpha/(1 + Omega^2 tau_chi^2),
a1 = -c_chi alpha tau_chi/(1 + Omega^2 tau_chi^2),
```

so the single-frequency quadrature is degenerate with a frequency-local
`dot F` Wilson coefficient.

Status: Proven. In the adiabatic regime, the dynamic state is reproduced
order-by-order by

```text
chi = alpha [F - tau_chi dot F + tau_chi^2 ddot F - ...].
```

Status: Proven. Exact reproduction of

```text
alpha/(1 + i Omega tau_chi)
```

over varying `Omega` requires an infinite derivative tower, the explicit
state `chi`, or an equivalent hereditary kernel. No finite polynomial in
`i Omega` equals this rational transfer when `alpha != 0` and
`tau_chi != 0`.

## Audit Table

| term | comparator | verdict |
| --- | --- | --- |
| static in-phase amplitude | instantaneous `F` | degenerate |
| static quadrature phase lag | instantaneous `F` | observable candidate |
| single-frequency quadrature | `{F, dot F}` | degenerate |
| adiabatic response | finite derivative EFT at chosen order | absorbed through truncation |
| finite-tau transfer relation | finite shared-coefficient polynomial derivative EFT | observable candidate |
| `tau_chi=0` | static EFT | null |
| `Omega=0` | static/secular EFT | null |
| `alpha c_chi=0` | any basis | null |
| linear two-frequency sidebands | sideband basis | null |

## Interpretation

Status: Counterexample candidate. The MVP does not yet establish novelty from
a single quadrature point if the comparator allows a free `dot F` coefficient.

Status: Counterexample candidate. The stronger surviving target is the
frequency-dependent rational transfer relation with common `tau_chi`, because
that relation is outside a finite static or finite-order polynomial derivative
basis.

Status: Proven. The linear two-frequency model still does not generate
sidebands. Sideband novelty requires a nonlinear drive, nonlinear readout,
or a higher-order internal mode.

## Machine Outputs

Status: Proven. The audit writes:

```text
outputs/tsv/dynamic_chi_basis_audit.tsv
outputs/json/dynamic_chi_basis_audit.json
```
