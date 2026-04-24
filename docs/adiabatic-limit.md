# Adiabatic Limit

Status: Proven. The first-order model can be written as

```text
(1 + tau_chi d/dt) chi = alpha F(t).
```

For slowly varying forcing, formally invert the operator:

```text
chi = alpha (1 + tau_chi d/dt)^(-1) F
    = alpha [F - tau_chi dot F + tau_chi^2 ddot F - tau_chi^3 F''' + ...].
```

Status: Proven. For `F(t)=F0 cos(Omega t)`, the periodic response is

```text
chi(t)
  = alpha F0 [cos(Omega t) + Omega tau_chi sin(Omega t)]
    / [1 + Omega^2 tau_chi^2].
```

Status: Proven. Expanding at small `x = Omega tau_chi`,

```text
chi(t)
  = alpha F0 [(1 - x^2 + O(x^4)) cos(Omega t)
              + (x - x^3 + O(x^5)) sin(Omega t)].
```

## Static Collapse Boundaries

Status: Proven. Exact instantaneous collapse occurs when `tau_chi = 0`:

```text
m/m0 = 1 + (c_Y + c_chi alpha) F.
```

Status: Proven. Exact zero-frequency collapse occurs when `Omega = 0`, because the quadrature term vanishes and the response is a static shift.

Status: Proven. If `Omega tau_chi << 1`, a local derivative EFT through order `N` can reproduce the response up to `O((Omega tau_chi)^(N+1))`.

Status: Proven. More sharply, for `z=i Omega` and `|Omega tau_chi| <= rho`,
the degree-`N` Taylor derivative residual is bounded by
`|alpha c_chi| rho^(N+1)`.

Status: Proven. The same operational collapse applies to the concrete
orbital-harmonic dictionary when all retained harmonics satisfy
`|k n tau_chi| <= rho` and the observable tolerance exceeds the Taylor
residual bound.

Status: Proven. In the range-like projection channel, the same adiabatic
collapse statement applies after deprojection when
`Gamma != 0` and `kappa^2 + z_k^2 != 0` on the sampled band.

Status: Proven. Exact reproduction of `alpha / (1 + i Omega tau_chi)` for arbitrary `Omega` requires either the explicit state `chi`, a hereditary kernel, or an infinite derivative tower.

## Boundary Statement

Status: Proven. The one-state relaxation model collapses to the static instantaneous sensitivity EFT only under `Omega tau_chi = 0`, `alpha c_chi = 0`, or when the measured sector is insensitive to quadrature.

Status: Proven. The same model collapses perturbatively to a finite local derivative EFT only after choosing a truncation order and restricting to `Omega tau_chi` below the corresponding error tolerance.

Status: Counterexample candidate. The theorem-breaking content is therefore not a static coefficient shift; it is the finite-relaxation transfer relation outside the adiabatic truncation boundary.
