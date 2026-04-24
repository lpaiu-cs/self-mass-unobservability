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
| Proven | `K <= floor((N+1)/2)` positive sampled frequencies with freely fitted real degree-`N` derivative coefficients | Even and odd channels can be interpolated at all sampled points. | Physical finite-sample degeneracy. |
| Proven | `K <= N + 1` sampled frequencies with freely fitted complex degree-`N` derivative coefficients | Polynomial interpolation can match the relaxation transfer at all sampled points. | Conservative complex finite-sample degeneracy. |
| Proven | Low-frequency sweep with `|Omega tau_chi| <= rho` and tolerance above `|alpha c_chi| rho^(N+1)` | Taylor truncation through order `N` is within tolerance over the band. | Operational derivative-EFT collapse. |
| Proven | Concrete orbital template with `e = 0`, `p = 0`, or `p = -3` at `O(e^2)` | The two-harmonic `n, 2n` sweep is absent or incomplete. | Insufficient forcing dictionary. |
| Proven | Observable projection with arbitrary complex `Lambda(Omega_k)` at each frequency | The projection nuisance can absorb `G(i Omega_k)` point by point. | Observable dictionary collapse. |
| Proven | Acceleration-like projection with `Gamma = 0` | `delta a_hat = Gamma q_hat` is blind to the body response. | Projection-channel collapse. |
| Proven | Range-like projection sampled at `kappa^2 + z_k^2 = 0` without a resonance model | Deprojection by `(kappa^2+z^2)/Gamma` is singular. | Projection-channel singularity. |
| Proven | Projection with a zero at `z=-1/tau_chi` | The observable projection cancels the relaxation pole. | Pole-cancellation degeneracy. |
| Proven | Triple bridge with `Omega_in = Omega_out` or `F_in F_out = 0` | The inner/outer carrier inventory has fewer than two usable frequency samples. | One-frequency derivative degeneracy. |
| Proven | Triple bridge against real degree `N >= 3` derivative EFT with only two carriers | The real odd channel has enough coefficients to interpolate both positive-frequency samples. | Need additional carrier or derivative-order prior. |
| Proven | Triple bridge against complex degree `N >= 1` polynomial comparator with only two carriers | Complex polynomial interpolation absorbs two samples. | Conservative complex-comparator degeneracy. |
| Proven | Three-carrier triple inventory with `Omega_in = 2 Omega_out` or `Omega_out = 2 Omega_in` | The combination carrier `|Omega_in-Omega_out|` duplicates one monopole carrier. | Falls back to two-carrier bridge. |
| Proven | Three-carrier triple inventory against real degree `N >= 5` or complex degree `N >= 2` comparator | Three samples are below the corresponding obstruction count. | Need sideband carrier, order prior, or projection prior. |
| Proven | Three-carrier triple inventory with arbitrary complex projection nuisance per carrier | The projection absorbs the shared transfer law point by point. | Observable bridge collapse. |
| Proven | Triple projection gate with arbitrary complex `Lambda_k` independently assigned to each carrier | Choosing `Lambda_k=O_k/(G(i Omega_k)F_k)` fits every carrier pointwise. | Projection-nuisance collapse. |
| Proven | Triple projection gate with zero, singular, or pole-cancelling projection | A carrier sample is blind, cannot be deprojected, or cancels the relaxation pole. | Carrier sample lost. |
| Proven | Three-carrier gate with a projection/comparator class whose effective dimension exceeds the three-sample pressure | The finite carrier inventory is below the obstruction count after nuisance parameters are admitted. | Need more carriers, order prior, or external projection prior. |
| Proven | Linear first-order `chi` plus linear two-frequency forcing and linear readout | Superposition leaves only the input frequencies; no `Omega_1 +/- Omega_2` terms appear. | No sideband loophole in the linear MVP. |
| Proven | Nonlinear sideband model with `beta_F2 = lambda_Fchi = lambda_chi2 = 0` | The model reduces to the linear no-sideband one-state response. | Nonlinear sideband collapse. |
| Proven | Mixed sideband test with `F1 F2 = 0` | Sum and difference sideband amplitudes vanish because only one harmonic is present. | Insufficient nonlinear forcing. |
| Proven | Orbital `3n` sideband test with `e = 0`, `p = 0`, `p = -3`, or `beta_F2 = 0` | The `n,2n` input pair or nonlinear drive coefficient vanishes. | Orbital sideband collapse. |

## Non-Collapse Candidate

Status: Counterexample candidate. At finite `Omega tau_chi`, the exact transfer

```text
H(Omega) = alpha / (1 + i Omega tau_chi)
```

has a relaxation pole and cannot be represented exactly by a finite instantaneous static sensitivity basis.

Status: Counterexample candidate. If the allowed comparator is a finite local derivative EFT, the loophole is not a single quadrature point; it is the frequency-dependent transfer relation across drives or the need for an infinite derivative tower to reproduce the pole exactly.

Status: Counterexample candidate. For a fixed derivative order `N`,
`floor((N+1)/2)+1` distinct positive frequencies give the first exact
obstruction to real shared-coefficient absorption. The obstruction is zero only
when `alpha c_chi = 0`, `tau_chi = 0`, or the frequencies are not distinct.

Status: Counterexample candidate. If complex derivative coefficients are
allowed, `N + 2` distinct frequencies give the conservative exact obstruction.
The obstruction is zero only when `alpha c_chi = 0`, `tau_chi = 0`, the
frequencies are not distinct, or a sampled point lies on the pole.

Status: Counterexample candidate. The current concrete forcing candidate is
`F(Y(t))=(a/r(t))^p` on a small-eccentricity orbit. It is a structural
observable dictionary, not an empirical detectability claim.

Status: Counterexample candidate. The hierarchical-triple shared-tau bridge is
the current clean carrier-inventory route: existing inner/outer carriers can
test real degree-`1` and degree-`2` shared derivative comparators without
requiring a new tensor harmonic.

Status: Counterexample candidate. Adding the existing GR outer-dipole
combination carrier `|Omega_in-Omega_out|` strengthens the same route to real
degree `N <= 4` and complex degree `N <= 1`, subject to nonresonance and finite
projection assumptions.

Status: Counterexample candidate. The current projection-channel audit shows
that a range-like linear readout has finite shared nuisance parameters
`Gamma,kappa`, not arbitrary frequency-local nuisance, but a named measurement
channel still needs its own projection justification.

Status: Counterexample candidate. The triple projection-nuisance realism gate
keeps the dynamic-chi branch alive only when the carrier projection is
calibrated or finite-dimensional and shared across carriers.

Status: Proven. The same gate collapses if the standard timing nuisance model
is allowed arbitrary per-carrier complex projection amplitudes.

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

Status: Proven. Adding `beta_F2 F^2` to the internal drive or
`lambda_Fchi F chi + lambda_chi2 chi^2` to the readout is sufficient to
generate sidebands under two-frequency forcing.

## Next Escalations

Status: Counterexample candidate. If the one-state model is judged insufficient against a derivative-EFT comparator, the next smallest escalations are:

1. hereditary kernel with non-rational memory response,
2. second-order internal mode with resonant phase behavior,
3. projection-specific sideband survival checks for named measurement channels,
4. nonanalytic threshold or hysteresis to break local analytic expansion.
