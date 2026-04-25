# chi-state loophole

## Minimal Action

- Status: Counterexample candidate. Introduce one explicit orbital-timescale internal coordinate `\chi_A(\tau)` with inertial scale `\mu_\chi`:

```math
S_A
=
\int d\tau \left[
-m_0 \left(1 + s_1 Y(z_A) + \frac12 s_2 Y(z_A)^2\right)
+ \frac{\mu_\chi}{2} \dot{\chi}_A^2
- \frac{\mu_\chi \omega_\chi^2}{2} \chi_A^2
+ g_\chi \chi_A Y(z_A)
\right].
```

- Status: Counterexample candidate. The action remains local and analytic in the fields kept explicitly.
- Status: Counterexample candidate. The only directly violated theorem assumption is A4: no orbital-timescale internal state variable.

## Explicit Dynamics

- Status: Counterexample candidate. The Euler-Lagrange equation is

```math
\ddot{\chi}_A + \omega_\chi^2 \chi_A = \frac{g_\chi}{\mu_\chi} Y(z_A).
```

- Status: Counterexample candidate. The formal solution is

```math
\chi_A
=
\chi_{A,\mathrm{hom}}
+ \frac{g_\chi}{\mu_\chi}\left(\omega_\chi^2 + D_\tau^2\right)^{-1} Y.
```

- Status: Counterexample candidate. The homogeneous piece `\chi_{A,\mathrm{hom}}` carries independent initial data.
- Status: Counterexample candidate. The free-fall response therefore depends on more than the instantaneous local invariant `Y`; it depends on the extra state and its initial conditions.

## Why This Escapes Instantaneous Sensitivity Collapse

- Status: Counterexample candidate. When `\omega_\chi \sim \Omega_{\rm orb}`, the resolvent `(\omega_\chi^2 + D_\tau^2)^{-1}` is not reducible to a short local derivative expansion over the orbital timescale.
- Status: Counterexample candidate. In that regime the effective monopole response is not an algebraic function of finitely many instantaneous sensitivity coordinates.
- Status: Counterexample candidate. This is a genuine loophole to the theorem candidate if A4 is false, even though locality and analyticity remain intact.

## Adiabatic-Collapse Limit

- Status: Counterexample candidate. If `\omega_\chi \gg \Omega_{\rm orb}` and the homogeneous mode is absent or decays, then

```math
\chi_A
=
\frac{g_\chi}{\mu_\chi \omega_\chi^2}
\left[
1 - \frac{D_\tau^2}{\omega_\chi^2}
+ O\!\left(\frac{D_\tau^4}{\omega_\chi^4}\right)
\right] Y.
```

- Status: Counterexample candidate. Substituting back gives a local expansion of the form

```math
\Delta L_{\rm eff}
=
\frac{g_\chi^2}{2\mu_\chi \omega_\chi^2} Y^2
+ \frac{g_\chi^2}{2\mu_\chi \omega_\chi^4} (\dot Y)^2
+ O\!\left(\omega_\chi^{-6}\right),
```

up to total derivatives.

- Status: Counterexample candidate. In this heavy-state limit the loophole collapses back into ordinary local sensitivity and Wilson-coefficient data.
- Status: Counterexample candidate. The `chi` loophole is therefore minimal in a precise sense: it fails only when the extra state survives on the orbital timescale.

## Explicit Boundary Between Collapse And Non-Collapse

- Status: Counterexample candidate. Define the adiabatic parameter

```math
\varepsilon_\chi := \frac{\Omega_{\rm orb}}{\omega_\chi}.
```

- Status: Counterexample candidate. `adiabatic-collapse` requires both `\varepsilon_\chi \ll 1` and absence or decay of the homogeneous mode `\chi_{A,\mathrm{hom}}`.
- Status: Counterexample candidate. `non-collapse` occurs whenever either `\varepsilon_\chi = O(1)` or larger, or the homogeneous mode survives with independent initial data.
- Status: Counterexample candidate. The theorem-relevant boundary is therefore not vague heaviness; it is the explicit failure or success of `\varepsilon_\chi \ll 1` plus homogeneous-mode suppression.

## Smallest-Model Reason

- Status: Counterexample candidate. Only one new degree of freedom is added, and the rest of the worldline EFT remains local and analytic.

## Nonlinear Sideband Comparator Boundary

- Status: Proven. A nonlinear extension of the `chi` response can create sidebands, but sideband existence alone is not a unique dynamic signature once static nonlinear local comparators are admitted.
- Status: Counterexample candidate. The stronger observable target is a shared relaxation pole across generated sidebands and the linear-frequency response.
- Status: Proven. The current comparator boundary is recorded in [`../../docs/nonlinear-comparator-audit.md`](../../docs/nonlinear-comparator-audit.md) and checked in [`../../symbolic/nonlinear_comparator_audit.py`](../../symbolic/nonlinear_comparator_audit.py).
- Status: Counterexample candidate. The current ratio theorem/no-go target is recorded in [`../../docs/shared-tau-ratio-theorem.md`](../../docs/shared-tau-ratio-theorem.md) and checked in [`../../symbolic/shared_tau_ratio_audit.py`](../../symbolic/shared_tau_ratio_audit.py).
- Status: Proven. The current minimum-sample boundary is recorded in [`../../docs/sample-budget-theorem.md`](../../docs/sample-budget-theorem.md) and checked in [`../../symbolic/sample_budget_audit.py`](../../symbolic/sample_budget_audit.py); the present orbital and two-tone dictionaries are underbudget for realistic first-derivative linear comparators.
- Status: Counterexample candidate. The bounded richer-orbital rescue route is recorded in [`../../docs/richer-orbital-forcing-budget-audit.md`](../../docs/richer-orbital-forcing-budget-audit.md) and checked in [`../../symbolic/orbital_harmonic_budget_audit.py`](../../symbolic/orbital_harmonic_budget_audit.py).
- Status: Counterexample candidate. The stronger second-order internal mode route is recorded in [`../../docs/second-order-internal-mode.md`](../../docs/second-order-internal-mode.md) and checked in [`../../symbolic/second_order_mode_response.py`](../../symbolic/second_order_mode_response.py).
- Status: Counterexample candidate. The second-order projection and resonant-comparator budgets are recorded in [`../../docs/second-order-projection-audit.md`](../../docs/second-order-projection-audit.md), [`../../docs/resonant-comparator-theorem.md`](../../docs/resonant-comparator-theorem.md), [`../../symbolic/second_order_projection_audit.py`](../../symbolic/second_order_projection_audit.py), and [`../../symbolic/resonant_comparator_audit.py`](../../symbolic/resonant_comparator_audit.py).
- Status: Counterexample candidate. The exact-in-e second-order forcing sample-provisioning audit is recorded in [`../../docs/exact-in-e-resonant-forcing-budget-audit.md`](../../docs/exact-in-e-resonant-forcing-budget-audit.md) and checked in [`../../symbolic/exact_in_e_resonant_forcing_audit.py`](../../symbolic/exact_in_e_resonant_forcing_audit.py).
- Status: Counterexample candidate. The amplitude-weighted second-order design theorem is recorded in [`../../docs/amplitude-weighted-resonant-design-theorem.md`](../../docs/amplitude-weighted-resonant-design-theorem.md) and checked in [`../../symbolic/amplitude_weighted_resonant_design.py`](../../symbolic/amplitude_weighted_resonant_design.py).
- Status: Counterexample candidate. The physical detectability map is recorded in [`../../docs/physical-detectability-map.md`](../../docs/physical-detectability-map.md) and checked in [`../../symbolic/physical_detectability_map.py`](../../symbolic/physical_detectability_map.py); it replaces the relative cutoff by an explicit SNR budget and minimum channel-scale boundary.
- Status: Counterexample candidate. The minimal nonlinear second-order branch is recorded in [`../../docs/nonlinear-second-order-mode.md`](../../docs/nonlinear-second-order-mode.md) and checked in [`../../symbolic/nonlinear_second_order_mode.py`](../../symbolic/nonlinear_second_order_mode.py); it asks whether generated lines and linear lines share the same quadratic denominator beyond static nonlinear comparator budgets.
- Status: Counterexample candidate. The nonlinear second-order generated-line detectability audit is recorded in [`../../docs/nonlinear-second-order-detectability.md`](../../docs/nonlinear-second-order-detectability.md) and checked in [`../../symbolic/nonlinear_second_order_detectability.py`](../../symbolic/nonlinear_second_order_detectability.py); it classifies whether exact-in-e `F^2` harmonics beat the static nonlinear budget under relative or SNR usable-sample criteria.
- Status: Counterexample candidate. The generated-component separability theorem is recorded in [`../../docs/component-separability-theorem.md`](../../docs/component-separability-theorem.md) and checked in [`../../symbolic/component_separability_audit.py`](../../symbolic/component_separability_audit.py); it tests whether the nonlinear amplitude `Q_beta` is locally identifiable when linear and generated terms overlap at the same harmonics.
