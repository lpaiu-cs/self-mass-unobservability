# Request 10.8e: Phase-Marginalized Estimator And Physical Drive Anchoring

Status: Counterexample candidate. The two remaining publication items from
the 10.8 review: (1) the zero-drive-phase convention (review attack A3), and
(2) unit-drive normalization (A4). This note pre-registers both fixes,
committed before either scan is run on the real data.

## E1: phase-marginalized estimator (resolves A3)

Status: Counterexample candidate. The unresolved convention is exactly ONE
parameter: the integrator-time phase origin `t_off`. The carrier phases are
then determined, `phi_w = phi_w0 + w * t_off` (arbitrary INDEPENDENT
per-carrier phases are not marginalized — that is the Request 10.4
"arbitrary complex projection = collapse" class, physically excluded by the
10.5 phase-locked manifold). Estimator: at each `(tau, t_off)` the response
templates are the measured quadrature columns recombined,

```text
T_cY(t_off)     = sum_w [ cos(w t_off) col_w_c + sin(w t_off) col_w_s ]
T_beta(tau,t_off)= sum_w [ g1 cos + g2 sin ; rotated by w t_off ] (same columns)
```

and the linear `(c_Y, beta)` fit runs as in 10.8b. Grids: the 10.8b `tau`
grid; `t_off` on 24 uniform points in `[0, P_in)` (the fast carriers wrap
`2 pi` across it). Quantities, fixed in advance:

- Detection: `Z = max |z|` over the FULL `(tau, t_off)` grid, calibrated by
  500 simulations through the SAME max statistic (trials included); the
  anti-causal control (`g2 -> -g2`) scanned identically.
- Limit: `u95_pm(tau) = max over t_off of u95(tau, t_off)` — the
  phase-marginalized (worst-phase) limit, valid whatever the true phase.
  Quoted at the K_dyn = 10 anchor with the Fisher floor and K = 934 bracket
  alongside; this SUPERSEDES the zero-phase 10.8b/c curve as the headline
  (the zero-phase curve remains recorded).
- Gates: carry over per the 10.7e argument (the templates are linear
  recombinations of the SAME live-gate-validated measured columns; no new
  nonlinear-model content).

## E2: dictionary-anchored physical limit (resolves A4; Conjectural)

Status: Conjectural (model-anchored; drive-invariant choice per Request
10.2b). Drive = the Request-8 clock invariant `F = U/c^2` (the free-fall
response carries NO 1/Omega clock-integration factor). Leading-order
dimensionless drive amplitudes with O(1) geometric factors set to 1:

```text
f_in  = [G m_B/(a_in  c^2)] e_in   = 4.27e-11
f_out = [G m_C/(a_out c^2)] e_out  = 1.21e-10
f_dif = [G m_C beta_A eps/(a_out c^2)] = 1.12e-11
```

(masses/eccentricities/axes as in 10.7e E2). Phases tasc-locked
(`ph_in = -Omega_in tasc_p - pi/2`, `ph_out = -Omega_out tasc_b - pi/2`,
`ph_dif = ph_in - ph_out`) with the common `t_off` profiled as in E1.
Template: `a_w = beta_phys * f_w e^{i ph_w}/(1 + i w tau)` (+ the `c_Y f_w`
instantaneous term) mapped through the measured columns. Quote:
`u95_phys(tau) = max over t_off` at the K_dyn = 10 anchor, labeled
Conjectural, with the note that `beta_phys` is the dimensionless coupling of
`Delta` to `U/c^2` and the bound scales inversely with the true drive
factors.

## Outputs

Status: Note. `sep_dynamic/sep_phase_marg_10_8e.json`,
`sep_dynamic/sep_limit_curve_10_8e.tsv` (columns: zero-phase, phase-marg,
physical), updated `joint_fit_summary.md` and review note.
