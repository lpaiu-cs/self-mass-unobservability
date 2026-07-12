# Request 10.7e: Physical Anchoring And Red-Noise Marginalization Of The 10.7d Limit

Status: Counterexample candidate. The 10.7d quoted limit is in unit-drive
convention (`Lambda_k F_k = 1`) under a white-noise model (EFAC 1.1225 +
global scale). This note pre-registers (E1) convention-free per-carrier
pole-amplitude bounds, (E2) a physically anchored dimensionless bound via the
Request 8 triple outer-potential clock dictionary, and (E3) a red-noise
marginalized re-derivation. Committed BEFORE any of the new pipeline modes is
run on the real data.

## E1 (pre-registered arithmetic): per-carrier pole-amplitude bounds

Status: Counterexample candidate. The quoted `u95(tau)` bounds the common
pre-transfer drive amplitude of the pole term. The implied 95% bound on the
OBSERVED timing amplitude of the pole term at carrier `k` is

```text
A_k^pole,95(tau) = u95(tau) / |1 + i Omega_k tau|    [us]
```

reported at the anchors for `k in {in, out, dif}`. Pure arithmetic on the
already-quoted curve; no new data look.

## E2 (pre-registered model anchor): Request-8 dictionary drive template

Status: Conjectural (model-anchored; the drive-invariant choice is the
Request 10.2b `Conjectural` step). Drive = the leading clock-sector
invariant of `notes/REQUEST8_TRIPLE_OUTER_POTENTIAL.md` (branch
`clock-timing-dictionary`): `F = U/c^2` acting through
`d tau/dt = ... - alpha_A U/c^2`, with the dynamic-chi transfer applied to
the fractional-rate readout. Leading-order carrier drive amplitudes, with all
O(1) geometric projection factors set to 1 (documented limitation):

```text
D_in  = [G m_B/(a_in  c^2)] * e_in  / Omega_in[rad/s] * 1e6   [us]
D_out = [G m_C/(a_out c^2)] * e_out / Omega_out       * 1e6   [us]
D_c   = [G m_C beta_A eps/(a_out c^2)]  / Omega_c     * 1e6   [us]
beta_A = m_B/(m_A+m_B),  eps = a_in/a_out,
```

with `m_A = 1.43781441, m_B = 0.19753639, m_C = 0.41010271 Msun`,
`P_in = 1.62939901 d, P_out = 327.25512703 d`,
`e_in = |(eta_p, kappa_p)| = 6.990e-4, e_out = |(eta_b, kappa_b)| = 3.529e-2`
(all from the baseline parfile / Nutimo construction), `a` from Kepler.
Numerically `D ~= (0.96, 545, 0.25) us` (in, out, dif) — the outer carrier
dominates. Phases are locked to the orbit: `ph_in = -Omega_in tasc_p - pi/2`,
`ph_out = -Omega_out tasc_b - pi/2`, `ph_c = ph_in - ph_out` (the 10.5
phase-locked manifold), `tasc_p = -575.13824 d`, `tasc_b = -262.10186 d`.

Status: Counterexample candidate. Estimator: the SAME pre-registered joint
pipeline (truncated span, seed, grids, u95 definition, detection rule,
controls) with the signal columns re-weighted:
`A_k(cY) = D_k e^{i ph_k}` and `A_k(beta) = D_k e^{i ph_k}/(1 + i Omega_k tau)`.
The fitted `beta_phys` is then the DIMENSIONLESS amplitude of the dynamic
response of `alpha_A` (the pole part of the sensitivity transfer), and
`u95_phys(tau)` is its 95% bound. D1/D2 validity carries over: the estimator
change is a re-weighting of exactly-linear signal columns; the nonlinear-model
content (the 28-parameter block) is unchanged and was gate-validated at
amplitudes bracketing the physical template's carrier weights.

## E3 (pre-registered): red-noise marginalization

Status: Counterexample candidate. Nuisance extension: `M = 30` Fourier pairs
`{cos, sin}(2 pi j t / T)`, `j = 1..30`, `T = 2987.8605 d` (the TOA span),
appended to the 29-column nuisance block before the SVD truncation (same
relative cut 1e-3). This marginalizes, without a spectral prior, all
low-frequency structure up to `f = 30/T ~= 1/(100 d)`; the pair nearest the
outer carrier (`j = 9`, `2 pi 9/T = 0.01893` vs `Omega_out = 0.01920` rad/d,
phase separation 0.81 rad over the span) is deliberately INCLUDED — the joint
fit charges the honest price. Anchors at small tau (inner/dif carriers,
`f ~ 1834/T`) are far above the marginalized band.

Status: Note. Gate carry-over argument: the added columns are exact sinusoids
(linear by construction); D1/D2 validated the only nonlinear ingredient (the
28-parameter model response). No new WSL gate is required. A periodogram of
the projected residuals is reported for context.

## Runs And Quote Rules

Status: Counterexample candidate. Four deterministic pipeline runs (same
seed): `10_7d` (already quoted: unit drive, white), `_rn` (unit drive, red-
marginalized), `_phys` (dictionary template, white), `_phys_rn` (dictionary
template, red-marginalized). Rules, fixed in advance:

- Detection is re-evaluated per run with the 10.7b rule (global `p < 0.003`
  AND controls quiet); any detection candidate suspends quoting.
- E1 bounds are quoted from the already-promoted `10_7d` curve.
- The red-marginalized unit-drive curve `u95_rn` REPLACES the 10.7d curve as
  the headline quoted limit if it differs by more than 10% anywhere on the
  quoted window (noise-model robustness dominates); otherwise the 10.7d curve
  stands with `u95_rn` reported as a robustness band.
- `u95_phys_rn` is quoted as the physical anchor (labeled Conjectural,
  model-anchored, O(1) geometric factors), with `u95_phys` as its white-noise
  cross-check.

## Outputs

Status: Note. `joint_fit_upper_limit_rn.json`, `_phys.json`, `_phys_rn.json`
(+ matching `beta_limit_curve_*.tsv`), `pole_amplitude_bounds.json`,
`residual_periodogram.json`, updated `joint_fit_summary.md`, under
`../request10_external/`.
