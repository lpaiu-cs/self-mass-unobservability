# Request 10.8b: Real-Data Fit Of The Dynamic Free-Fall SEP Channel

Status: Counterexample candidate. The 10.8 feasibility gates passed (F2
survival up to 43.8%; F3 anchored projection 6.1e-8 <= 1e-5), promoting the
channel to a real-data experiment. This note pre-registers that experiment.

Status: Note. Pre-registration discipline: committed BEFORE any real-data
`beta_ff` statistic is computed. Only design-matrix quantities (survivals,
Fisher/anchored sigmas) and the 10.7-chain real-data statistics are known.

## Model And Templates

Status: Counterexample candidate. Signal: the shared one-state pole transfer
applied to the SEP parameter in the free-fall sector,

```text
Delta(t) = Delta_static + Re sum_w [ (c_Y + beta_ff/(1 + i w tau_chi)) e^{i w t} ],
w in {Omega_in, Omega_out, Omega_dif},  unit drive per carrier, zero phase.
```

Templates from the T2 MEASURED exact response columns (no kernel assumption:
the columns are the integrator's response):

```text
T_beta(tau) = sum_w [ g1 col_w_c + g2 col_w_s ],   g1 = 1/(1+w^2 tau^2),
T_cY        = sum_w col_w_c,                        g2 = w tau/(1+w^2 tau^2).
```

## Pre-Registered Estimator

Status: Counterexample candidate. Joint weighted fit, FWL block elimination:

- Nuisance block: 28 normalized Jacobian columns (v2) + offset + 30 Fourier
  pairs + the NORMALIZED STATIC SEP column `j_Delta` (static-leakage guard);
  SVD truncation at relative cut 1e-3. Noise scale `s = sqrt(chi2_red)` of the
  nuisance-only fit; stop if `s > 1.75`.
- Signal block: `[T_cY, T_beta(tau)]`; `tau` grid: 61 log points on
  `[2, 500] d` plus anchors `{2, 5, 18, 52, 200} d`.
- Detection: `Z = max|z|` over the quoted window; 500-simulation null
  (seed 20260710); anti-causal control (`g2 -> -g2`, an unphysical
  advanced-response template) scanned identically with its own null — a
  detection candidate requires global `p < 0.003` AND a quiet anti-causal
  control, else upper limits.
- `u95(tau)`: flat-prior Gaussian posterior, quoted in TWO bracket curves
  fixed in advance: HEADLINE = anchored (`K = 934` times the Fisher sigma;
  conservative, carrying the unexplained static Fisher-vs-published gap) and
  FLOOR = raw Fisher (statistical-only). The anchored curve is THE quoted
  limit; resolving K is deferred to 10.8c.
- Linear injection-recovery at all anchors at the anchored u95 amplitudes
  (coverage <= 0.5 sigma_Fisher... reported in Fisher units).

## Pre-Registered Live Gates (D-gate analogue; WSL, real model)

Status: Counterexample candidate.

- G1 (live null): one GN fit of the real residuals (truncated-pinv updates on
  the 28-parameter block, <= 4 iterations, 0.5% break); the 10.8b estimator
  measured on the remainder must track the linear pipeline:
  `|z_GN - z_lin| <= 0.3` at anchors `{18, 52} d`.
- G2 (estimator calibration): inject `beta_inj = u95_anchored(tau)` of
  `T_beta(tau)` into the real residuals at anchors `{18, 52} d`, GN-absorb
  with the live model, recover through the 10.8b estimator:
  `|beta_hat_GN - beta_hat_lin_null - beta_inj| <= 0.5 sigma_Fisher(tau)`.

## Quote Rules (fixed in advance)

- No detection + G1/G2 pass + anti-causal control quiet -> QUOTE the anchored
  u95 curve over `[2, 500] d` as the 10.8b upper limit (Fisher floor
  reported alongside).
- Detection candidate -> suspend quoting; run the detection protocol
  (systematics review, anti-causal and null-sim scrutiny) before any claim.
- Any gate failure or control anomaly -> no quote; 10.8c required.

## Outputs

Status: Note. `sep_dynamic/sep_joint_fit_10_8b.json`,
`sep_dynamic/sep_beta_limit_curve_10_8b.tsv`, `sep_dynamic/sep_gateG.json`,
updated `joint_fit_summary.md`.
