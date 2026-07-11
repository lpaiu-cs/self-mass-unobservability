# Request 10.8c: Resolving The K = 934 Anchor Factor

Status: Counterexample candidate. The 10.8b headline carries K = 934 = 
(published static marginal 1.8e-6)/(our Fisher static 1.922e-9) as a
conservative multiplier of unknown origin. This note pre-registers a
data-independent diagnosis and a re-anchoring rule, BEFORE the diagnostic
ratio spectrum is computed.

Status: Note. Provenance correction recorded here: `steps.npz` carries
`frac = 1.0`, i.e. the 10.7a Jacobian steps were 1.0 x sigma_MCMC (the
`jacobian_runner.py` docstring said 0.3x; the stored frac is authoritative).
`sigma_MCMC,j := abs_steps_j` (from the released chain; the chain file itself
is not in the runtime tree, so steps.npz is the source of record,
`sep_dynamic/steps_provenance.json`). The 10.8 secant critique is unaffected
(tasc_extra1 step = 8739 d = 1.0 sigma_MCMC = 2.8 planet periods).

## Diagnostic (pre-registered before computation)

Status: Counterexample candidate. Hypotheses for the 934x gap:

- H-global: our Fisher machinery (noise model / linear span) is globally
  optimistic -> the gap should appear for ALL well-measured parameters.
- H-direction: the gap is specific to nonlinearly-degenerate directions
  (the static Delta trades against masses/periods through turn-wrap-scale
  excursions; the published posterior width reflects that nonlinear manifold,
  which our point Fisher cannot see). The 10.8b live gates PROVED this
  mechanism does not operate on the dynamic templates (the real model could
  not absorb them at any tested amplitude, up to ~48,000 whitened units).

Discriminator: the ratio spectrum `rho_j = sigma_MCMC,j / sigma_Fisher,j`
over the CORE parameter set (the 20 fitted parameters excluding the 7
`*_extra1` planet parameters and `delta_i`, which are known
nonlinear/wrap-dominated). `sigma_Fisher,j` = full-rank leave-one-out
projection marginal (`s / ||P_perp(others) a_j||`) in the 10.8b nuisance
context (all other Jacobian columns + offset + 30 Fourier pairs + static SEP
column); NO SVD truncation (the comparison target is their full posterior).

## Re-Anchoring Rule (fixed in advance)

- `median(rho_core) <= 3` AND `max(rho_core) <= 10` -> H-direction confirmed:
  the machinery is globally calibrated and the static gap is
  direction-specific. Re-anchor `K_dyn := max(10, ceil(max(rho_core)))`
  (safety floor 10) and quote the re-anchored 10.8b curve.
- `3 < median(rho_core) <= 30` -> partial global gap:
  `K_dyn := 3 x max(rho_core)`, quote re-anchored curve.
- `median(rho_core) > 30` -> H-global unresolved: K stays 934; no re-quote.
- In all branches the planet-block and SEP ratios are reported as mechanism
  diagnostics, and the Fisher floor remains reported alongside.

## Outputs

Status: Note. `sep_dynamic/anchor_resolution_10_8c.json`, re-anchored curve
columns appended to `sep_dynamic/sep_beta_limit_curve_10_8b.tsv` companion
(`sep_beta_limit_curve_10_8c.tsv`), updated `joint_fit_summary.md`, plus a
referee-style review note of the 10.8 chain.
