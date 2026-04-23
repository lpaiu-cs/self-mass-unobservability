# Request 6: Low-Side Covariance-Proxy Push

## Goal

This stage executes the one additional low-side push that was still worth
doing inside the currently accessible data level.

It goes beyond the earlier simple surrogate stage by carrying some of the
published timing-analysis covariance structure directly into the low-side
likelihood, without claiming a full TOA-level fit.

The practical question is:

- if the low-side branch is upgraded from a simple PK-summary surrogate to a
  timing-summary covariance proxy,
- does `|kappa_*|_95` finally move toward the `10^-2` regime,
- or does Request 6 remain a local-audit branch only?

The answer is that it remains a local-audit branch only.

## What Changed Relative To The Previous Low-Side Stage

The previous low-side stage treated both `J1141` and `J1906` with a single
shared `gamma` observable per source and simple inclination branches.

This new stage makes two stronger changes.

### `PSR J1141-6545`

`J1141` is upgraded from a single scintillation-based inclination branch to an
analysis-branch mixture:

- `scintillation_independent`: `i = 76 ± 2.5 deg`
- `timing_xdot_geometry`: `i = 71 ± 2 deg`

The second branch is kept only as a nuisance-aware timing branch, not as a
clean independent source. Both branches use the same published `gamma`
measurement from the 2008 timing paper.

### `PSR J1906+0746`

`J1906` is the more important upgrade.

Instead of using one shared `gamma = 4.59(2) × 10^-4 s` and then letting that
measurement kill the `xdot` branch, this stage carries the published
`xdot-gamma` covariance at the timing-summary level by using two explicit
analysis branches:

- `dd_external_geometry`: `i = 45 ± 3 deg`, `gamma = 4.59(2) × 10^-4 s`
- `ddgr_xdot_correlated`: `i = 50 ± 1 deg`, `gamma = 4.1(1) × 10^-4 s`

This is still not a full TOA-level covariance model, but it is a much more
honest approximation to the published statement that fitting `xdot` shifts
`gamma` and the individual masses together.

## Primary Sources

- `PSR J1141-6545`:
  [Bhat, Bailes, Verbiest 2008](https://arxiv.org/abs/0804.0956),
  [Venkatraman Krishnan et al. 2020](https://arxiv.org/abs/2001.11405)
- `PSR J1906+0746`:
  [van Leeuwen et al. 2026](https://arxiv.org/abs/2602.05947)

## Main Numerical Result

Running

```bash
python3 request6_low_side_covariance_proxy.py
```

gives, at the original Request 6 basis point `s* = 0.134196`:

- baseline after `B1913`: `|kappa_*|_95 = 4.877e-2`
- previous simple low-side stage (`J1141 + J1906`): `4.861e-2`
- `+ J1141` covariance proxy: `4.875e-2`
- `+ J1906` covariance proxy: `4.877e-2`
- `+ both` covariance proxies: `4.875e-2`

So the stronger covariance-proxy push does **not** improve the slope audit.
Relative to the previous simple low-side stage, it slightly worsens it.

This is the main result.

## Source-Level Diagnostics

### `PSR J1141-6545`

Effective row:

- `sbar = 0.119334`
- `X_c = 0.443877`
- `sigma_delta = 1.940e-2`
- `delta_gamma_obs = -3.94e-5`

GR-side branch weights:

- `scintillation_independent = 0.536`
- `timing_xdot_geometry = 0.464`

Interpretation:

- `J1141` remains a genuine low-side source,
- but even after branch-aware timing treatment it is still far too broad to
  move `kappa_*` in a meaningful way.

### `PSR J1906+0746`

Effective row:

- `sbar = 0.153588`
- `X_c = 0.459032`
- `sigma_delta = 4.130e-2`
- `delta_gamma_obs = 7.97e-4`

GR-side branch weights:

- `dd_external_geometry = 0.072`
- `ddgr_xdot_correlated = 0.928`

Interpretation:

- once the published `xdot-gamma` covariance is carried through, the
  `ddgr_xdot_correlated` branch dominates,
- that pushes the pulsar mass up to about `1.41 M_sun`,
- and `J1906` ceases to behave like a useful low-side source at all.

So the stronger treatment does not rescue the low-side lever arm. It instead
shows more clearly that `J1906` is an analysis-systematic mixture between a
mildly low-side and a positive-side timing solution.

## Combined Interpretation

This stage sharpens the project-level decision.

At the covariance-proxy level:

- `J1141` stays low-side but weak,
- `J1906` stops being cleanly low-side once the published covariance structure
  is respected,
- and the combined low-side update leaves `|kappa_*|_95` essentially unchanged
  from the `B1913` baseline.

Around the combined weighted reference point, the clock-only branch still
constrains mostly the local amplitude direction:

- `s_ref = 0.155864`
- `|eta_*|_95 = 2.086e-3`
- `|kappa_*|_95 = 4.622e-2`

So even after the stronger low-side push, Request 6 still does **not** become a
competitive slope-measurement engine.

## Programmatic Consequence

This result is strong enough to change the project status.

At the published-summary level, the one additional covariance-aware low-side
push has already been tried and it failed to move the branch toward the
`10^-2` regime. That means:

- the current Request 6 ceiling is now numerically clear,
- additional source-chasing at the same level is low-return,
- and only a genuine TOA-level low-side timing fit could still materially
  change the verdict.

Absent that TOA-level escalation, the practical stop rule should now be treated
as met:

- Request 6 gets frozen as a support / local clock-sector audit,
- and the tied-vs-decoupled bottom line belongs to the joint
  free-fall-plus-clock analysis.

## Files

- [request6_low_side_covariance_proxy.py](request6_low_side_covariance_proxy.py#L1)
- [request6_low_side_covariance_proxy_summary.json](request6_low_side_covariance_proxy_summary.json#L1)
- [request6_low_side_covariance_proxy_summary.svg](request6_low_side_covariance_proxy_summary.svg)
- [request6_low_side_covariance_proxy_posterior_j1141.tsv](request6_low_side_covariance_proxy_posterior_j1141.tsv#L1)
- [request6_low_side_covariance_proxy_posterior_j1906.tsv](request6_low_side_covariance_proxy_posterior_j1906.tsv#L1)
- [request6_low_side_covariance_proxy_posterior_both.tsv](request6_low_side_covariance_proxy_posterior_both.tsv#L1)
