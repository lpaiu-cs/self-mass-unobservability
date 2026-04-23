# Request 6: Low-Side Source Follow-Up (`J1141-6545` and `J1906+0746`)

## Goal

This stage executes the next two source additions identified by the scout,
after the `B1913+16` covariance-aware extension was already in place.

The question is narrow:

- if both available low-side follow-ups are added,
- does the decoupled clock-only Request 6 program finally open the local slope
  direction `kappa_*`,
- or does the low-side bottleneck remain?

The answer, in the current staging framework, is that the bottleneck remains.

## What Was Implemented

The code path is in `request6_low_side_extensions.py`.

Baseline:

- current Request 6 clock-only leave-one-out fit,
- plus the nuisance-aware `B1913+16` extension from
  `request6_b1913_covariance.py`.

Two low-side follow-ups were then added.

### `PSR J1141-6545`

This is treated as the cleaner but less precise low-side source.

- non-`gamma` mass inference uses the independent scintillation-geometry
  inclination prior,
- the WD-driven `xdot` / spin-orbit sector is left explicit in the metadata as
  a nuisance caveat rather than folded into the likelihood,
- the implementation is therefore useful and honest, but still conditional.

### `PSR J1906+0746`

This is treated as the more precise but more nuisance-sensitive source.

- one branch uses the external inclination quoted in the timing paper,
- one branch uses the `DDGR + xdot` inclination branch,
- the latter is explicitly marked as correlated with `gamma`, because the paper
  states that fitting `xdot` shifts the masses by about `3.5 sigma`.

So this source is implemented as a branch-mixed surrogate, not as a full
covariance-closed timing fit.

## Primary Inputs

- `PSR J1141-6545`:
  [Bhat, Bailes, Verbiest 2008](https://arxiv.org/abs/0804.0956),
  [Venkatraman Krishnan et al. 2020](https://arxiv.org/abs/2001.11405)
- `PSR J1906+0746`:
  [van Leeuwen et al. 2026](https://arxiv.org/abs/2602.05947)

## Main Numerical Result

Running

```bash
python3 request6_low_side_extensions.py
```

currently gives, at the original Request 6 local basis point
`s* = 0.134196`:

- baseline after `B1913`: `|kappa_*|_95 = 4.877e-2`
- `+ J1141`: `4.874e-2`
- `+ J1906`: `4.864e-2`
- `+ both`: `4.861e-2`

So the combined low-side follow-up only improves the current staged
`kappa_*` bound by about `0.3%` relative to the `B1913` baseline.

That is the main result of this stage.

## Source-Level Diagnostics

### `PSR J1141-6545`

Effective row:

- `sbar = 0.120135`
- `X_c = 0.441328`
- `sigma_delta = 1.787e-2`
- `delta_gamma_obs = 7.520e-3`

GR-side branch-selected metadata:

- `m_p = 1.278868 +- 0.008243 M_sun`
- `m_c = 1.010255 +- 0.008251 M_sun`
- `i = 74.799 +- 1.734 deg`

Interpretation:

- this is a real low-side source in the sense that its effective `sbar` sits
  below the current local pivot,
- but its compressed observable is too broad to move the slope audit much.

### `PSR J1906+0746`

Effective row:

- `sbar = 0.128679`
- `X_c = 0.496944`
- `sigma_delta = 6.163e-3`
- `delta_gamma_obs = 2.813e-4`

GR-side branch-selected metadata:

- `m_p = 1.314632 +- 0.004248 M_sun`
- `m_c = 1.298658 +- 0.004248 M_sun`
- `i = 44.693 +- 0.185 deg`

Branch weights at GR:

- `external_i = 1.000000`
- `xdot_conditioned_i = 4.96e-21`

Interpretation:

- nominally this source is more precise than `J1141`,
- but the `gamma` likelihood almost completely suppresses the
  `xdot`-conditioned branch,
- so the source does not become the decisive low-side lever-arm system that
  the broader Request 6 program still lacks.

## Combined Interpretation

This stage settles something important about the current Request 6 program.

It does **not** show that the clock-only branch is uninformative. The branch is
already informative about the local amplitude direction. After combining the
baseline with both low-side follow-ups, the posterior around its own weighted
reference point is

- `s_ref = 0.155832`
- `|eta_*|_95 = 2.106e-3`
- `|kappa_*|_95 = 4.628e-2`

But it **does** show that the slope direction remains effectively closed at the
present PK-summary / surrogate-`s_p` level.

So the correct reading is:

- `B1913` remains the precision anchor on the positive side,
- `J1141` remains the physically cleaner low-side source,
- `J1906` remains the more precise but nuisance-entangled low-side source,
- and even taking both low-side sources together does not materially change the
  Request 6 ceiling.

That means Request 6 is still a **local clock-sector audit**, not the arena in
which the tied-vs-decoupled question is finally decided.

## Programmatic Consequence

The next meaningful upgrade is no longer “add another published `gamma`
system.” The current stage already shows that this will not solve the core
problem by itself.

The next meaningful upgrade has to be one of:

1. a genuinely covariance-aware low-side PK likelihood,
2. a TOA-level low-side timing fit with the nuisance sectors carried directly,
3. and that next upgrade should be treated as the decision point for the
   Request 6 branch,
4. if even that stronger low-side treatment still does not move
   `|kappa_*|_95` into the `10^-2` range, then further effort here is
   low-return and Request 6 should be explicitly frozen as a
   support/local-audit section,
5. the main tied-vs-decoupled verdict should then be taken from the joint
   free-fall-plus-clock program rather than from additional clock-only source
   chasing.

## Branch Decision Rule

This stage makes the branch logic concrete.

- `B1913 + J1141 + J1906` already show that adding more published
  PK-summary-level `gamma` systems is not enough.
- So only one more class of upgrade is worth pursuing here:
  covariance-aware or TOA-level low-side treatment.
- If that still fails to pull `|kappa_*|_95` into the `10^-2` regime, then the
  correct move is to stop pushing Request 6 as a main-engine novelty claim.

At that point, Request 6 should be kept in the paper or note, but only as a
documented local clock-sector audit. The project's tied-vs-decoupled
bottom-line should come from the joint free-fall-plus-clock consistency
analysis.

## Files

- [request6_low_side_extensions.py](request6_low_side_extensions.py#L1)
- [request6_low_side_extensions_summary.json](request6_low_side_extensions_summary.json#L1)
- [request6_low_side_extensions_summary.svg](request6_low_side_extensions_summary.svg)
- [request6_low_side_posterior_j1141.tsv](request6_low_side_posterior_j1141.tsv#L1)
- [request6_low_side_posterior_j1906.tsv](request6_low_side_posterior_j1906.tsv#L1)
- [request6_low_side_posterior_both.tsv](request6_low_side_posterior_both.tsv#L1)
