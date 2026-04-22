# Request 6 Source Scout

## Goal

Find actual `gamma`-measured binary-pulsar systems that can extend the decoupled
clock-only Request 6 analysis in a meaningful way, not just by adding more
systems but by adding real `s_p` lever arm around the current local basis

`s* = 0.134196`.

The scout is intentionally tied to the current Request 6 machinery:

- the same surrogate `s_p(m_p)` prior from `request6_clock_sector.py`
- the same linearized `eta_* , kappa_*` basis around `s*`
- the same Fisher-style `kappa_*` audit from `request6_lever_arm_audit.py`

This is therefore a **source-selection audit for the current staging analysis**,
not an exhaustive catalogue of all relativistic pulsars.

## What The Audit Actually Wants

From `request6_lever_arm_audit`, the practical single-source target for opening
the slope direction is roughly

- `|Delta s| ~ 0.054`
- `sigma_delta ~ 2e-4`

where `sigma_delta` is the fractional uncertainty of the compressed clock
observable `delta_gamma = eta(s)/(1+X_c)`.

This target is useful because it lets us interpret real sources cleanly:

- a source can have the right `Delta s` but useless precision
- or excellent precision but too little `Delta s`
- or sit on the wrong side of `s*`

## Main Findings

### 1. The strongest published single addition in the scanned set is `PSR B1913+16`

`PSR B1913+16` is not the largest-lever-arm source, but it is by far the most
precise one checked.

- `Delta s = +0.0258`
- `sigma_delta = 1.86e-4`
- adding it alone moves `|kappa_*|_95` from `8.42` to `5.20e-2`

So the real published sample already contains a source with the **precision**
needed by the audit, but not with the ideal `|Delta s| ~ 0.054` offset.

### 2. `PSR J1913+1102` gives the right high-side lever arm, but not the precision

This is the best high-`s` source found in the scout.

- `Delta s = +0.0583`
- `sigma_delta = 3.18e-2`
- adding it alone only gets to `|kappa_*|_95 = 1.53`

So `J1913+1102` is the right **directional** source, but it is not yet a slope
driver by itself in the current PK-summary approximation.

### 3. The usable low-side candidate is `PSR J1141-6545`, but only conditionally

Among the sources checked, `J1141-6545` is the most useful opposite-side source
that is not immediately disqualified.

- `Delta s = -0.0154`
- `sigma_delta = 1.42e-2`
- adding it alone gets to `|kappa_*|_95 = 2.44`

This is not enough numerically, and the source is only conditionally usable
because the WD-driven `xdot` / spin-orbit structure must be modelled explicitly.

### 4. `PSR J1906+0746` is interesting but not front-line

Numerically, `J1906+0746` is the best negative-side source in the raw scout
because its `gamma` precision is better than `J1141-6545`.

- `Delta s = -0.0052`
- `sigma_delta = 4.36e-3`
- adding it alone gets to `|kappa_*|_95 = 2.23`

But the current 2026 timing paper explicitly reports that fitting `xdot` shifts
the inferred masses by about `3.5 sigma` because `xdot` and `gamma` are
correlated. For Request 6 this makes it a **tricky follow-up target**, not the
first clean source to add.

### 5. Clean Shapiro systems are methodologically nice, but they do not add real `s` leverage

Both `PSR J1757-1854` and `PSR J1756-2251` are attractive because their masses
can be reconstructed from non-`gamma` PK information with Shapiro support.

But for the current clock-sector problem, they sit too close to the existing
`s*` direction:

- `PSR J1757-1854`: `Delta s = +0.0004`, `|kappa_*|_95 = 8.41`
- `PSR J1756-2251`: `Delta s = +0.0011`, `|kappa_*|_95 = 8.16`

So they are good **consistency systems**, not good **lever-arm systems**.

## Ranked Candidate Summary

| System | Role in Request 6 | `Delta s` | `sigma_delta` | Single-source `|kappa_*|_95` | Status |
| --- | --- | ---: | ---: | ---: | --- |
| `PSR B1913+16` | best precision anchor | `+0.0258` | `1.86e-4` | `5.20e-2` | moderate |
| `PSR J1913+1102` | best high-side lever arm | `+0.0583` | `3.18e-2` | `1.53` | moderate |
| `PSR J1141-6545` | best usable low-side source | `-0.0154` | `1.42e-2` | `2.44` | conditional |
| `PSR J1906+0746` | precise but correlated low-side source | `-0.0052` | `4.36e-3` | `2.23` | tricky |
| `PSR J1756-2251` | clean Shapiro cross-check | `+0.0011` | `7.84e-3` | `8.16` | clean |
| `PSR J1757-1854` | clean Shapiro cross-check | `+0.0004` | `3.35e-3` | `8.41` | clean |

## Pair-Level Result

The best pair among the scanned published sources is still dominated by
`PSR B1913+16`.

- best pair, unrestricted:
  - `PSR B1913+16 + PSR J1757-1854`
  - `|kappa_*|_95 = 5.17e-2`
- best cross-side pair:
  - `PSR B1913+16 + PSR J1906+0746`
  - `|kappa_*|_95 = 5.17e-2`
- best cross-side non-tricky pair:
  - `PSR B1913+16 + PSR J1141-6545`
  - `|kappa_*|_95 = 5.20e-2`

This means the scanned published set still misses the audit target badly:

- target: `|kappa_*|_95 < 1e-2`
- best actual non-tricky cross-side pair: `5.20e-2`

So the problem is now very clear:

- the **positive/high-side** is covered by `B1913+16` for precision and
  `J1913+1102` for lever arm
- the **negative/low-side** remains the bottleneck

## Practical Next Step

For the next Request 6 upgrade, the most defensible order is:

1. Add `PSR B1913+16` first in a covariance-aware PK fit.
2. Add one low-side system with explicit nuisance handling.
   - If prioritizing physical cleanliness: `PSR J1141-6545`
   - If prioritizing raw timing precision: `PSR J1906+0746`, but only with
     explicit `xdot-gamma` covariance treatment
3. Keep `PSR J1757-1854` and `PSR J1756-2251` as clean cross-check systems,
   not as primary slope drivers.

## Source Notes

- `PSR J1756-2251` and `PSR J1141-6545` source tables label `gamma` in `ms`
  while giving numerical values on the `10^-3 s` scale. The scout interprets
  those entries in seconds because that is the only interpretation consistent
  with the GR Einstein-delay scale and the rest of the papers.
- `PSR J1141-6545` uses the 2008 timing paper for the `gamma` proxy and the
  2020 `xdot` / spin-orbit paper for the current caveat and mass scale.

## Primary Sources Checked

- `PSR B1913+16`: [Weisberg, Nice, Taylor 2010](https://arxiv.org/abs/1011.0718)
- `PSR J1757-1854`: [Cameron et al. 2017](https://arxiv.org/abs/1711.07697)
- `PSR J1756-2251`: [Ferdman et al. 2014](https://arxiv.org/abs/1406.5507)
- `PSR J1913+1102`: [Ferdman et al. 2020](https://arxiv.org/abs/2007.04175)
- `PSR J1141-6545`: [Bhat, Bailes, Verbiest 2008](https://arxiv.org/abs/0804.0956), [Venkatraman Krishnan et al. 2020](https://arxiv.org/abs/2001.11405)
- `PSR J1906+0746`: [van Leeuwen et al. 2026](https://arxiv.org/abs/2602.05947)
