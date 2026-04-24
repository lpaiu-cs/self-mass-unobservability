# Problem Statement

## Question

Does the decoupled clock EFT predict a binary pulsar timing observable beyond
Einstein delay and trivial secular spin-redshift renormalization?

## Sector Split

The MVP assumes:

- free-fall sector: GR
- propagation sector: GR
- clock sector: decoupled EFT

The only modified ingredient is the clock-rate map

```text
d tau_A / dt
= 1 - v_A^2/(2 c^2)
  - (1 + zeta_1 s_A + zeta_2 s_A^2) U_ext / c^2 .
```

No orbital equations of motion, light-propagation model, or external timing
estimator is modified in this repo.

## Binary MVP

The binary branch must answer, at minimum:

- which clock-sector terms enter the timing model at leading post-Keplerian
  order,
- which binary harmonics they project onto,
- whether each term is observable, absorbed, degenerate, or null,
- which nuisance or standard timing parameter absorbs a non-observable term,
- whether a gamma-only theorem is supported or whether a counterexample exists.

## Expansion Discipline

Every term must state:

- expansion variable,
- post-Keplerian order,
- eccentricity order,
- periodic or secular character,
- distinguishability from standard timing nuisance parameters.

The first skeleton uses a Newtonian Kepler ellipse and expands the clock-rate
terms at `O(c^-2)` through the declared eccentricity order. It is a derivation
scaffold, not a final theorem.

## Stop Rule

If the binary dictionary closes with only Einstein-delay renormalization and
secular spin-redshift renormalization, the binary branch should stop and write
that no-go result rather than expanding into a source-driven search.
