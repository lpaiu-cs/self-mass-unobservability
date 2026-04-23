# Request 4: APOLLO Baseline Nuisance Scaffold

## Scope

This stage still does **not** perform a GR residual fit.

It implements the narrow nuisance/design scaffold that should sit between
`real-data ingest` and any actual APOLLO-only baseline fit. The goal is to
make the next estimator target explicit and reproducible without pretending
that a production LLR residual model already exists.

## Why This Stage Exists

After the APOLLO ingest step, the next real bottleneck is no longer data
access. It is estimator design.

The current scaffold is built to encode three practical choices:

- APOLLO-only is a useful first weak-field branch, but still a single-station
  branch,
- batch boundaries may coincide with reduction-version boundaries and therefore
  need explicit nuisance treatment,
- the first baseline target should be narrow and controlled, not a fake
  all-purpose "production LLR estimator".

## Implementation

The code is in
[request4_llr_apollo_baseline_scaffold.py](request4_llr_apollo_baseline_scaffold.py#L1).

It reads the parsed APOLLO normal-point table and writes a basis/design table
with:

- batch indicators:
  `group_b`, `group_c_replacement`, `group_d_complete`,
  `group_e_complete`, `group_f`
- reflector indicators:
  `Apollo11`, `Apollo14`, `Luna17`, `Luna21`
  with `Apollo15` as the reference reflector
- centered time terms:
  linear and quadratic in centered tropical years
- meteo / observational surrogates:
  `pressure_z`, `temperature_z`, `humidity_z`, `np_duration_z`,
  `log_photon_count_z`, `uncertainty_z`
- harmonic surrogate blocks:
  annual, synodic-month, sidereal-month, and anomalistic-month `sin/cos`

These are not yet a physical Earth-Moon model. They are only the nuisance side
of the first narrow baseline target.

## Outputs

- basis table:
  [request4_llr_apollo_baseline_basis.tsv](request4_llr_apollo_baseline_basis.tsv#L1)
- scaffold summary:
  [request4_llr_apollo_baseline_scaffold_summary.json](request4_llr_apollo_baseline_scaffold_summary.json#L1)
- scaffold summary figure:
  [request4_llr_apollo_baseline_scaffold_summary.svg](request4_llr_apollo_baseline_scaffold_summary.svg)

## Current Role In Request 4

This scaffold is the right intermediate layer if the next APOLLO-only target is
a narrow baseline model with:

- Earth-Moon initial-state freedom on the physics side,
- reflector-specific bias terms,
- batch-level nuisance terms,
- a small meteo/geometry surrogate block,
- and a small lunar-periodic surrogate block.

That is deliberately smaller than a real production LLR estimator.

## What It Does Not Do

This scaffold does **not** yet provide:

- a GR prediction for the normal points,
- a residual series,
- a fit to physical nuisance parameters,
- or any EFT parameter inference.

So it should be read as a design/input layer only.

## Stop Rule

The stop rule for the self-built APOLLO-only branch should now be explicit.

If a narrow APOLLO-only baseline GR residual model cannot be stabilized on top
of this scaffold, then the correct next move is **not** to keep growing a
bespoke estimator indefinitely. The correct next move is either:

1. attach to an existing LLR codebase / estimation stack, or
2. pivot to the CRD/ILRS canonical data path.

That is the practical boundary between a healthy scaffold and a scope trap.

That follow-on probe has now been carried out in
[REQUEST4_LLR_BASELINE_FIT.md](REQUEST4_LLR_BASELINE_FIT.md#L1),
and the current result falls on the `pivot` side of that boundary rather than
the `keep growing the bespoke surrogate` side.
