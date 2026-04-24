# AGENTS.md

## Repository Scope

This repository derives the timing-observable dictionary for the decoupled
clock sector of the self-mass-unobservability EFT.

The active scientific question is:

```text
Does the decoupled clock EFT predict anything beyond Einstein delay and
trivial secular spin-redshift renormalization?
```

Assume unless a later note explicitly changes it:

- free-fall sector = GR
- propagation sector = GR
- clock sector = decoupled EFT only
- binary pulsar timing must close before any hierarchical-triple extension

## Active EFT

Use the clock-rate map:

```text
d tau_A / dt
= 1 - v_A^2/(2 c^2)
  - (1 + zeta_1 s_A + zeta_2 s_A^2) U_ext / c^2
```

For the binary MVP:

```text
r = a (1 - e cos E)
U_ext = G m_B / r
v_A^2 = (m_B / M)^2 G M (2 / r - 1 / a)
```

## Required Discipline

For every derivation, state:

- expansion variable,
- post-Keplerian order,
- eccentricity order,
- whether each term is periodic or secular,
- whether each term is distinguishable from standard timing nuisance
  parameters.

Every projected term must be classified as exactly one of:

- `observable`
- `absorbed`
- `degenerate`
- `null`

Do not leave a dictionary row unclassified.

## Stop Rules

If the binary dictionary closes with only:

- Einstein-delay renormalization, plus
- secular spin-redshift renormalization,

then stop expanding the binary branch and write that as the main binary result.

If a later triple extension also collapses into standard nuisance structure,
write a no-go result rather than continuing source-driven exploration.

## Non-Goals

Do not:

- start triple-system work before the binary deliverables exist,
- source-chase additional gamma systems,
- build or provision external estimators,
- run PEP, MLRS, Tempo2, Nutimo, or legacy timing toolchains,
- use Request 3 scale checks as evidence,
- claim a final tied-vs-decoupled verdict from this repo alone,
- state a final gamma-only theorem before the binary derivation is closed.

## Preferred Artifacts

Use small, reversible edits and keep derivation separate from numerics.

Preferred paths:

- `notes/REQUEST8_*.md`
- `scripts/run_*.py`
- `outputs/json/*.json`
- `outputs/tsv/*.tsv`
- `outputs/figures/*.svg`

The current MVP is binary-only.
