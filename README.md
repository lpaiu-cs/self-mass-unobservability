# SMU Clock Timing Dictionary

This repository derives the timing-observable dictionary for the decoupled
clock sector of the self-mass-unobservability EFT.

The working scientific question is:

> Does the decoupled clock EFT predict anything beyond Einstein delay and
> trivial secular spin-redshift renormalization?

## Scope

This repo is theory/observability-first. It does not build or revive external
timing runtimes.

Assumptions for the MVP:

- free-fall sector = GR
- propagation sector = GR
- clock sector = decoupled EFT only
- binary pulsar timing is closed before any hierarchical-triple extension

The clock EFT used here is

```text
d tau_A / dt
= 1 - v_A^2/(2 c^2)
  - (1 + zeta_1 s_A + zeta_2 s_A^2) U_ext / c^2
```

where `s_A` is a sensitivity-like body parameter.

## MVP Deliverables

The binary MVP is limited to:

- one derivation memo
- one symbolic derivation script
- one observable-classification table
- one synthetic sanity-check script
- one theorem-or-counterexample memo

This initial pass creates the first three scaffolding pieces only. It does not
state a final theorem, does not inspect triple systems, and does not use
external estimators.

## Current Entry Point

Generate the first binary dictionary table:

```bash
python scripts/run_binary_dictionary.py
```

The script writes:

- `outputs/tsv/binary_pk_dictionary.tsv`
- `outputs/json/binary_pk_dictionary.json`
- `outputs/tsv/binary_pk_dictionary_exact_e.tsv`
- `outputs/json/binary_pk_dictionary_exact_e.json`

The table schema is:

```text
term | expansion_order | harmonic | observable_status | absorption_target | notes
```

## Repository Layout

```text
docs/
  problem_statement.md
  notation.md
notes/
  REQUEST8_CLOCK_BINARY_DICTIONARY.md
src/
  smu_clock/
    binary_pk_dictionary.py
scripts/
  run_binary_dictionary.py
outputs/
  json/
  tsv/
```

Existing historical request folders are treated as background context only.
They are not source material for this Request 8 MVP unless explicitly promoted
into the new binary dictionary notes.
