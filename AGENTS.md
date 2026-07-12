# AGENTS.md

## Mission

This repository is dynamic-loophole-first.

Status: Imported from prior work. The earlier theorem work left A4 as the active assumption to attack: no orbital-timescale internal state variable in the free-fall sector.

Status: Conjectural. This repository asks whether the smallest such state, `chi_A`, creates a real observable beyond a static finite-dimensional sensitivity-manifold EFT.

## Primary Scope

Work on the free-fall-style response model first.
Do not start with the clock sector.
Do not reopen LLR / MLRS / PEP / Nutimo / runtime / build-environment work.
Do not make static primitive-family audits the mainline unless they are strictly needed to define the drive basis `Y`.

## Core Target

Status: Counterexample candidate. The first target is

```text
tau_chi * d chi_A / dt + chi_A = alpha * F(Y)
m_A(Y, chi_A) = m_A^(0) * [1 + c_Y F(Y) + c_chi chi_A]
```

Status: Conjectural. The target observable must be one of:

- phase-lagged quadrature relative to a static drive,
- sidebands or mixed-frequency response,
- frequency-dependent transfer not absorbable into static coefficients,
- a sharply stated no-go boundary showing collapse of the minimal model.

## Allowed Outputs

Codex may produce:

- theorem or no-go statements,
- assumption ledgers,
- short analytic derivations,
- symbolic response checks,
- explicit counterexample candidates,
- failure ledgers with exact collapse conditions.

Codex must not:

- do more empirical LLR / MLRS / PEP work,
- do more pulsar / Nutimo runtime work,
- invent empirical claims without derivation or source,
- silently strengthen assumptions,
- count a static parameter redefinition as novelty.

## Claim Labels

Every substantive scientific claim must be tagged as exactly one of:

- Proven
- Imported from prior work
- Conjectural
- Counterexample candidate

## Maintained Files

Maintain these files during dynamic-chi work:

- `docs/model-definition.md`
- `docs/observable-targets.md`
- `docs/adiabatic-limit.md`
- `docs/nonadiabatic-regime.md`
- `docs/failure-ledger-dynamic-chi.md`

If a model fails to escape collapse, record the exact failing step and the minimal missing assumption in `docs/failure-ledger-dynamic-chi.md` (the theorem track owns `docs/failure-ledger.md`).

## MVP Done Rule

A task is done only if:

- the relevant markdown notes are updated,
- symbolic checks run without error,
- and the result is classified as theorem progress or loophole progress.

For this MVP, done means:

1. the one-state `chi_A` model is defined,
2. the monochromatic response is solved,
3. the adiabatic collapse boundary is written,
4. the non-adiabatic observable classification is written,
5. either a genuine observable candidate is isolated or the exact no-go boundary is stated.

## Git Discipline

Before each major task:

- create a checkpoint commit.

After each major task:

- create a checkpoint commit with a one-line scientific summary.
