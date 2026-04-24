# AGENTS.md

## Mission

This repository pursues a theorem-level follow-up:

Prove or refute that self-mass-unobservability hypotheses collapse to a
finite-dimensional sensitivity-manifold EFT under explicit assumptions.

If the theorem fails, the required output is not hand-waving.
It is the smallest explicit loophole model and the exact assumption it violates.

## Primary scope

Work on the free-fall sector first.
Do NOT start with the clock sector.
Do NOT reopen LLR / MLRS / PEP / Nutimo / build-environment work in this repo.

## Core theorem target

Target theorem:

Given a quasi-static, nearly spherical, nonspinning compact body with
self-bound equilibrium, a local worldline EFT, no orbital-timescale internal
state variable, and analytic coupling to external fields, the leading
body-dependent deviation must collapse to sensitivity coordinates
(s_A,I, s_A,IJ, ...) plus higher-multipole Wilson coefficients.

The 1D scalar-s EFT is only a corollary if manifold reduction is separately justified.

## Allowed outputs

Codex may produce:
- theorem statements
- assumption ledgers
- proof skeletons
- symbolic derivations
- explicit counterexample candidates
- comparison notes to the current project artifacts

Codex must NOT:
- do more real-data LLR work
- do more pulsar runtime/build work
- invent empirical claims without a derivation or source
- silently strengthen assumptions

## Working rules

Every nontrivial claim must be labeled as one of:
- Proven
- Imported from prior work
- Conjectural
- Counterexample candidate

Maintain these files at all times:
- docs/assumptions-ledger.md
- docs/failure-ledger.md
- docs/roadmap.md

If a proof attempt fails, record the exact failing step and the minimal missing assumption.
Do not patch over it with prose.

## First milestone

Milestone M1 is complete only when all of the following exist:
1. Theorem A draft for free-fall sensitivity collapse
2. Lemma note for internal-structure no-go
3. Lemma note for COM decoupling
4. Generic worldline action note
5. One explicit loophole candidate

## Done-when rule

A task is done only if:
- the relevant markdown notes are updated,
- symbolic checks run without error,
- and the result is classified as theorem progress or loophole progress.

## Git discipline

Before each major task:
- create a checkpoint commit

After each major task:
- create a checkpoint commit with a one-line scientific summary