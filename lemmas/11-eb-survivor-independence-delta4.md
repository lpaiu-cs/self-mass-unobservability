# Lemma 11: `E/B` Survivor Independence At `Delta<=4` (Corrected)

## Question

- Status: Proven. Is the raw `E/B` survivor list from the `Delta<=4` magnetic-family audit linearly independent under the currently allowed rules?

## Raw `E/B` Survivor Candidates

- Status: Proven. The raw `E/B` survivor candidate list from [`../symbolic/eb_sector_delta4.py`](../symbolic/eb_sector_delta4.py) is

```math
\{E2,\ B2,\ E3,\ EB2,\ E2^2,\ B2^2,\ dotE2,\ dotB2,\ EBDtB,\ E2B2,\ EB\_sq,\ TrE2B2,\ EBEB,\ gradE2,\ gradB2,\ divB2,\ mixedGradB2\}.
```

- Status: Note. This raw list contains `17` candidate survivors.
- Status: Note (corrected 2026-07-12). The pre-correction raw list additionally carried `divE2` and `mixedGradE2`; those were artifacts of the generic gradient model. Under the STF-3 kinematics of `\nabla E` (Schwarz total symmetry + vacuum trace-free), `divE2 = 0` and `mixedGradE2 = gradE2`, so neither is a candidate operator. The magnetic gradient triple `\{gradB2, divB2, mixedGradB2\}` is retained: `B` is an independent primitive with no assumed potential structure. See [`07-gradient-sector-audit.md`](07-gradient-sector-audit.md).

## Exact Rank Check

- Status: Proven. [`../symbolic/eb_survivor_rank_check.py`](../symbolic/eb_survivor_rank_check.py) expands the raw `E/B` survivor candidates into explicit STF component polynomials (with `\nabla E` on its 7-parameter STF-3 parametrization) and computes the exact coefficient rank.
- Status: Proven. The raw `17`-element list has rank `16`.
- Status: Proven. Therefore the raw `E/B` survivor list is not linearly independent.

## First Explicit Dependence Relation

- Status: Proven. The nullspace is one-dimensional and yields the exact mixed quartic relation

```math
EBEB + 2 TrE2B2 - EB\_sq - \frac12 E2B2 = 0.
```

- Status: Proven. Equivalently,

```math
EBEB = EB\_sq + \frac12 E2B2 - 2 TrE2B2.
```

- Status: Note. This is now promoted into the explicit reduction rules as the mixed quartic STF identity.

## Corrected `E/B` Basis

- Status: Proven. After eliminating `EBEB`, the corrected `E/B` basis is

```math
\{E2,\ B2,\ E3,\ EB2,\ E2^2,\ B2^2,\ dotE2,\ dotB2,\ EBDtB,\ E2B2,\ EB\_sq,\ TrE2B2,\ gradE2,\ gradB2,\ divB2,\ mixedGradB2\}.
```

- Status: Proven. The corrected `16`-element list has exact rank `16`.
- Status: Proven. Therefore the corrected `E/B` basis is linearly independent as an operator basis modulo the explicitly stated rules.

## Result

- Status: Proven. The raw `E/B` survivor list is not independent.
- Status: Proven. The first explicit dependence relation sits in the mixed quartic algebraic sector.
- Status: Proven. After explicit reduction by that relation, the corrected `E/B` basis is linearly independent.
