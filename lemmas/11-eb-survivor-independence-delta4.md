# Lemma 11 Draft: `E/B` Survivor Independence At `Delta<=4`

## Question

- Status: Conjectural. Is the raw `E/B` survivor list from the `Delta<=4` magnetic-family audit linearly independent under the currently allowed rules?

## Raw `E/B` Survivor Candidates

- Status: Proven. The raw `E/B` survivor candidate list from [`../symbolic/eb_sector_delta4.py`](../symbolic/eb_sector_delta4.py) is

```math
\{E2,\ B2,\ E3,\ EB2,\ E2^2,\ B2^2,\ dotE2,\ dotB2,\ EBDtB,\ E2B2,\ EB\_sq,\ TrE2B2,\ EBEB,\ gradE2,\ divE2,\ mixedGradE2,\ gradB2,\ divB2,\ mixedGradB2\}.
```

- Status: Proven. This raw list contains `19` candidate survivors.

## Exact Rank Check

- Status: Proven. [`../symbolic/eb_survivor_rank_check.py`](../symbolic/eb_survivor_rank_check.py) expands the raw `E/B` survivor candidates into explicit STF component polynomials and computes the exact coefficient rank.
- Status: Proven. The raw `19`-element list has rank `18`.
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

- Status: Proven. This is now promoted into the explicit reduction rules as the mixed quartic STF identity.

## Corrected `E/B` Basis

- Status: Proven. After eliminating `EBEB`, the corrected `E/B` basis is

```math
\{E2,\ B2,\ E3,\ EB2,\ E2^2,\ B2^2,\ dotE2,\ dotB2,\ EBDtB,\ E2B2,\ EB\_sq,\ TrE2B2,\ gradE2,\ divE2,\ mixedGradE2,\ gradB2,\ divB2,\ mixedGradB2\}.
```

- Status: Proven. The corrected `18`-element list has exact rank `18`.
- Status: Proven. Therefore the corrected `E/B` basis is linearly independent as an operator basis modulo the explicitly stated rules.

## Result

- Status: Proven. The raw `E/B` survivor list is not independent.
- Status: Proven. The first explicit dependence relation sits in the mixed quartic algebraic sector.
- Status: Proven. After explicit reduction by that relation, the corrected `E/B` basis is linearly independent.
