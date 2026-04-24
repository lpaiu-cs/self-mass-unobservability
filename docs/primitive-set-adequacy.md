# Primitive-Set Adequacy

## Purpose

- Status: Proven. This note separates two statements that must not be conflated:
- Status: Proven. Statement 1: theorem candidate for the exact current primitive set.
- Status: Proven. Statement 2: theorem candidate for a physically justified minimal free-fall sector.

## What Makes The Current Catalog Minimal

- Status: Proven. The current catalog is minimal only in the weak operational sense that it keeps the smallest electric-sector family already needed for the free-fall `Delta<=4` contraction audit:
  `E_ij`, `D_tau E_ij`, `nabla_k E_ij`, `D_tau^2 E_ij`, and the test-EOM block `a_i`.
- Status: Proven. It is minimal relative to the current reduction problem because every retained block either survives or is needed to test an allowed reduction channel explicitly.
- Status: Proven. It is not yet minimal in the stronger physical sense of being the complete parity-even nonspinning free-fall sector.

## What Is Still Only A Working Restriction

- Status: Conjectural. The current catalog excludes any additional primitive family not generated from the electric tidal STF `E_ij`.
- Status: Conjectural. In particular, it excludes the magnetic tidal STF family `B_ij`.
- Status: Conjectural. That exclusion is currently a working restriction of the exact-current-set theorem, not a justified consequence of the stated symmetry assumptions alone.

## Smallest Primitive-Family Attack

- Status: Proven. The smallest physically motivated adequacy attack is to add the magnetic tidal STF family `B_ij` while leaving the rest of the exact current primitive set unchanged.
- Status: Conjectural. The magnetic tidal STF family `B_ij` is the most conservative next family to test because it is parity-compatible with a nonspinning scalar action once it enters quadratically.
- Status: Proven. [`../symbolic/primitive_family_attack.py`](../symbolic/primitive_family_attack.py) performs this one-family extension at `Delta<=4`.
- Status: Proven. [`magnetic-family-ordering.md`](magnetic-family-ordering.md) now records the explicit ordering question rather than leaving `B_ij` tacitly excluded.
- Status: Proven. [`../symbolic/eb_sector_delta4.py`](../symbolic/eb_sector_delta4.py) computes the full `E/B`-expanded `Delta<=4` survivor list under the currently allowed rules.

## Adequacy Result

- Status: Proven. The magnetic-family extension produces new surviving scalars.
- Status: Proven. The smallest explicit new survivor is `B2`.
- Status: Proven. In the smallest `B`-only spot attack, additional new surviving candidates already appear at the same fixed order, including `EB2`, `B2DtE`, `E2B2`, `EB_sq`, `TrE2B2`, `EBEB`, and `B2^2`.
- Status: Proven. In the fuller `E/B`-expanded audit, the mixed time-derivative content is more naturally represented by the surviving class `EBDtB`, while `B2DtE` becomes reducible by an explicit total-derivative relation.
- Status: Proven. Therefore the exact current primitive set is not yet adequate as a physically justified minimal sector unless one adds a new explicit assumption excluding the magnetic family.
- Status: Proven. The `E/B` enlargement itself still closes on a finite explicit `Delta<=4` survivor list, so adequacy is the live issue rather than fixed-order blowup.

## What This Gains

- Status: Proven. The M5 bottleneck is now sharply localized.
- Status: Proven. Survivor independence is no longer the active risk for the exact current primitive set.
- Status: Proven. Primitive-set adequacy, sharpened to magnetic-family ordering, is the active risk for any stronger theorem that aims to describe a physically justified minimal free-fall sector.
