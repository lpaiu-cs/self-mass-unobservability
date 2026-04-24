# Lemma 06 Draft: Pre-M4 Normal-Form Completeness Attempt At Delta_max = 4

## Statement

- Status: Proven. This note records the pre-M4 internal reduction attempt before the full contraction-level audit was added.
- Status: Conjectural. In the minimal free-fall sector with the pre-M4 catalog, every listed `Delta<=4` scalar candidate reduces to the normal-form target set

```math
\{E2,\ E3,\ E2^2,\ dotE2,\ gradE2\}
```

modulo total derivatives, lower-order equations of motion, and the traceless `3x3` Cayley-Hamilton identity.

- Status: Conjectural. This is the concrete M3 attack on normal-form reduction before catalog exhaustiveness was checked explicitly.

## Explicit Delta<=4 Reduction Path

- Status: Proven. The reduction script [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) verifies the catalog-internal reductions used at that stage.
- Status: Proven. M4 later showed that the pre-M4 catalog itself omitted `divE2` and `mixedGradE2`, so this lemma is not the final contraction-level statement.

## Proof Skeleton

1. Status: Proven. Start from the pre-M4 `Delta<=4` catalog.
2. Status: Proven. Use total-derivative reduction to remove `E_DtE`, `Dt2_E2`, and to rewrite `E_Dt2E` in terms of `dotE2`.
3. Status: Proven. Use the lower-order worldline EOM to remove acceleration insertions.
4. Status: Proven. Use the traceless `3x3` Cayley-Hamilton identity to reduce `E4` to `E2^2 / 2`.
5. Status: Conjectural. Conclude that the listed pre-M4 catalog collapses to the old five-element target, provided the catalog itself was exhaustive.

## Exact Historical Gap

- Status: Proven. The exact missing substep identified here was catalog exhaustiveness.
- Status: Proven. M4 resolved that omission by explicit contraction enumeration and found the omitted surviving gradient operators `divE2` and `mixedGradE2`.
- Status: Conjectural. The corrected contraction-level statement now lives in [`07-gradient-sector-audit.md`](07-gradient-sector-audit.md), [`08-mixed-time-derivative-audit.md`](08-mixed-time-derivative-audit.md), and [`../docs/theorem-A-freefall.md`](../docs/theorem-A-freefall.md).
