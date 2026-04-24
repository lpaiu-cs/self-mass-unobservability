# Lemma 06 Draft: Normal-Form Completeness At Delta_max = 4

## Statement

- Status: Conjectural. In the minimal free-fall sector with the exact catalog in [`../docs/primitive-catalog.md`](../docs/primitive-catalog.md), every listed `Delta<=4` scalar candidate reduces to the normal-form target set

```math
\{E2,\ E3,\ E2^2,\ dotE2,\ gradE2\}
```

modulo total derivatives, lower-order equations of motion, and the traceless `3x3` Cayley-Hamilton identity.

- Status: Conjectural. This is the concrete M3 attack on `normal-form completeness modulo total derivatives and lower-order equations of motion`.

## Explicit Delta<=4 Reduction Path

- Status: Proven. The reduction script [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) verifies the following catalog-internal reductions:

| Candidate | Status | Reduction channel | Normal-form image |
| --- | --- | --- | --- |
| `E_DtE` | Proven | total derivative | `0` |
| `Dt2_E2` | Proven | total derivative | `0` |
| `E_Dt2E` | Proven | total derivative | `-dotE2` |
| `a2` | Proven | lower-order worldline EOM | `0` |
| `aEa` | Proven | lower-order worldline EOM | `0` |
| `E4` | Proven | traceless `3x3` Cayley-Hamilton identity | `E2^2 / 2` |

- Status: Conjectural. No explicit obstruction has been found among the currently enumerated `Delta<=4` candidates in the exact primitive catalog.

## Proof Skeleton

1. Status: Proven. Start from the exact `Delta<=4` catalog in [`../docs/primitive-catalog.md`](../docs/primitive-catalog.md).
2. Status: Proven. Use total-derivative reduction to remove `E_DtE`, `Dt2_E2`, and to rewrite `E_Dt2E` in terms of `dotE2`.
3. Status: Proven. Use the lower-order worldline EOM to remove acceleration insertions `a2` and `aEa`.
4. Status: Proven. Use the traceless `3x3` Cayley-Hamilton identity to reduce `E4` to `E2^2 / 2`.
5. Status: Conjectural. Conclude that the catalog collapses to the normal-form target list above, provided the primitive catalog is exhaustive at `Delta<=4`.

## Exact Remaining Burden

- Status: Conjectural. The exact remaining burden is `normal-form completeness modulo total derivatives and lower-order equations of motion` for the exact `Delta<=4` primitive catalog.
- Status: Conjectural. Concretely, the missing substep is catalog exhaustiveness: prove that no additional admissible `Delta<=4` scalar contraction exists outside the listed candidates once the minimal-sector symmetries and lower-order EOM are imposed.
- Status: Counterexample candidate. A single additional admissible `Delta<=4` scalar that does not reduce to the current target list would be the smallest explicit obstruction.
