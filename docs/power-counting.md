# Power Counting

## M2 Fixed-Order Rule

- Status: Proven. The current M2 symbolic enumerator uses a fixed truncation order `\Delta_{\max} = 4`.
- Status: Conjectural. The same counting rule is intended to scale to higher `\Delta_{\max}`, but only the `\Delta_{\max} = 4` catalog is exercised in code at present.
- Status: Proven. Let the external field content be a finite set of primitive parity-even tensor species `{\cal F} = \{X_a\}_{a=1}^{r}`.
- Status: Proven. Assign each primitive species a positive intrinsic weight `w(X_a) \ge 1`.
- Status: Proven. Assign `w(D_\tau) = 1` and `w(\nabla) = 1`.
- Status: Proven. A decorated building block `D_\tau^{p}\nabla^{q}X_a` has weight

```math
w\!\left(D_\tau^{p}\nabla^{q}X_a\right)
=
w(X_a) + p + q.
```

- Status: Proven. A composite operator has weight equal to the sum of the weights of its factors, and only operators with total weight `\Delta \le \Delta_{\max}` are retained.
- Status: Proven. In the minimal sector, only nonspinning, nearly spherical, parity-even worldline scalars contribute to the monopole response, while finite-rank parity-even tensors contribute to higher multipoles.

## Candidate Scalar Catalog Used In Code

- Status: Conjectural. The current code path uses a candidate scalar generator list built from one parity-even electric-type tidal tensor family `E_{ij}`:

| Generator | Status | Weight | Intended meaning |
| --- | --- | --- | --- |
| `E2` | Conjectural | 2 | `E_{ij}E^{ij}` |
| `E3` | Conjectural | 3 | `E_i{}^j E_j{}^k E_k{}^i` |
| `dotE2` | Conjectural | 4 | `(D_\tau E_{ij})(D_\tau E^{ij})` |
| `gradE2` | Conjectural | 4 | `(\nabla_k E_{ij})(\nabla^k E^{ij})` |

- Status: Conjectural. This table is a candidate primitive catalog for the M2 basis-closure theorem candidate, not an assumed final basis.
- Status: Proven. Given a fixed candidate catalog with positive weights, the set of monomials of weight `\le \Delta_{\max}` is finite.

## What Is Proven Versus Unresolved

- Status: Proven. Fixed order plus finite external field content implies a finite set of decorated primitive building blocks.
- Status: Proven. Positive weights imply a finite set of candidate monomials built from those primitives.
- Status: Conjectural. The unresolved step is physical completeness: whether every admissible minimal-sector local scalar operator reduces, modulo total derivatives and lower-order equations of motion, to the chosen candidate monomial list.
