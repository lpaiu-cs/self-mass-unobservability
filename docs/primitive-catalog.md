# Primitive Catalog

## Exact M3 Delta<=4 Catalog

- Status: Proven. This document fixes the exact primitive field/operator content used for the `Delta_max = 4` normal-form completeness attack.
- Status: Proven. It is not an extra theorem assumption beyond finite external field content; it is the explicit catalog whose normal-form completeness must be checked.

## Primitive External Field Content

| Symbol | Status | Weight | Tensor type | Role |
| --- | --- | --- | --- | --- |
| `E_ij` | Conjectural | 1 | symmetric trace-free rank-2, parity-even | Minimal external tidal field species kept in the M3 attack. |
| `D_tau E_ij` | Proven | 2 | worldline derivative of `E_ij` | Needed for derivative corrections through `Delta<=4`. |
| `nabla_k E_ij` | Proven | 2 | spatial gradient of `E_ij` | Needed for gradient corrections through `Delta<=4`. |
| `D_tau^2 E_ij` | Proven | 3 | second worldline derivative of `E_ij` | Admitted as a reducible candidate block, not as a normal-form primitive. |
| `a_i = D_tau u_i` | Proven | 1 | worldline acceleration | Admitted only to test lower-order EOM reduction; not retained in normal form. |

## Delta<=4 Candidate Scalar Operators

### Irreducible normal-form targets

| Operator | Status | Weight | Meaning |
| --- | --- | --- | --- |
| `E2` | Conjectural | 2 | `E_ij E^ij` |
| `E3` | Conjectural | 3 | `E_i{}^j E_j{}^k E_k{}^i` |
| `E2^2` | Conjectural | 4 | quartic scalar built from `E2` |
| `dotE2` | Conjectural | 4 | `(D_tau E_ij)(D_tau E^ij)` |
| `gradE2` | Conjectural | 4 | `(nabla_k E_ij)(nabla^k E^ij)` |

### Reducible Delta<=4 candidates

| Operator | Status | Weight | Intended reduction channel |
| --- | --- | --- | --- |
| `E_DtE` | Proven | 3 | total derivative of `E2/2` |
| `Dt2_E2` | Proven | 4 | second total derivative of `E2` |
| `E_Dt2E` | Proven | 4 | total-derivative reduction to `-dotE2` |
| `E4` | Proven | 4 | algebraic reduction to `E2^2/2` by the traceless `3x3` Cayley-Hamilton identity |
| `a2` | Proven | 2 | lower-order worldline EOM reduction |
| `aEa` | Proven | 3 | lower-order worldline EOM reduction |

## What This Catalog Claims

- Status: Proven. The reduction script [`../symbolic/normal_form_reduce.py`](../symbolic/normal_form_reduce.py) handles the reducible operators listed above.
- Status: Conjectural. The remaining theorem burden is `normal-form completeness modulo total derivatives and lower-order equations of motion` for this exact `Delta<=4` catalog.
- Status: Conjectural. No additional admissible scalar operator has yet been proven necessary inside this catalog, but exhaustiveness is not yet a theorem.
