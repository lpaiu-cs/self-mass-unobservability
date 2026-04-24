# Reduction Rules

This note lists every reduction rule currently allowed in the `Delta<=4` minimal free-fall audit.

## Allowed Total-Derivative Reductions

- Status: Proven. Only total derivatives along the worldline parameter `\tau` are quotiented out.
- Status: Proven. Spatial total derivatives are not discarded in the worldline action.

| Rule | Status | Consequence |
| --- | --- | --- |
| `D_\tau(E2) = 2 E_DtE` | Proven | `E_DtE` is total-derivative reducible. |
| `D_\tau(E_DtE) = dotE2 + E_Dt2E` | Proven | `E_Dt2E` reduces to `-dotE2`. |
| `D_\tau^2(E2)` is a total derivative | Proven | `Dt2_E2` is discarded. |
| `D_\tau(E3) = 3 TrE2DtE` | Proven | `TrE2DtE` is total-derivative reducible. |
| `D_\tau(B2) = 2 B_DtB` | Proven | `B_DtB` is total-derivative reducible once the magnetic family is admitted. |
| `D_\tau(B_DtB) = dotB2 + B_Dt2B` | Proven | `B_Dt2B` reduces to `-dotB2` once the magnetic family is admitted. |
| `D_\tau(EB2) = B2DtE + 2 EBDtB` | Proven | `B2DtE` is reducible to `-2 EBDtB` if `EBDtB` is chosen as the mixed `E/B` normal-form representative. |

## Allowed Lower-Order EOM Reductions

- Status: Proven. The lower-order free-fall worldline equation is geodesic/test-body motion, so explicit acceleration insertions are removed.

| Rule | Status | Consequence |
| --- | --- | --- |
| `a_i = 0` modulo lower-order worldline EOM | Proven | Every operator with at least one explicit `a_i` factor is reducible. |

- Status: Proven. This covers `a2`, `aEa`, `aDivE`, `aDtEa`, `a2E2`, `aE2a`, `a4`, `aEGradE_1`, `aEGradE_2`, and `aEGradE_3`.

## Allowed Algebraic Identities

| Rule | Status | Consequence |
| --- | --- | --- |
| `Tr(E^4) = \frac12 (Tr(E^2))^2` for a symmetric trace-free `3x3` tensor `E` | Proven | `E4` reduces to `E2^2 / 2`. |
| `Tr(B^4) = \frac12 (Tr(B^2))^2` for a symmetric trace-free `3x3` tensor `B` | Proven | `B4` reduces to `B2^2 / 2` once the magnetic family is admitted. |

## Not Allowed Unless Added Explicitly

- Status: Proven. No transversality condition such as `\nabla_i E^{ij} = 0` is used in the current audit.
- Status: Proven. No vacuum/Bianchi relation is used to collapse `divE2` or `mixedGradE2` into `gradE2`.
- Status: Proven. No parity-odd epsilon-tensor identity is used in the current audit.
- Status: Proven. No extra mixed `E/B` trace identity is used to collapse `E2B2`, `EB_sq`, `TrE2B2`, or `EBEB`.

## Operational Summary

- Status: Proven. Under the currently allowed rules, `TrE2DtE` is reducible.
- Status: Proven. Under the currently allowed rules, `B_DtB`, `B_Dt2B`, and `B2DtE` are reducible once the magnetic family is admitted.
- Status: Proven. Under the currently allowed rules, every acceleration-bearing scalar is reducible.
- Status: Conjectural. Under the currently allowed rules, `divE2` and `mixedGradE2` remain surviving gradient candidates.
- Status: Conjectural. Under the currently allowed rules, `divB2`, `mixedGradB2`, and `EBDtB` remain surviving `E/B`-sector candidates once the magnetic family is admitted.
