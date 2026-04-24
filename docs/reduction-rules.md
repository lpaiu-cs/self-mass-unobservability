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
| `D_\tau(S) = DtS` | Proven | `DtS` is total-derivative reducible once the scalar-like family is admitted. |
| `D_\tau(S2) = 2 S_DtS` | Proven | `S_DtS` is total-derivative reducible once the scalar-like family is admitted. |
| `D_\tau(DtS) = Dt2S` | Proven | `Dt2S` is total-derivative reducible once the scalar-like family is admitted. |
| `D_\tau(S DtS) = dotS2 + S_Dt2S` | Proven | `S_Dt2S` reduces to `-dotS2` if `dotS2` is chosen as the scalar-family time-derivative normal-form representative. |
| `D_\tau(S E2) = DtS_E2 + 2 S E_DtE` | Proven | `SE_DtE` is reducible to `-DtS_E2 / 2` once the scalar-like family is admitted. |
| `D_\tau(S B2) = DtS_B2 + 2 S B_DtB` | Proven | `S_BDtB` is reducible to `-DtS_B2 / 2` once the scalar-like family is admitted. |

## Allowed Lower-Order EOM Reductions

- Status: Proven. The lower-order free-fall worldline equation is geodesic/test-body motion, so explicit acceleration insertions are removed.

| Rule | Status | Consequence |
| --- | --- | --- |
| `a_i = 0` modulo lower-order worldline EOM | Proven | Every operator with at least one explicit `a_i` factor is reducible. |

- Status: Proven. This covers `a2`, `aEa`, `aDivE`, `aDtEa`, `a2E2`, `aE2a`, `a4`, `aEGradE_1`, `aEGradE_2`, and `aEGradE_3`.
- Status: Proven. Once the scalar-like family is admitted, this also covers `aGradS`, `a2S`, `a2DtS`, `aSEa`, `aEGradS`, `SaDivE`, `SaGradS`, and `a2S2`.

## Allowed Algebraic Identities

| Rule | Status | Consequence |
| --- | --- | --- |
| `Tr(E^4) = \frac12 (Tr(E^2))^2` for a symmetric trace-free `3x3` tensor `E` | Proven | `E4` reduces to `E2^2 / 2`. |
| `Tr(B^4) = \frac12 (Tr(B^2))^2` for a symmetric trace-free `3x3` tensor `B` | Proven | `B4` reduces to `B2^2 / 2` once the magnetic family is admitted. |
| `Tr(EBEB) = (E:B)^2 + \frac12 Tr(E^2) Tr(B^2) - 2 Tr(E^2 B^2)` for symmetric trace-free `3x3` tensors `E, B` | Proven | `EBEB` reduces to `EB_sq + E2B2/2 - 2 TrE2B2` once the mixed quartic STF identity is admitted explicitly. |

## Not Allowed Unless Added Explicitly

- Status: Proven. No transversality condition such as `\nabla_i E^{ij} = 0` is used in the current audit.
- Status: Proven. No vacuum/Bianchi relation is used to collapse `divE2` or `mixedGradE2` into `gradE2`.
- Status: Proven. No parity-odd epsilon-tensor identity is used in the current audit.
- Status: Proven. No extra mixed `E/B` trace identity beyond the explicit quartic STF relation above is used to collapse `E2B2`, `EB_sq`, or `TrE2B2`.
- Status: Proven. No scalar shift symmetry, derivative-only scalar rule, or scalar-background restriction is used unless stated explicitly in [`scalar-family-ordering.md`](scalar-family-ordering.md).
- Status: Proven. In the derivative-only scalar audit, no non-shift-symmetric preimage such as `D_\tau(S X)` is used as a reduction rule because bare `S` is excluded from that test catalog.

## Operational Summary

- Status: Proven. Under the currently allowed rules, `TrE2DtE` is reducible.
- Status: Proven. Under the currently allowed rules, `B_DtB`, `B_Dt2B`, and `B2DtE` are reducible once the magnetic family is admitted.
- Status: Proven. Under the currently allowed rules, `DtS`, `S_DtS`, `Dt2S`, `S_Dt2S`, `SE_DtE`, and `S_BDtB` are reducible once the scalar-like family is admitted.
- Status: Proven. Under the currently allowed rules, `EBEB` is reducible once the explicit mixed quartic STF identity is admitted.
- Status: Proven. Under the currently allowed rules, every acceleration-bearing scalar is reducible.
- Status: Conjectural. Under the currently allowed rules, `divE2` and `mixedGradE2` remain surviving gradient candidates.
- Status: Conjectural. Under the currently allowed rules, `divB2`, `mixedGradB2`, and `EBDtB` remain surviving `E/B`-sector candidates once the magnetic family is admitted.
- Status: Conjectural. Under the currently allowed rules, `S`, `SE2`, `SB2`, `DtS_E2`, `DtS_B2`, `gradS2`, and `divEGradS` remain surviving scalar-family candidates once the scalar-like family is admitted.
- Status: Conjectural. Under the derivative-only scalar rule, `dotS2`, `gradS2`, `DtS_E2`, `DtS_B2`, and `divEGradS` remain surviving candidates even after bare `S` has been excluded.
