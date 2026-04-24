# Lemma 14: Derivative-Only Scalar Audit At `\Delta \le 4`

## Setup

- Status: Proven. Start from the corrected `E/B` exact-current-set audit.
- Status: Proven. Adjoin only the derivative blocks `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S`.
- Status: Proven. Exclude bare `S` from the primitive catalog to model a shift-symmetric or derivative-only scalar family.
- Status: Proven. In this audit, no non-shift-symmetric total-derivative preimage such as `D_\tau(S X)` is used to reduce operators, because that would reintroduce the forbidden bare scalar into the reduction step.

## Symbolic Audit Result

- Status: Proven. [`../symbolic/shift_scalar_sector_delta4.py`](../symbolic/shift_scalar_sector_delta4.py) finds `52` parity-even scalar classes in the `E/B+derivative-only-scalar` sector at `\Delta \le 4`.
- Status: Proven. The surviving list has `23` elements:

```math
\{B2,\ E2,\ EB2,\ E3,\ B2^2,\ DtS\_B2,\ E2B2,\ EB\_sq,\ TrE2B2,\ EBDtB,\ dotB2,\ dotE2,\ dotS2,\ DtS\_E2,\ E2^2,\ divB2,\ gradB2,\ mixedGradB2,\ divE2,\ gradE2,\ mixedGradE2,\ divEGradS,\ gradS2\}.
```

- Status: Proven. The first new survivors beyond the corrected `E/B` basis appear only at weight `4`.
- Status: Proven. The full smallest-weight new-survivor set is
  `\{DtS_B2,\ dotS2,\ DtS_E2,\ divEGradS,\ gradS2\}`.
- Status: Proven. A canonical explicit obstruction is `dotS2`.

## Consequence

- Status: Proven. Shift symmetry removes the bare-scalar obstruction `S`.
- Status: Proven. Shift symmetry does not rescue minimal-sector uniqueness by itself.
- Status: Proven. The derivative-only scalar family still enlarges the corrected `E/B` survivor list at the same fixed truncation order.
- Status: Proven. Therefore a stronger minimal-sector theorem still needs an explicit family-suppression assumption even after the derivative-only scalar test.

## Boundary

- Status: Proven. This audit is evidence against minimal-sector uniqueness, not against the broader finite-family collapse program.
- Status: Conjectural. The remaining positive question is whether the finite-family collapse theorem can be stated for arbitrary explicit admitted catalogs rather than just the audited examples.
