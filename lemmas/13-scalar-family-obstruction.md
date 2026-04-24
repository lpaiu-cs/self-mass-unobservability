# Lemma 13: Scalar-Family Obstruction At `\Delta \le 4`

## Setup

- Status: Proven. Start from the corrected `E/B` exact-current-set basis already isolated in [`11-eb-survivor-independence-delta4.md`](11-eb-survivor-independence-delta4.md).
- Status: Proven. Adjoin one parity-even scalar-like external family with primitive blocks
  `S`, `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S`.
- Status: Proven. Use only the currently allowed reduction rules: worldline total derivatives, lower-order free-fall EOM, single-family STF identities, and the explicit mixed quartic `E/B` identity already recorded in [`../docs/reduction-rules.md`](../docs/reduction-rules.md).

## Symbolic Audit Result

- Status: Proven. [`../symbolic/es_sector_delta4.py`](../symbolic/es_sector_delta4.py) finds `72` parity-even scalar classes in the `E/B+scalar` sector at `\Delta \le 4`.
- Status: Proven. The surviving list has `33` elements:

```math
\{S,\ B2,\ E2,\ S2,\ EB2,\ SB2,\ E3,\ SE2,\ S3,\ B2^2,\ DtS\_B2,\ E2B2,\ EB\_sq,\ TrE2B2,\ SEB2,\ S2B2,\ EBDtB,\ dotB2,\ dotE2,\ dotS2,\ DtS\_E2,\ E2^2,\ SE3,\ S2E2,\ divB2,\ gradB2,\ mixedGradB2,\ divE2,\ gradE2,\ mixedGradE2,\ divEGradS,\ gradS2,\ S4\}.
```

- Status: Proven. The first explicit new survivor beyond the corrected `E/B` basis is `S`.
- Status: Proven. [`../symbolic/es_survivor_rank_check.py`](../symbolic/es_survivor_rank_check.py) finds rank `33` out of `33`, so no new linear dependence appears at this stage.

## Obstruction Statement

- Status: Proven. The scalar-family extension falsifies any attempt to treat the corrected `E/B` exact-current-set theorem as though it already described a physically justified minimal free-fall sector without extra assumptions.
- Status: Proven. The smallest explicit obstruction is the bare scalar survivor `S`.
- Status: Proven. This obstruction is sharper than a vague completeness worry because it identifies a specific additional primitive family and the first operator it contributes.
- Status: Proven. The obstruction does not presently falsify the broader fixed-order collapse program, because the enlarged `E/B+scalar` sector still closes on a finite linearly independent basis at `\Delta \le 4`.

## Boundary Of What Is Still Open

- Status: Conjectural. A stronger minimal-sector theorem could still be salvaged if one adds an explicit scalar-ordering assumption, derivative-only rule, or background restriction.
- Status: Conjectural. Without such an added assumption, the corrected `E/B` exact-current-set theorem should not be advertised as a physically complete minimal-sector statement.
