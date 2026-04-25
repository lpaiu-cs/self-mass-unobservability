# Family-Admission Theorem

## Sharp No-Go Target

- Status: Proven. The sharp negative target in this repo is no longer a vague adequacy warning.
- Status: Proven. It is the family-admission no-go for minimal-sector uniqueness across explicitly stated family classes.
- Status: Proven. For the audited family classes below, unsuppressed admission generates a new low-order survivor and therefore obstructs any stronger theorem that tries to identify a unique physically justified minimal free-fall sector without extra suppression assumptions.

## Family Classes Covered

- Status: Proven. Class R2: unsuppressed STF rank-2 families admitted into the parity-even scalar sector at weight `1`.
- Status: Proven. Class R0a: unsuppressed scalar-like rank-0 families with bare source `S`.
- Status: Proven. Class R0b: derivative-only or shift-symmetric scalar families whose primitive catalog excludes bare `S` but admits `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S`.
- Status: Proven. Class R1: genuine local parity-even rank-1 vector families `V_i`, excluding derivative-generated vectors already belonging to audited scalar or rank-2 families.
- Status: Proven. Class Rodd+: genuine local parity-even fully symmetric trace-free rank-3 tensor families `T_{ijk}`, excluding derivative-generated rank-3 blocks already belonging to audited scalar, STF rank-2, or vector families.
- Status: Proven. Class Reven4+: genuine local parity-even fully symmetric trace-free rank-4 tensor families `Q_{ijkl}`, excluding trace descendants reducible to audited lower-rank families and derivative-generated rank-4 blocks already belonging to audited scalar, STF rank-2, vector, or rank-3 families.
- Status: Proven. Class Rodd5+: genuine local parity-even fully symmetric trace-free rank-5 tensor families `U_{ijklm}`, excluding trace descendants reducible to audited lower-rank families and derivative-generated rank-5 blocks already belonging to audited scalar, STF rank-2, vector, rank-3, or rank-4 families.

## Theorem Statement

- Status: Proven. For Class R2, the quadratic invariant `X2 = Tr(X^2)` is a new parity-even low-order survivor unless an explicit extra identity, ordering rule, or background restriction is imposed.
- Status: Proven. For Class R0a, the bare scalar `S` itself is the smallest new survivor.
- Status: Proven. For Class R0b, the smallest new survivors appear at weight `4`; a canonical witness is `dotS2`, and the full smallest-weight audited witness set is `\{DtS_B2,\ dotS2,\ DtS_E2,\ divEGradS,\ gradS2\}`.
- Status: Proven. For Class R1, the first self witness is `V2` at weight `2`, while the first mixed witness is `EVV` at weight `3`.
- Status: Proven. For Class Rodd+, the first self witness is `T2` at weight `2`, while the first mixed witness is `ETT` at weight `3`.
- Status: Proven. For Class Reven4+, the first self witness is `Q2` at weight `2`, while the first mixed witness is `EQQ` at weight `3`.
- Status: Proven. For Class Rodd5+, the first self witness is `U2` at weight `2`, while the first mixed witness is `EUU` at weight `3`.
- Status: Proven. Therefore minimal-sector uniqueness is not stable under unsuppressed admission of these audited family classes.

## What The Theorem Does Not Say

- Status: Proven. This is not a universal theorem for every imaginable primitive family.
- Status: Proven. The no-go is only claimed for the explicitly stated family classes above.
- Status: Proven. The theorem does not say that broad finite-family fixed-order collapse fails.
- Status: Proven. The theorem does not say that every new family causes fixed-order blowup; it only says that uniqueness fails generically across the audited classes unless extra suppressions are written down.

## Operational Consequence

- Status: Proven. The positive theorem target must therefore be finite-family fixed-order collapse, not minimal-sector uniqueness.
- Status: Proven. New family-admission audits are classification steps: a new witness reinforces the uniqueness no-go unless it also destroys fixed-order finiteness.
- Status: Proven. The audited-set composition question is now closed positively for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+}` at `\Delta \le 4`.
- Status: Proven. The current live bottleneck is no longer audited-set composition, but family-envelope completeness and the next smallest unaudited family class `Reven6+`.
- Status: Proven. [`../symbolic/family_witness_map.py`](../symbolic/family_witness_map.py) records the current witness table linking each audited family class to the exact theorem layer it obstructs.
- Status: Proven. [`suppression-budget-theorem.md`](suppression-budget-theorem.md) and [`../symbolic/witness_threshold_map.py`](../symbolic/witness_threshold_map.py) now promote those witnesses into necessary suppression thresholds for the audited classes.
