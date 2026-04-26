# Family-Admission Theorem

## Sharp No-Go Target

- Status: Proven. The sharp negative target in this repo is no longer a vague adequacy warning.
- Status: Proven. It is the family-admission no-go for minimal-sector uniqueness across explicitly stated family classes.
- Status: Proven. For the audited family classes below, unsuppressed admission generates a new low-order survivor and therefore obstructs any stronger theorem that tries to identify a unique physically justified minimal free-fall sector without extra suppression assumptions.

## Family Classes Covered

- Status: Proven. STF rank-2 families are kept as a separate special class.
- Status: Proven. The higher-rank STF branch is now split into a universal self-witness threshold theorem for genuine parity-even STF primitive families with `L \ge 3` and a separate mixed-pattern classification.
- Status: Proven. Class R0a: unsuppressed scalar-like rank-0 families with bare source `S`.
- Status: Proven. Class R0b: derivative-only or shift-symmetric scalar families whose primitive catalog excludes bare `S` but admits `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S`.
- Status: Proven. Class R1: genuine local parity-even rank-1 vector families `V_i`, excluding derivative-generated vectors already belonging to audited scalar or rank-2 families.
- Status: Proven. The audited higher-rank STF instances now run through `L = 3, 4, 5, 6`, represented by `T_{ijk}`, `Q_{ijkl}`, `U_{ijklm}`, and `Z_{ijklmn}`.

## Theorem Statement

- Status: Proven. For the rank-2 STF special class, the quadratic invariant `X2 = Tr(X^2)` is a new parity-even low-order survivor unless an explicit extra identity, ordering rule, or background restriction is imposed.
- Status: Proven. The rank-2 STF class remains mixed-aware because the mixed quadratic witness `EX` already exists.
- Status: Proven. For Class R0a, the bare scalar `S` itself is the smallest new survivor.
- Status: Proven. For Class R0b, the smallest new survivors appear at weight `4`; a canonical witness is `dotS2`, and the full smallest-weight audited witness set is `\{DtS_B2,\ dotS2,\ DtS_E2,\ divEGradS,\ gradS2\}`.
- Status: Proven. For Class R1, the first self witness is `V2` at weight `2`, while the first mixed witness is `EVV` at weight `3`.
- Status: Proven. For every genuinely new parity-even STF primitive family with rank `L \ge 3`, the first unavoidable self witness is the quadratic norm `Y2`, no linear scalar witness exists, and no mixed quadratic scalar `E Y_L` exists.
- Status: Proven. Therefore the higher-rank STF branch obeys a universal self-only sharp threshold `w_Y \ge 3` at `\Delta \le 4`.
- Status: Proven. The mixed-pattern story remains separate: audited ranks `3`, `5`, and `6` support the `Y2 / EYY` first-witness pattern, while audited rank `4` is an explicit mixed exception because `EEQ` survives alongside `EQQ`.
- Status: Proven. Therefore minimal-sector uniqueness is not stable under unsuppressed admission of these audited family classes.

## What The Theorem Does Not Say

- Status: Proven. This is not a universal theorem for every imaginable primitive family.
- Status: Proven. The no-go is only claimed for the explicitly stated family classes above.
- Status: Proven. The theorem does not say that broad finite-family fixed-order collapse fails.
- Status: Proven. The theorem does not say that every new family causes fixed-order blowup; it only says that uniqueness fails generically across the audited classes unless extra suppressions are written down.

## Operational Consequence

- Status: Proven. The positive theorem target must therefore be finite-family fixed-order collapse, not minimal-sector uniqueness.
- Status: Proven. New family-admission audits are classification steps: a new witness reinforces the uniqueness no-go unless it also destroys fixed-order finiteness.
- Status: Proven. The audited-set composition question is now closed positively for the enlarged audited family set `{R2, R0a, R0b, R1, Rodd+, Reven4+, Rodd5+, Reven6+}` at `\Delta \le 4`.
- Status: Proven. The former live bottleneck `"prove or refute irreducible family-envelope closure beyond the audited scalar/vector/STF classes"` is now resolved positively.
- Status: Proven. The former live bottleneck `"prove or refute the positive finite-family collapse theorem once the irreducible scalar/vector/STF envelope is fixed"` is now resolved positively.
- Status: Proven. The negative uniqueness no-go now sits alongside a closed positive finite-family fixed-order collapse theorem inside the current theorem domain.
- Status: Proven. [`../symbolic/family_witness_map.py`](../symbolic/family_witness_map.py) records the current witness table linking each audited family class to the exact theorem layer it obstructs.
- Status: Proven. [`suppression-budget-theorem.md`](suppression-budget-theorem.md) and [`../symbolic/witness_threshold_map.py`](../symbolic/witness_threshold_map.py) now promote those witnesses into necessary suppression thresholds for the audited classes.
