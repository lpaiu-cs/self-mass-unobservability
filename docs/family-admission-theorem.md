# Family-Admission Theorem

## Sharp No-Go Target

- Status: Proven. The sharp negative target in this repo is no longer a vague adequacy warning.
- Status: Proven. It is the family-admission no-go for minimal-sector uniqueness across explicitly stated family classes.
- Status: Proven. For the audited family classes below, unsuppressed admission generates a new low-order survivor and therefore obstructs any stronger theorem that tries to identify a unique physically justified minimal free-fall sector without extra suppression assumptions.

## Family Classes Covered

- Status: Proven. Class R2: unsuppressed STF rank-2 families admitted into the parity-even scalar sector at weight `1`.
- Status: Proven. Class R0a: unsuppressed scalar-like rank-0 families with bare source `S`.
- Status: Proven. Class R0b: derivative-only or shift-symmetric scalar families whose primitive catalog excludes bare `S` but admits `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S`.

## Theorem Statement

- Status: Proven. For Class R2, the quadratic invariant `X2 = Tr(X^2)` is a new parity-even low-order survivor unless an explicit extra identity, ordering rule, or background restriction is imposed.
- Status: Proven. For Class R0a, the bare scalar `S` itself is the smallest new survivor.
- Status: Proven. For Class R0b, the smallest new survivors appear at weight `4`; a canonical witness is `dotS2`, and the full smallest-weight audited witness set is `\{DtS_B2,\ dotS2,\ DtS_E2,\ divEGradS,\ gradS2\}`.
- Status: Proven. Therefore minimal-sector uniqueness is not stable under unsuppressed admission of these audited family classes.

## What The Theorem Does Not Say

- Status: Proven. This is not a universal theorem for every imaginable primitive family.
- Status: Proven. The no-go is only claimed for the explicitly stated family classes above.
- Status: Proven. The theorem does not say that broad finite-family fixed-order collapse fails.
- Status: Proven. The theorem does not say that every new family causes fixed-order blowup; it only says that uniqueness fails generically across the audited classes unless extra suppressions are written down.

## Operational Consequence

- Status: Proven. The positive theorem target must therefore be finite-family fixed-order collapse, not minimal-sector uniqueness.
- Status: Proven. New family-admission audits are classification steps: a new witness reinforces the uniqueness no-go unless it also destroys fixed-order finiteness.
- Status: Proven. The audited-set composition question is now closed positively for the already audited family classes at `\Delta \le 4`.
- Status: Proven. The next live bottleneck is family-envelope completeness or the next smallest unaudited family obstruction.
- Status: Proven. [`../symbolic/family_witness_map.py`](../symbolic/family_witness_map.py) records the current witness table linking each audited family class to the exact theorem layer it obstructs.
- Status: Proven. [`suppression-budget-theorem.md`](suppression-budget-theorem.md) and [`../symbolic/witness_threshold_map.py`](../symbolic/witness_threshold_map.py) now promote those witnesses into necessary suppression thresholds for the audited classes.
