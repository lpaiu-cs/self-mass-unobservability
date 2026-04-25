# Assumptions Ledger

This ledger records the active assumptions behind the free-fall sensitivity-collapse target. Each scientific entry carries one of the required status labels.

## Active Theorem Assumptions

| ID | Status | Assumption | Why it is present | If dropped |
| --- | --- | --- | --- | --- |
| A1 | Conjectural | The body is quasi-static on the orbital timescale. | Needed to eliminate fast internal readout variables from the free-fall EFT. | Time-dependent internal modes can appear as explicit state variables. |
| A2 | Conjectural | The minimal M4 sector is nearly spherical, nonspinning, and parity-even. | Keeps the operator search space inside the minimal free-fall sector requested for this theorem candidate. | Spin, parity-odd structures, or large asphericity add new operator families at the same order. |
| A3 | Conjectural | A local worldline EFT exists after integrating out short-distance structure. | Supplies the operator classification language used by Theorem A. | Nonlocal kernels can survive directly in the effective action. |
| A4 | Conjectural | There is no orbital-timescale internal state variable in the free-fall sector. | This is the sharp assumption that rules out `chi`-type hidden coordinates. | A `chi` state can carry body memory that is not reducible to instantaneous sensitivities. |
| A5 | Conjectural | Couplings to the relevant external invariants are analytic near the reference background. | Needed to define a finite Taylor jet of sensitivity coordinates at fixed order. | Threshold or cusp behavior can evade a polynomial sensitivity jet. |
| A6 | Conjectural | The object admits a self-bound equilibrium before external perturbations are applied. | Separates body formation from later passive coupling to external gravity. | Otherwise the theorem can mix equilibrium failure with observational coupling. |
| A7 | Conjectural | The theorem is stated at fixed order in the operator counting rule of [`power-counting.md`](power-counting.md). | This prevents the proof target from silently becoming an all-orders closure claim. | Without fixed-order truncation, finite-dimensionality needs an extra closure theorem. |
| A8 | Conjectural | The minimal sector uses finite external field content: only finitely many primitive external tensor species are admitted at the chosen order. | This is the explicit hypothesis that replaces the older hidden finite-basis assumption. | An uncontrolled field catalog can generate infinitely many primitive directions before normal-form reduction starts. |

## Imported Exclusions From Earlier Work

| ID | Status | Imported statement | Source | Use in the new theorem-first repo |
| --- | --- | --- | --- | --- |
| E1 | Imported from prior work | Literal internal self-unobservability does not produce a viable compact self-bound equilibrium. | [`request2/REQUEST2_INTERNAL_STRUCTURE.md`](../request2/REQUEST2_INTERNAL_STRUCTURE.md) | Removes the old "change the stellar structure equations directly" branch from Theorem A. |
| E2 | Imported from prior work | COM self-subtraction does not create a new monopole force in the free-fall sector. | [`request1/REQUEST1_COM_DECOUPLING.md`](../request1/REQUEST1_COM_DECOUPLING.md) | Prevents the proof from smuggling in a self-force through coordinate bookkeeping. |
| E3 | Imported from prior work | The strongest previous global result is a provisional consistency scaffold, not a theorem. | [`request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md`](../request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md) | Justifies pivoting the new repo toward proof-or-loophole classification. |

## Non-Assumption For M4

- Status: Proven. Finite invariant basis closure is not listed as an assumption in M4.
- Status: Proven. No transversality condition such as `\nabla_i E^{ij} = 0` is assumed in the current `Delta<=4` audit.
- Status: Proven. No vacuum/Bianchi identity is assumed to collapse `divE2` or `mixedGradE2` into `gradE2`.
- Status: Proven. Contraction-level exhaustiveness for the exact current primitive blocks is attacked by explicit enumeration rather than by absence arguments.

## Non-Assumption For M5

- Status: Proven. Linear survivor independence is now attacked directly by an exact polynomial rank check rather than inferred from non-reducibility.
- Status: Proven. The exact-current-set theorem does not claim that the present primitive catalog already exhausts every physically admissible minimal-sector family.
- Status: Proven. Excluding the magnetic tidal family `B_ij` is currently a working restriction of the exact-current-set theorem, not a derived completeness statement.

## Non-Assumption For M6

- Status: Proven. No magnetic-ordering assumption is currently active that suppresses `B_ij` relative to `E_ij`.
- Status: Proven. The repository does not currently assume `B_ij = 0` on the allowed backgrounds.
- Status: Proven. Therefore any stronger electric-only minimal-sector theorem would require a new explicit assumption, recorded rather than implied.

## Explicit Rule For M7

- Status: Proven. The mixed quartic `E/B` STF identity used to eliminate `EBEB` is now explicit in [`reduction-rules.md`](reduction-rules.md); it is not being used tacitly.
- Status: Proven. No stronger magnetic-ordering salvage assumption is currently active.

## Non-Assumption For M8

- Status: Proven. No scalar-family ordering assumption is currently active that suppresses a parity-even rank-0 external family `S`.
- Status: Proven. No derivative-only or shift-symmetric scalar rule is being used tacitly.
- Status: Proven. No scalar-background restriction such as `S = 0`, `S = const.`, or `\nabla_i S = 0` is currently assumed.

## Explicit Rule For M9

- Status: Proven. The derivative-only scalar audit is a separate family-admission test, not a new global theorem assumption.
- Status: Proven. In that audit, non-shift-symmetric total-derivative preimages such as `D_\tau(S X)` are not used to reduce operators, because bare `S` is excluded by construction.
- Status: Proven. No cross-family orthogonality or alignment condition such as `Tr(EX) = 0`, `SE2 = 0`, or `DtS_E2 = 0` is assumed in the mixed-witness threshold audit.
- Status: Proven. The mixed-witness threshold notes do not infer mixed-sector harmlessness from self-witness analysis alone; self and mixed survivors are checked separately.

## Non-Assumption For M12

- Status: Proven. The repository does not currently assume that the audited family set `{R2, R0a, R0b, R1}` already exhausts the parity-even local MVP primitive-family envelope.
- Status: Proven. No explicit assumption restricts admitted parity-even local primitive families to scalar and rank-2 classes only.
- Status: Proven. Therefore local parity-even higher-rank tensor families remain open family-envelope targets unless they are audited directly or excluded by a new explicit assumption.

## Non-Assumption For M13

- Status: Proven. No current MVP assumption says that every local parity-even vector block must be derivative-generated from an audited scalar or rank-2 family.
- Status: Proven. Derivative-generated vectors such as `\nabla_i S` or `\nabla_j X_i{}^j` are assigned to their parent families and are not counted as new primitive-family admissions.
- Status: Proven. No vector-family ordering, absorption rule, or vector-background restriction is currently active that would remove a genuine primitive family `V_i` from the theorem domain.

## Explicit Representative Choice For M14

- Status: Proven. No current MVP assumption says that every local parity-even rank-3 tensor block must be derivative-generated from an audited lower-rank family.
- Status: Proven. Derivative-generated rank-3 blocks such as `\nabla E`, `\nabla B`, `\nabla V`, or scalar-derivative descendants are assigned to their parent families and are not counted as new primitive-family admissions.
- Status: Proven. The `Rodd+` audit uses a fully symmetric trace-free rank-3 representative explicitly, because trace descendants reduce to the already audited vector gate `R1`, while parity-odd or spin-like antisymmetric pieces are already outside the current MVP theorem domain by `A2`.
- Status: Proven. No rank-3 ordering, absorption rule, or rank-3 background restriction is currently active that would remove a genuine primitive family `T_{ijk}` from the theorem domain.

## Explicit Representative Choice For M15

- Status: Proven. No current MVP assumption says that every local parity-even rank-5 tensor block must be derivative-generated from an audited lower-rank family.
- Status: Proven. Derivative-generated rank-5 blocks such as `\nabla Q`, `\nabla T`, higher-gradient vector descendants, or scalar descendants dressed only by derivatives are assigned to their parent families and are not counted as new primitive-family admissions.
- Status: Proven. The `Rodd5+` audit uses a fully symmetric trace-free rank-5 representative explicitly, because trace descendants reduce to already audited lower-rank classes, while antisymmetric or spin-carrying pieces are already outside the current MVP theorem domain by `A2`.
- Status: Proven. No rank-5 ordering, absorption rule, or rank-5 background restriction is currently active that would remove a genuine primitive family `U_{ijklm}` from the theorem domain.
