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
