# Assumptions Ledger

This ledger records the active assumptions behind the free-fall sensitivity-collapse target. Each scientific entry carries one of the required status labels.

## Active Theorem Assumptions

| ID | Status | Assumption | Why it is present | If dropped |
| --- | --- | --- | --- | --- |
| A1 | Conjectural | The body is quasi-static on the orbital timescale. | Needed to eliminate fast internal readout variables from the free-fall EFT. | Time-dependent internal modes can appear as explicit state variables. |
| A2 | Conjectural | The body is nearly spherical and nonspinning in the MVP sector. | Keeps the leading body dependence in the monopole sector plus higher multipoles. | Shape or spin data can enter at the same order as the putative sensitivities. |
| A3 | Conjectural | A local worldline EFT exists after integrating out short-distance structure. | Supplies the operator classification language used by Theorem A. | Nonlocal kernels can survive directly in the effective action. |
| A4 | Conjectural | There is no orbital-timescale internal state variable in the free-fall sector. | This is the sharp assumption that rules out `chi`-type hidden coordinates. | A `chi` state can carry body memory that is not reducible to instantaneous sensitivities. |
| A5 | Conjectural | Couplings to the relevant external invariants are analytic near the reference background. | Needed to define a finite Taylor jet of sensitivity coordinates at fixed order. | Threshold or cusp behavior can evade a polynomial sensitivity jet. |
| A6 | Conjectural | The object admits a self-bound equilibrium before external perturbations are applied. | Separates body formation from later passive coupling to external gravity. | Otherwise the theorem can mix equilibrium failure with observational coupling. |
| A7 | Conjectural | The theorem is stated at fixed order in the local derivative and multipole expansion. | This prevents the proof target from silently becoming an all-orders closure claim. | Without fixed-order truncation, finite-dimensionality needs an extra closure theorem. |

## Imported Exclusions From Earlier Work

| ID | Status | Imported statement | Source | Use in the new theorem-first repo |
| --- | --- | --- | --- | --- |
| E1 | Imported from prior work | Literal internal self-unobservability does not produce a viable compact self-bound equilibrium. | [`request2/REQUEST2_INTERNAL_STRUCTURE.md`](../request2/REQUEST2_INTERNAL_STRUCTURE.md) | Removes the old "change the stellar structure equations directly" branch from Theorem A. |
| E2 | Imported from prior work | COM self-subtraction does not create a new monopole force in the free-fall sector. | [`request1/REQUEST1_COM_DECOUPLING.md`](../request1/REQUEST1_COM_DECOUPLING.md) | Prevents the proof from smuggling in a self-force through coordinate bookkeeping. |
| E3 | Imported from prior work | The strongest previous global result is a provisional consistency scaffold, not a theorem. | [`request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md`](../request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md) | Justifies pivoting the new repo toward proof-or-loophole classification. |

## Current Minimal Missing Assumption

- Status: Conjectural. The present draft still needs a clean statement that the chosen invariant basis is finite at the truncation order of interest; without that bookkeeping assumption, "finite-dimensional sensitivity manifold" risks becoming ambiguous.
