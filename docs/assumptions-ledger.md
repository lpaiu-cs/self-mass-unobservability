# Assumptions Ledger

This ledger records the active assumptions behind the free-fall sensitivity-collapse target.

**Status label for `A1`-`A9`: `Assumption`.** These are posited hypotheses that
define the theorem domain -- they are neither `Proven` (they are not derived) nor
`Conjectural` (they are not open questions); they are the premises the positive
theorem is conditional on. The `If dropped` column records the counterexample or
failure mode that follows when a premise is removed. (Earlier revisions tagged them
`Conjectural`, which conflicted with `theorem-package.md` listing the same premises
under `Proven`; both now read `Assumption` for these; the representation premise `A9` (added in the review-response revision) carries the same `Assumption` label.)

## Active Theorem Assumptions

| ID | Status | Assumption | Why it is present | If dropped |
| --- | --- | --- | --- | --- |
| A1 | Assumption | The body is quasi-static on the orbital timescale. | Needed to eliminate fast internal readout variables from the free-fall EFT. | Time-dependent internal modes can appear as explicit state variables. |
| A2 | Assumption | The minimal M4 sector is nearly spherical, nonspinning, and parity-even. | Keeps the operator search space inside the minimal free-fall sector requested for this theorem candidate. | Spin, parity-odd structures, or large asphericity add new operator families at the same order. |
| A3 | Assumption | A local worldline EFT exists after integrating out short-distance structure. | Supplies the operator classification language used by Theorem A. | Nonlocal kernels can survive directly in the effective action. |
| A4 | Assumption | There is no orbital-timescale internal state variable in the free-fall sector. | This is the sharp assumption that rules out `chi`-type hidden coordinates. | A `chi` state can carry body memory that is not reducible to instantaneous sensitivities. |
| A5 | Assumption | Couplings to the relevant external invariants are analytic near the reference background. | Needed to define a finite Taylor jet of sensitivity coordinates at fixed order. | Threshold or cusp behavior can evade a polynomial sensitivity jet. |
| A6 | Assumption | The object admits a self-bound equilibrium before external perturbations are applied. | Separates body formation from later passive coupling to external gravity. | Otherwise the theorem can mix equilibrium failure with observational coupling. |
| A7 | Assumption | The theorem is stated at fixed order in the operator counting rule of [`power-counting.md`](power-counting.md). | This prevents the proof target from silently becoming an all-orders closure claim. | Without fixed-order truncation, finite-dimensionality needs an extra closure theorem. |
| A8 | Assumption | The admitted primitive-family spectrum is locally finite below the fixed theorem cutoff: only finitely many primitive-family species have intrinsic weight `w \le \Delta_{\max}`. | This is the weakest hypothesis that replaces the older hidden finite-basis assumption at fixed order. | An infinite low-weight tower can generate infinitely many primitive directions before normal-form reduction starts. |
| A9 | Assumption | The external tidal field enters through a leading-Newtonian harmonic scalar potential: `E_{ij} = \partial_i\partial_j \Phi_ext` with `\nabla^2\Phi_ext = 0` (purely electric, leading order). | Fixes the tidal representation. It -- not tracelessness of a generic electric-Weyl tidal tensor -- makes `\nabla_k E_{ij}` totally symmetric (Schwarz) and trace-free on every index pair, i.e. an STF-3 octupole, which licenses the gradient-sector reduction. | A generic electric-Weyl tidal tensor restores the divergence and curl gradient pieces fixed by the Bianchi / gravitoelectromagnetic constraint equations (Danehkar 2022); the pre-correction generic-gradient model returns three quadratic gradient invariants. |

## Imported Exclusions From Earlier Work

| ID | Status | Imported statement | Source | Use in the new theorem-first repo |
| --- | --- | --- | --- | --- |
| E1 | Imported from prior work | Literal internal self-unobservability does not produce a viable compact self-bound equilibrium. | [`request2/REQUEST2_INTERNAL_STRUCTURE.md`](../request2/REQUEST2_INTERNAL_STRUCTURE.md) | Removes the old "change the stellar structure equations directly" branch from Theorem A. |
| E2 | Imported from prior work | COM self-subtraction does not create a new monopole force in the free-fall sector. | [`request1/REQUEST1_COM_DECOUPLING.md`](../request1/REQUEST1_COM_DECOUPLING.md) | Prevents the proof from smuggling in a self-force through coordinate bookkeeping. |
| E3 | Imported from prior work | The strongest previous global result is a provisional consistency scaffold, not a theorem. | [`request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md`](../request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md) | Justifies pivoting the new repo toward proof-or-loophole classification. |

## Non-Assumption For M4

- Status: Note. Finite invariant basis closure is not listed as an assumption in M4.
- Status: Note (corrected 2026-07-12; representation premise made explicit in the review response). Given the leading-Newtonian harmonic-scalar-potential representation (premise `A9`: `E_{ij} = \partial_i\partial_j\Phi`, `\nabla^2\Phi = 0`), the gradient-block collapse `mixedGradE2 = gradE2`, `divE2 = 0` is **not a separate assumption**: it then follows from the Schwarz total symmetry of `\nabla_k E_{ij} = \partial^3\Phi` (kinematics) and the same external-vacuum condition `\nabla^2\Phi = 0` already in force for `E` itself. What *is* a premise is the representation `A9` — the choice of a harmonic scalar potential rather than a generic electric-Weyl tidal tensor, now recorded in the active ledger above; for the generic Weyl gradient the divergence and curl pieces survive (Danehkar 2022). The pre-correction ledger instead listed the collapse identities as optional reductions that were declined, without recording the representation they rest on; that bookkeeping was internally inconsistent and is superseded — see [`reduction-rules.md`](reduction-rules.md) and [`../lemmas/07-gradient-sector-audit.md`](../lemmas/07-gradient-sector-audit.md).
- Status: Note. Contraction-level exhaustiveness for the exact current primitive blocks is attacked by explicit enumeration rather than by absence arguments.

## Non-Assumption For M5

- Status: Note. Linear survivor independence is now attacked directly by an exact polynomial rank check rather than inferred from non-reducibility.
- Status: Note. The exact-current-set theorem does not claim that the present primitive catalog already exhausts every physically admissible minimal-sector family.
- Status: Note. Excluding the magnetic tidal family `B_ij` is currently a working restriction of the exact-current-set theorem, not a derived completeness statement.

## Non-Assumption For M6

- Status: Note. No magnetic-ordering assumption is currently active that suppresses `B_ij` relative to `E_ij`.
- Status: Note. The repository does not currently assume `B_ij = 0` on the allowed backgrounds.
- Status: Note. Therefore any stronger electric-only minimal-sector theorem would require a new explicit assumption, recorded rather than implied.

## Explicit Rule For M7

- Status: Note. The mixed quartic `E/B` STF identity used to eliminate `EBEB` is now explicit in [`reduction-rules.md`](reduction-rules.md); it is not being used tacitly.
- Status: Note. No stronger magnetic-ordering salvage assumption is currently active.

## Non-Assumption For M8

- Status: Note. No scalar-family ordering assumption is currently active that suppresses a parity-even rank-0 external family `S`.
- Status: Note. No derivative-only or shift-symmetric scalar rule is being used tacitly.
- Status: Note. No scalar-background restriction such as `S = 0`, `S = const.`, or `\nabla_i S = 0` is currently assumed.

## Explicit Rule For M9

- Status: Note. The derivative-only scalar audit is a separate family-admission test, not a new global theorem assumption.
- Status: Note. In that audit, non-shift-symmetric total-derivative preimages such as `D_\tau(S X)` are not used to reduce operators, because bare `S` is excluded by construction.
- Status: Note. No cross-family orthogonality or alignment condition such as `Tr(EX) = 0`, `SE2 = 0`, or `DtS_E2 = 0` is assumed in the mixed-witness threshold audit.
- Status: Note. The mixed-witness threshold notes do not infer mixed-sector harmlessness from self-witness analysis alone; self and mixed survivors are checked separately.

## Non-Assumption For M12

- Status: Note. The repository does not currently assume that the audited family set `{R2, R0a, R0b, R1}` already exhausts the parity-even local MVP primitive-family envelope.
- Status: Note. No explicit assumption restricts admitted parity-even local primitive families to scalar and rank-2 classes only.
- Status: Note. Therefore local parity-even higher-rank tensor families remain open family-envelope targets unless they are audited directly or excluded by a new explicit assumption.

## Non-Assumption For M13

- Status: Note. No current MVP assumption says that every local parity-even vector block must be derivative-generated from an audited scalar or rank-2 family.
- Status: Note. Derivative-generated vectors such as `\nabla_i S` or `\nabla_j X_i{}^j` are assigned to their parent families and are not counted as new primitive-family admissions.
- Status: Note. No vector-family ordering, absorption rule, or vector-background restriction is currently active that would remove a genuine primitive family `V_i` from the theorem domain.

## Explicit Representative Choice For M14

- Status: Note. No current MVP assumption says that every local parity-even rank-3 tensor block must be derivative-generated from an audited lower-rank family.
- Status: Note. Derivative-generated rank-3 blocks such as `\nabla E`, `\nabla B`, `\nabla V`, or scalar-derivative descendants are assigned to their parent families and are not counted as new primitive-family admissions.
- Status: Note. The `Rodd+` audit uses a fully symmetric trace-free rank-3 representative explicitly, because trace descendants reduce to the already audited vector gate `R1`, while parity-odd or spin-like antisymmetric pieces are already outside the current MVP theorem domain by `A2`.
- Status: Note. No rank-3 ordering, absorption rule, or rank-3 background restriction is currently active that would remove a genuine primitive family `T_{ijk}` from the theorem domain.

## Explicit Representative Choice For M15

- Status: Note. No current MVP assumption says that every local parity-even rank-5 tensor block must be derivative-generated from an audited lower-rank family.
- Status: Note. Derivative-generated rank-5 blocks such as `\nabla Q`, `\nabla T`, higher-gradient vector descendants, or scalar descendants dressed only by derivatives are assigned to their parent families and are not counted as new primitive-family admissions.
- Status: Note. The `Rodd5+` audit uses a fully symmetric trace-free rank-5 representative explicitly, because trace descendants reduce to already audited lower-rank classes, while antisymmetric or spin-carrying pieces are already outside the current MVP theorem domain by `A2`.
- Status: Note. No rank-5 ordering, absorption rule, or rank-5 background restriction is currently active that would remove a genuine primitive family `U_{ijklm}` from the theorem domain.

## Non-Assumption For The High-Rank Exhaustiveness Patch

- Status: Note. The repository does not assume that the current manual rank-3, rank-4, or rank-5 survivor lists are already exhaustive.
- Status: Note. High-rank survivor exhaustiveness is now checked by explicit contraction generation in [`high-rank-audit-methodology.md`](high-rank-audit-methodology.md), [`../symbolic/high_rank_family_enumerator.py`](../symbolic/high_rank_family_enumerator.py), and [`../symbolic/high_rank_diff_report.py`](../symbolic/high_rank_diff_report.py).

## Explicit Rule For The Irreducible-Envelope Theorem

- Status: Note. Primitive-family distinctness is tested after linear `O(3)` irrep decomposition rather than before it.
- Status: Note. Trace descendants are not counted as new primitive-family admissions once their lower-rank scalar/vector/STF outputs are already present in the audited catalog.
- Status: Note. Epsilon-dual pseudo sectors are excluded from the current theorem domain by the parity-even nonspinning assumption `A2`.
- Status: Note. Even-dual mixed-symmetry sectors are reduced to ordinary lower-rank tensors and then decomposed again into STF plus trace-descendant pieces rather than counted as new primitive families.

## Explicit Rule For The Positive Collapse Bridge

- Status: Note. The scalar normal-form coordinates `Y^I` refer only to the reduced parity-even scalar operator basis at fixed order; they are not a placeholder for higher-multipole Wilson data.
- Status: Note. The higher-multipole Wilson coefficients `C_{A,\alpha}` remain separate from the monopole sensitivity jet and are not merged into the same coefficient family.
- Status: Note. Fixed-order operator finiteness uses `A7` and the sharpened local weight-spectrum condition `A8` together with positive intrinsic weights; it does not use analyticity `A5`.
- Status: Note. Monopole jet collapse uses analyticity `A5` only after the finite scalar normal-form basis `Y^I` has already been established.
- Status: Note. The sensitivity/Wilson split uses locality `A3` and the finite reduced operator set guaranteed by `A7` plus the sharpened `A8`; it is not being smuggled in as an automatic consequence of family-envelope closure alone.

## Explicit Rule For The A5 Boundary Stress Test

- Status: Note. The `A5` boundary-stress program drops analyticity only; it keeps locality `A3`, fixed-order truncation `A7`, and the sharpened local weight-spectrum condition `A8`.
- Status: Note. Nonanalytic monopole-response data are not merged with higher-multipole Wilson coefficients.
- Status: Note. Failure of the analytic jet does not by itself count as failure of irreducible family-envelope closure or fixed-order operator finiteness.

## Explicit Rule For The A3 Boundary Stress Test

- Status: Note. The `A3` boundary-stress program drops locality only after the finite primitive-family envelope, fixed-order local operator finiteness, and local scalar normal-form quotient have already been fixed.
- Status: Note. Instantaneous analyticity is kept on the local-coordinate side whenever possible; the intended stress point is hereditary memory, not a mixed `A3+A5` escape.
- Status: Note. Exponential or finite multi-exponential kernels are not counted as sharp `A3` counterexamples if they are finitely Markovianizable and can be rewritten as finite auxiliary-state extensions.
- Status: Note. Genuine hereditary kernel or spectral data are kept separate from higher-multipole Wilson coefficients.

## Explicit Rule For The A4 Boundary Stress Test

- Status: Note. The `A4` boundary-stress program drops only the no-orbital-timescale-state assumption while keeping locality `A3`, analyticity `A5`, fixed-order truncation `A7`, and the sharpened local weight-spectrum condition `A8`.

## Explicit Rule For The A8 Sharpness Theorem

- Status: Note. Finite total admitted primitive-family content implies the sharpened `A8` condition, but it is stronger than necessary at fixed order.
- Status: Note. The fixed-order positive theorem only needs local family-spectrum finiteness below `\Delta_{\max}`, not a globally finite family catalog.
- Status: Note. Infinite admitted family count is harmless by itself if only finitely many species survive with intrinsic weight `w \le \Delta_{\max}`.
- Status: Note. The sharp `A8` failure is the infinite low-weight tower in which infinitely many species survive below the fixed cutoff and produce infinitely many pre-reduction scalar witnesses.
- Status: Note. Adiabatically eliminable or algebraically slaved local state variables do not count as sharp `A4` counterexamples, because they collapse back into a `Y`-only monopole model.
- Status: Note. The sharp `A4` counterexample must keep finitely many genuinely dynamical local state variables explicit and distinct from the scalar normal-form coordinates `Y^I`.
- Status: Note. Explicit local state variables `\chi^a` are not merged with higher-multipole Wilson coefficients.
- Status: Note. The surviving positive salvage statement in the finite-state branch is a finite-dimensional local state-space theorem; an additional coefficient-counting truncation in `(Y^I,\chi^a)` is not smuggled in automatically.
