# Lemma 07: Gradient-Sector Audit (Corrected: Kinematic STF-3 Collapse)

## Question

- Status: Proven. Is `gradE2` the only surviving `Delta<=4` gradient scalar, or do alternative contractions survive under the allowed reduction rules?

## Corrected Result

- Status: Proven. `gradE2` is the only surviving `Delta<=4` gradient invariant. The gradient sector is one-dimensional.
- Status: Proven. The gradient block is `\nabla_k E_{ij} = \partial_k\partial_i\partial_j \Phi_{\rm ext}`: it is **totally symmetric** in `(k,i,j)` by equality of mixed partials (Schwarz — a kinematic identity of the potential, not an optional reduction rule), and **trace-free on every index pair** in the external vacuum (`\nabla^2\Phi_{\rm ext} = 0` — the same condition already used to make `E_{ij}` itself trace-free). It is therefore an STF-3 octupole with 7 components, not a generic `(STF\text{-}2)\otimes\text{vector}` object with 15.
- Status: Assumption. This audit presumes the leading-Newtonian harmonic-scalar-potential representation `E_{ij} = \partial_i\partial_j\Phi_{\rm ext}`, `\nabla^2\Phi_{\rm ext} = 0` — theorem-domain premise `A9` in [`../docs/assumptions-ledger.md`](../docs/assumptions-ledger.md). The STF-3 collapse is a consequence of that representation, not of tracelessness of a generic electric-Weyl tidal tensor: for the electric part of the Weyl tensor in the fully relativistic theory the covariant gradient `\nabla_k E_{ij}` additionally carries divergence and curl pieces fixed by the Bianchi / gravitoelectromagnetic constraint equations (Danehkar 2022), which the harmonic-scalar-potential restriction removes. Relaxing `A9` therefore restores the generic gradient sector; the Newtonian sourced-background case (`\nabla^2\Phi \ne 0`) below is the sub-case in which only the divergence piece returns.
- Status: Proven. Consequently, of the three contraction classes of the generic model:

| Operator | Meaning | Fate under the STF-3 kinematics |
| --- | --- | --- |
| `gradE2` | `(\nabla_k E_{ij})(\nabla^k E^{ij})` | survives — the unique gradient invariant |
| `divE2` | `(\nabla_i E^{ij})(\nabla^k E_{kj})` | `= 0` (vacuum trace) |
| `mixedGradE2` | `(\nabla_k E_{ij})(\nabla^i E^{k j})` | `= gradE2` (Schwarz total symmetry) |

- Status: Proven. Both identities are verified exactly: symbolically on the explicit 7-parameter STF-3 parametrization ([`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py)) and numerically to machine precision over random STF-3 samples ([`../verification/verify_survivors.py`](../verification/verify_survivors.py)); the corrected contraction enumerator ([`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py)) and an independent brute-force signature sweep ([`../verification/verify_completeness.py`](../verification/verify_completeness.py)) both find exactly one `(GradE, GradE)` class.

## History Of The Error (Recorded, Not Erased)

- Status: Note. The M4-era audit modeled `GradE` as symmetric/trace-free **only on the `E` indices** with a free gradient index (`sym_groups=((1,2),)` in the enumerator's block table). Under that generic model, three contraction classes are genuinely independent (the generic-model control in `verify_survivors.py` still reproduces rank 3), and the M4 audit accordingly promoted `divE2` and `mixedGradE2` to survivors, moving the electric normal-form dimension from the pre-M4 5 to 7.
- Status: Note. That step conflated two distinct facts: (i) *no allowed reduction rule* removes `divE2`/`mixedGradE2` — true under the generic model; and (ii) the objects are independent — false for the actual gradient of an actual tidal field, because Schwarz symmetry is kinematics (there is no "option" to decline it) and the vacuum condition was already in force for `E` itself. Declining transversality for `\nabla E` while imposing it for `E` was internally inconsistent.
- Status: Note. The correction was triggered by an external adversarial review (2026-07-12) and verified against the repository sources; the corrected electric survivor set returns to the pre-M4 five-element basis, now with contraction-level exhaustiveness *and* the correct kinematics both in force.

## Why The Remaining Invariant Survives

- Status: Proven. `gradE2` contains no explicit acceleration factor, is not a worldline total derivative, is not removed by the STF-3 trace conditions, and is not a consequence of the Cayley–Hamilton identity. It is the unique quadratic invariant of an STF-3 tensor (one `SO(3)` irrep ⇒ one quadratic invariant).

## What Would Change The Count

- Status: Note. Off the vacuum background (matter at the worldline, `\nabla^2\Phi \ne 0`), the trace `\nabla_i E^{ij}` is sourced and `divE2` returns as an independent invariant tied to the local source density; `mixedGradE2 = gradE2` continues to hold by Schwarz regardless. Any such extension is a change of theorem domain (a sourced background), not a repair of the vacuum count.
