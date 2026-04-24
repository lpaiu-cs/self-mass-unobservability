# Lemma 07 Draft: Gradient-Sector Audit

## Question

- Status: Conjectural. Is `gradE2` the only surviving `Delta<=4` gradient scalar, or do alternative contractions survive under the currently allowed reduction rules?

## Enumerated Gradient Contractions

- Status: Proven. The exhaustive contraction generator [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) finds exactly three parity-even scalar contraction classes built from two `GradE` blocks.

| Operator | Classification | Status | Meaning |
| --- | --- | --- | --- |
| `gradE2` | Surviving candidate | Conjectural | `(\nabla_k E_{ij})(\nabla^k E^{ij})` |
| `divE2` | Surviving candidate | Conjectural | `(\nabla_i E^{ij})(\nabla^k E_{kj})` |
| `mixedGradE2` | Surviving candidate | Conjectural | `(\nabla_k E_{ij})(\nabla^i E^{k j})` |

## Audit Result

- Status: Proven. `gradE2` is not the only surviving gradient invariant under the currently allowed reduction rules.
- Status: Proven. No currently allowed total-derivative, lower-order EOM, or Cayley-Hamilton reduction removes `divE2`.
- Status: Proven. No currently allowed total-derivative, lower-order EOM, or Cayley-Hamilton reduction removes `mixedGradE2`.
- Status: Proven. Therefore the minimal corrected `Delta<=4` gradient sector is at least three-dimensional.

## Why The Extra Invariants Survive

- Status: Proven. `divE2` and `mixedGradE2` contain no explicit acceleration factor, so lower-order worldline EOM reduction does not apply.
- Status: Proven. They are not worldline total derivatives.
- Status: Proven. They are not consequences of the traceless `3x3` Cayley-Hamilton identity, which acts only on purely algebraic quartic powers of `E`.

## What Would Be Needed To Reduce Further

- Status: Conjectural. Eliminating `divE2` would require a new explicit identity such as a transversality/vacuum condition on `\nabla_i E^{ij}`.
- Status: Conjectural. Eliminating `mixedGradE2` would require a stronger gradient-sector identity than any currently assumed rule.
- Status: Proven. No such identity is used in the present `Delta<=4` audit.
