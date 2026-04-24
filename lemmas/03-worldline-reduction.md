# Lemma 03 Draft: Generic Worldline Reduction

## Draft Statement

Status: Conjectural. At fixed order in a local derivative expansion, the free-fall sector of a quasi-static compact body can be organized into a monopole response function plus higher-multipole Wilson operators.

## Generic Action Template

Status: Conjectural. The intended worldline form is

```math
S_A
=
\int d\tau \left[
-m_A(Y)
-\frac{1}{2} C_A^{ij} {\cal E}_{ij}
-\frac{1}{6} C_A^{ijk} \nabla_k {\cal E}_{ij}
+ \cdots
\right],
```

where:

- Status: Conjectural. `Y^I` is a finite basis of scalarized local invariants built from the external environment at the chosen order.
- Status: Conjectural. `m_A(Y)` contains monopole body dependence.
- Status: Conjectural. `C_A^{ij}`, `C_A^{ijk}`, and higher tensors are higher-multipole Wilson data.

## Reduction Logic

1. Status: Conjectural. Symmetry removes operators incompatible with near-sphericity and the no-spin MVP.
2. Status: Conjectural. The no-internal-state assumption removes explicit orbital-timescale coordinates from the free-fall action.
3. Status: Conjectural. Analyticity turns `m_A(Y)` into a Taylor jet whose coefficients are the sensitivity coordinates.

## Current Gap

- Status: Conjectural. The repository still needs a complete bookkeeping rule for the invariant basis `Y^I` at each fixed order; until that is written sharply, this note remains a reduction draft rather than a finished lemma.
