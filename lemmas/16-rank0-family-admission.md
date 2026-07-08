# Lemma 16: Rank-0 Family Admission

## Family Class

- Status: Proven. Consider a scalar-like rank-0 primitive family in the free-fall MVP sector.
- Status: Proven. Two audited subclasses matter separately:
  `R0a` unsuppressed bare-scalar admission with primitive `S`,
  and `R0b` derivative-only or shift-symmetric admission with primitives `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S` but no bare `S`.

## Subclass `R0a`: Unsuppressed Bare Scalar

- Status: Proven. The operator `S` itself is a parity-even weight-`1` scalar.
- Status: Proven. Under the stated reduction rules, it is not removed by total derivatives, lower-order EOM, or any listed algebraic identity.
- Status: Proven. Therefore `S` is the smallest explicit witness obstructing minimal-sector uniqueness for the unsuppressed scalar-family subclass.

## Subclass `R0b`: Derivative-Only Or Shift-Symmetric Scalar

- Status: Proven. Removing bare `S` by hand does eliminate the weight-`1` witness.
- Status: Proven. [`14-derivative-only-scalar-audit.md`](14-derivative-only-scalar-audit.md) shows that this does not make the scalar family harmless.
- Status: Proven. The first new survivors now appear at weight `4`.
- Status: Proven. The full smallest-weight audited witness set is
  `\{DtS_B2,\ dotS2,\ DtS_E2,\ divEGradS,\ gradS2\}`.
- Status: Proven. A canonical explicit derivative-only scalar witness is `dotS2`.

## Consequence

- Status: Proven. Scalar-family admission obstructs minimal-sector uniqueness both with and without bare `S`.
- Status: Proven. Shift symmetry changes the first witness but does not rescue uniqueness by itself.
- Status: Proven. Therefore the scalar-family no-go is part of the family-class structure, not just a side audit for one special model.

## Boundary

- Status: Note. This lemma does not say that finite-family fixed-order collapse fails.
- Status: Proven. It only says that rank-0 family admission is not harmless without additional explicit suppression assumptions.
