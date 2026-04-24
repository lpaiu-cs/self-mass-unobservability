# Lemma 09 Draft: Survivor Independence At `Delta<=4`

## Question

- Status: Conjectural. Are the corrected seven surviving `Delta<=4` scalars linearly independent modulo the currently allowed reduction rules?

## Independence Target

- Status: Proven. The target list is

```math
\{E2,\ E3,\ E2^2,\ dotE2,\ gradE2,\ divE2,\ mixedGradE2\}.
```

- Status: Proven. This lemma asks for linear operator independence over constant coefficients after quotienting by the allowed total-derivative, lower-order-EOM, and Cayley-Hamilton reductions.
- Status: Proven. This lemma does not claim algebraic functional independence.
- Status: Proven. In particular, `E2^2` is algebraically dependent on `E2` as a product, but it can still be linearly independent as a separate weight-4 operator.

## Exact Rank Check

- Status: Proven. [`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py) constructs explicit STF component parametrizations for `E_ij`, `D_tau E_ij`, and `nabla_k E_ij` and expands every survivor as a polynomial in independent tensor components.
- Status: Proven. The exact polynomial coefficient matrix has total rank `7`.
- Status: Proven. The `E`-sector subset `{E2, E3, E2^2}` has rank `3`.
- Status: Proven. The single `DtE` survivor `{dotE2}` has rank `1`.
- Status: Proven. The gradient-sector subset `{gradE2, divE2, mixedGradE2}` has rank `3`.

## Why The Rank Result Is Enough

- Status: Proven. A vanishing linear dependence

```math
c_1 E2 + c_2 E3 + c_3 E2^2 + c_4 dotE2 + c_5 gradE2 + c_6 divE2 + c_7 mixedGradE2 \equiv 0
```

would force the coefficient vector `(c_1,\dots,c_7)` into the nullspace of that exact polynomial coefficient matrix.
- Status: Proven. Since the computed rank is `7`, that nullspace is trivial.
- Status: Proven. Therefore no nontrivial linear dependence survives under the currently allowed rules for the exact current primitive set.

## Result

- Status: Proven. The corrected seven-scalar `Delta<=4` normal-form list is linearly independent as an operator basis for the exact current primitive set.
- Status: Proven. The current remaining bottleneck is therefore not survivor independence.
- Status: Conjectural. The main remaining burden is primitive-set adequacy: whether the exact current primitive set matches a physically justified minimal sector.
