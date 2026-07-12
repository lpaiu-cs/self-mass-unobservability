# Lemma 09: Survivor Independence At `Delta<=4`

## Question

- Status: Proven. Are the corrected five surviving `Delta<=4` scalars linearly independent modulo the currently allowed reduction rules?

## Independence Target

- Status: Proven. The target list is

```math
\{E2,\ E3,\ E2^2,\ dotE2,\ gradE2\}.
```

- Status: Note. The M4-era seven-element target additionally listed `divE2` and `mixedGradE2`; those are not independent operators once the gradient block is modeled with its actual kinematics (STF-3: Schwarz total symmetry + vacuum trace-free) — see [`07-gradient-sector-audit.md`](07-gradient-sector-audit.md).

- Status: Note. This lemma asks for linear operator independence over constant coefficients after quotienting by the allowed total-derivative, lower-order-EOM, and Cayley-Hamilton reductions.
- Status: Note. This lemma does not claim algebraic functional independence.
- Status: Proven. In particular, `E2^2` is algebraically dependent on `E2` as a product, but it can still be linearly independent as a separate weight-4 operator.

## Exact Rank Check

- Status: Proven. [`../symbolic/survivor_rank_check.py`](../symbolic/survivor_rank_check.py) constructs explicit STF component parametrizations for `E_ij` and `D_tau E_ij` (5 components each) and the explicit 7-parameter STF-3 parametrization for `nabla_k E_ij`, and expands every survivor as a polynomial in independent tensor components.
- Status: Proven. The exact polynomial coefficient matrix has total rank `5`.
- Status: Proven. The `E`-sector subset `{E2, E3, E2^2}` has rank `3`.
- Status: Proven. The single `DtE` survivor `{dotE2}` has rank `1`.
- Status: Proven. The gradient sector `{gradE2}` has rank `1`; on the STF-3 parametrization `divE2` vanishes identically and `mixedGradE2` coincides with `gradE2`.

## Why The Rank Result Is Enough

- Status: Proven. A vanishing linear dependence

```math
c_1 E2 + c_2 E3 + c_3 E2^2 + c_4 dotE2 + c_5 gradE2 \equiv 0
```

would force the coefficient vector `(c_1,\dots,c_5)` into the nullspace of that exact polynomial coefficient matrix.
- Status: Proven. Since the computed rank is `5`, that nullspace is trivial.
- Status: Proven. Therefore no nontrivial linear dependence survives under the currently allowed rules for the exact current primitive set.

## Result

- Status: Proven. The corrected five-scalar `Delta<=4` normal-form list is linearly independent as an operator basis for the exact current primitive set.
- Status: Proven. The current remaining bottleneck is therefore not survivor independence.
- Status: Conjectural. The main remaining burden is primitive-set adequacy: whether the exact current primitive set matches a physically justified minimal sector.
