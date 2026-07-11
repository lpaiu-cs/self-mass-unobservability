# Lemma 59: Nonanalytic Sensitivity Replacement

## Statement

- Status: Proven. Once the smooth flat local counterexample of [`58-nonanalytic-jet-failure.md`](58-nonanalytic-jet-failure.md) is admitted, the finite Taylor sensitivity jet of [`55-monopole-jet-collapse.md`](55-monopole-jet-collapse.md) no longer exists.
- Status: Proven. What replaces it is non-Taylor monopole germ data.
- Status: Proven. That replacement data remains separate from the higher-multipole Wilson coefficients.

## Explicit Replacement Data

- Status: Proven. For the explicit smooth flat model class

```math
m_A(Y)=m_0+\alpha\,\phi_{\mathrm{flat}}(Y),
```

the minimal replacement bookkeeping is the tuple

```math
(m_0,\ \alpha,\ \phi_{\mathrm{flat}}).
```

- Status: Proven. For the threshold model

```math
m_A(Y)=m_0+\alpha\,\Theta(Y-Y_c)(Y-Y_c)^{1/2},
```

the minimal replacement bookkeeping is the tuple

```math
(m_0,\ \alpha,\ Y_c,\ p=1/2,\ \text{branch label}).
```

## Interpretation

- Status: Note. These tuples are not finite Taylor sensitivity coordinates.
- Status: Proven. They are finite generalized branch data only because the model class has been fixed explicitly in advance.
- Status: Proven. If the model class is not fixed, locality plus finite operator closure alone do not force a universal finite replacement; one may need the full nonanalytic function germ.

## Separation From Wilson Data

- Status: Proven. The nonanalytic monopole germ data belong to the monopole response sector.
- Status: Proven. Wilson coefficients still multiply the residual higher-multipole operators in the action and are not absorbed into the nonanalytic branch data.
- Status: Proven. Therefore the counterexample breaks the analytic sensitivity-jet theorem, not the operator/Wilson bookkeeping split at the action level.
