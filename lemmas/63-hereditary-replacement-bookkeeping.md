# Lemma 63: Hereditary Replacement Bookkeeping

## Statement

- Status: Proven. Once the genuine hereditary power-law counterexample of [`62-hereditary-collapse-failure.md`](62-hereditary-collapse-failure.md) is admitted, the finite local sensitivity jet of [`55-monopole-jet-collapse.md`](55-monopole-jet-collapse.md) no longer exists.
- Status: Proven. What replaces it is hereditary memory data.
- Status: Proven. That hereditary data remains separate from higher-multipole Wilson coefficients.

## Control Case: Finite-State Memory

- Status: Proven. For the exponential-memory control, the minimal replacement bookkeeping is finite-state data:

```math
(m_0,\ \alpha,\ \lambda_{\exp},\ T_h,\ \chi_A(\tau_0)),
```

with `\chi_A` the auxiliary local memory state.
- Status: Note. This is not a pure hereditary replacement because it can be absorbed into an explicit local state extension.

## Genuine Hereditary Case

- Status: Proven. For the one-coordinate power-law model class, the minimal fixed-class bookkeeping is

```math
(m_0,\ \alpha,\ \lambda_{\gamma},\ \gamma,\ T_h,\ \text{kernel branch label}).
```

- Status: Proven. If the model class is not fixed in advance, the natural replacement data are the full causal kernel `K(s)` or, equivalently when available, a continuous relaxation-spectrum density `\rho(\mu)` satisfying

```math
K(s)=\int_0^{\infty} d\mu\, \rho(\mu)e^{-\mu s}.
```

- Status: Note. This is generally not a finite local sensitivity jet and not a finite auxiliary-state list.

## Separation From Wilson Data

- Status: Proven. The hereditary kernel data belong to the monopole-response sector.
- Status: Proven. Wilson coefficients still multiply residual higher-multipole local operators and are not absorbed into the hereditary memory data.
- Status: Proven. Therefore the genuine hereditary counterexample breaks the local sensitivity-jet side of the positive theorem, not the distinct Wilson bookkeeping sector.
