# Proposition 05: Abstract Finite Candidate-Set Generation At Fixed Order

## Statement

- Status: Proven. Fix the counting rule of [`../docs/power-counting.md`](../docs/power-counting.md), a truncation order `\Delta \le \Delta_{\max}`, and finite external field content.
- Status: Proven. Once a candidate primitive catalog is fixed, the abstract set of parity-even scalar monomials and parity-even higher-rank tensor contractions with total weight `\le \Delta_{\max}` is finite.
- Status: Proven. This proposition establishes abstract candidate-set finiteness only; it does not establish `normal-form completeness modulo total derivatives and lower-order equations of motion`.

## Proof Of Abstract Finiteness

1. Status: Proven. Finite external field content gives only finitely many primitive tensor species `X_a`.
2. Status: Proven. Positive weights and a fixed cutoff `\Delta_{\max}` imply that only finitely many decorated building blocks `D_\tau^p \nabla^q X_a` satisfy `w(X_a) + p + q \le \Delta_{\max}`.
3. Status: Proven. Each admissible monomial has exponent vector `n = (n_1, \dots, n_r)` constrained by

```math
\sum_{i=1}^{r} n_i w_i \le \Delta_{\max},
```

so each `n_i` is bounded above by `\Delta_{\max} / w_i`.
4. Status: Proven. Therefore only finitely many exponent vectors occur.
5. Status: Proven. Each decorated block has finite tensor rank, the number of factors is bounded at fixed order, and therefore the number of parity-even complete contractions is finite.

## What This Proposition Does Not Yet Prove

- Status: Conjectural. This proposition does not prove that the chosen primitive catalog is complete for the physical EFT.
- Status: Conjectural. The remaining burden lives in [`06-normal-form-completeness-delta4.md`](06-normal-form-completeness-delta4.md): `normal-form completeness modulo total derivatives and lower-order equations of motion`.

## Minimal-Sector Check In Code

- Status: Proven. [`../symbolic/enumerate_basis.py`](../symbolic/enumerate_basis.py) enumerates the current candidate scalar monomials for the sample `\Delta_{\max} = 4` catalog.
- Status: Proven. The code check demonstrates finite abstract candidate-set size at the chosen order; it does not by itself prove normal-form completeness.
