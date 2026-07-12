# Release Note

- Status: Note. This note is the short publication and handoff summary of the frozen free-fall theorem package.
- Status: Proven. It restates only closed results already established elsewhere in the repository.

## Theorem Domain

- Status: Proven. The closed theorem domain is the parity-even, nonspinning, local MVP free-fall sector at fixed order `Delta <= 4`.
- Status: Proven. The active domain assumptions are the `A1`-`A9` ledger assumptions in [`assumptions-ledger.md`](assumptions-ledger.md): quasi-static regime, nearly spherical body, local worldline EFT, no orbital-timescale internal state variable, analytic monopole response, fixed-order truncation, local weight-spectrum finiteness below the cutoff, and a leading-Newtonian harmonic-scalar-potential tidal representation.
- Status: Proven. Within that domain, the irreducible primitive-family envelope closes on the audited scalar, vector, rank-2 STF, and genuine rank-`L >= 3` STF classes.

## Main Positive Theorem

- Status: Proven. Once the irreducible scalar/vector/STF family envelope is fixed and the admitted primitive-family spectrum is locally finite below the cutoff, the candidate parity-even local scalar operator space at `Delta <= 4` is finite.
- Status: Proven. After the explicit total-derivative, lower-order-EOM, algebraic, and linear-dependence reductions are imposed, the reduced scalar operator space is finite-dimensional and admits a finite normal-form basis `Y^I`.
- Status: Proven. Under locality `A3`, no-state `A4`, and analyticity `A5`, the monopole response collapses to a finite Taylor jet in those finitely many scalar coordinates.
- Status: Proven. The remaining higher-multipole sector is carried by finitely many Wilson coefficients that remain separate from the monopole sensitivity data.
- Status: Proven. Therefore the positive finite-family collapse theorem closes inside the stated theorem domain.

## Negative Uniqueness No-Go

- Status: Proven. The repo does not prove minimal-sector uniqueness.
- Status: Proven. Across the audited unsuppressed primitive-family classes, admission of a genuinely new primitive family yields a low-order witness and therefore obstructs any theorem that tries to identify a unique physically justified minimal sector without further suppression assumptions.
- Status: Proven. The negative uniqueness result is a class-limited family-admission no-go, not a refutation of the positive finite-family collapse theorem.

## Sharp Boundary Escapes

| Boundary | Status | Smallest explicit counterexample | Exact theorem layer broken | Replacement bookkeeping |
| --- | --- | --- | --- | --- |
| `A5` dropped | Proven | Smooth-flat local monopole model `m_A(Y)=m_0+\alpha e^{-1/Y^2}\Theta(Y)` | Finite analytic Taylor-jet step | Non-Taylor monopole germ data |
| `A3` dropped | Proven | One-coordinate causal power-law hereditary kernel | Local reduction to an instantaneous monopole function | Memory kernel or spectral data |
| `A4` dropped | Proven | One-state local analytic `chi` model | Y-only monopole reduction | Finite local state-space data `(Y^I, chi^a)` |
| `A8` dropped | Proven | Infinite low-weight STF tower | Candidate operator-space finiteness before reduction | No finite pre-reduction primitive-family catalog remains |

## Explicit Non-Claims

- Status: Proven. No theorem is claimed here for parity-odd, spinning, clock, nonlocal, or orbital-timescale-state sectors.
- Status: Proven. No universal mixed-pattern theorem is claimed for all higher-rank tensor families.
- Status: Proven. No all-orders closure theorem is claimed beyond the fixed cutoff.
- Status: Proven. No empirical weak-field estimator, strong-field runtime, or TOA-performance claim is made in this theorem repo.

## Freeze Status

- Status: Proven. The mathematical content of this repository is now frozen unless a direct contradiction is found during cleanup or later review.
- Status: Proven. Further work belongs to empirical branches, runtime branches, or paper writing, not to new theorem expansion inside this repo.
