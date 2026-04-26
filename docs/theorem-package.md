# Theorem Package

- Status: Proven. This is the authoritative front-door theorem document for the free-fall theorem repo.
- Status: Proven. It is paper-facing, non-historical, and limited to the currently closed theorem and counterexample statements already established in the repository.
- Status: Proven. The mathematical content of this repository is now frozen unless a direct contradiction is found during cleanup or later review.
- Status: Proven. Further work belongs to empirical branches, runtime branches, or paper writing rather than new theorem expansion inside this repo.

## Theorem Domain

- Status: Proven. The closed positive theorem domain is the free-fall sector only.
- Status: Proven. The domain assumptions are the active `A1`-`A8` ledger assumptions in [`assumptions-ledger.md`](assumptions-ledger.md): quasi-static regime, nearly spherical and nonspinning parity-even sector, local worldline EFT, no orbital-timescale internal state variable, analytic monopole response, self-bound equilibrium, fixed operator cutoff, and local weight-spectrum finiteness below that cutoff.
- Status: Proven. Within that domain, the irreducible primitive-family envelope closes on the audited scalar, vector, rank-2 STF, and genuine rank-`L >= 3` STF classes.

## Main Positive Theorem

- Status: Proven. Inside the stated theorem domain at fixed order `Delta <= 4`, the admissible parity-even local free-fall scalar operator space is finite.
- Status: Proven. After the explicit reduction rules already recorded in the repo are imposed, the reduced scalar operator space is finite-dimensional and admits a finite normal-form basis `Y^I`.
- Status: Proven. Under locality `A3` and analyticity `A5`, the monopole response collapses to a finite Taylor jet in those finitely many scalar coordinates.
- Status: Proven. The residual higher-multipole sector is then carried by finitely many Wilson coefficients that remain separate from the monopole sensitivity coefficients.
- Status: Proven. Therefore the positive finite-family collapse theorem closes inside the current theorem domain.

## Negative Uniqueness No-Go

- Status: Proven. Minimal-sector uniqueness is not the positive theorem.
- Status: Proven. Across the audited unsuppressed family classes, admission of a genuine new primitive family yields a low-order witness and therefore obstructs any theorem that tries to identify a unique physically justified minimal sector without further suppression assumptions.
- Status: Proven. The negative branch is therefore a class-limited family-admission no-go for minimal-sector uniqueness, not a refutation of fixed-order collapse.

## Sharp Boundary Escapes

| Boundary | Status | Smallest explicit counterexample | Exact theorem layer broken | Replacement bookkeeping |
| --- | --- | --- | --- | --- |
| `A5` dropped | Proven | Smooth-flat local monopole model `m_A(Y)=m_0+\alpha e^{-1/Y^2}\Theta(Y)` | Analytic monopole Taylor-jet step | Non-Taylor monopole germ data |
| `A3` dropped | Proven | One-coordinate causal power-law kernel | Reduction to a local monopole function of instantaneous normal-form coordinates | Kernel or spectral memory data |
| `A4` dropped | Proven | One-state local analytic `chi` model | Y-only monopole reduction | Finite local state-space data `(Y^I, chi^a)` |
| `A8` dropped | Proven | Infinite low-weight STF tower | Candidate operator-space finiteness before reduction | No finite pre-reduction catalog remains |

## What Is Not Claimed

- Status: Proven. No theorem is claimed here for parity-odd, spinning, nonlocal, or orbital-timescale-state sectors.
- Status: Proven. No uniqueness theorem is claimed for a physically privileged minimal sector.
- Status: Proven. No all-orders closure theorem is claimed beyond the fixed cutoff.
- Status: Proven. No theorem is claimed here for clock observables.
- Status: Proven. Raw survivor counts, nullity counts, and rank-specific bookkeeping are not theorem statements unless explicitly promoted elsewhere.
