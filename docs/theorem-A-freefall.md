# Theorem A: Free-Fall Collapse Backbone

- Status: Proven. This note is the compressed technical backbone of the free-fall theorem package.
- Status: Proven. The authoritative front-door statement is [`theorem-package.md`](theorem-package.md).

## Closed Theorem Domain

- Status: Proven. Work in the parity-even, nonspinning, local MVP free-fall sector only.
- Status: Proven. The active theorem assumptions are `A1`-`A8` in [`assumptions-ledger.md`](assumptions-ledger.md): quasi-static regime, local worldline EFT, no orbital-timescale internal state variable, analytic monopole response, fixed order `Delta <= 4`, and local weight-spectrum finiteness below the cutoff.
- Status: Proven. Within that domain, the irreducible primitive-family envelope closes on the audited scalar, vector, rank-2 STF, and genuine rank-`L >= 3` STF classes.

## Positive Theorem Statement

- Status: Proven. Once the irreducible primitive-family envelope is fixed and the admitted primitive-family spectrum is locally finite below `Delta <= 4`, the candidate parity-even local scalar operator space is finite before reduction.
- Status: Proven. After the explicit total-derivative, lower-order-EOM, algebraic, and extracted linear-dependence relations are imposed, the reduced scalar operator space is finite-dimensional and admits a finite normal-form basis `Y^I`.
- Status: Proven. Under locality `A3` and analyticity `A5`, the monopole response collapses to a finite Taylor jet in those finitely many scalar coordinates.
- Status: Proven. The remaining higher-multipole sector is carried by finitely many Wilson coefficients `C_{A,alpha}` that remain separate from the monopole sensitivity coefficients.
- Status: Proven. Therefore the positive finite-family collapse theorem closes inside the current theorem domain.

## Negative Uniqueness No-Go

- Status: Proven. Minimal-sector uniqueness is a separate claim and it fails across the audited unsuppressed family classes.
- Status: Proven. The audited no-go covers the rank-2 STF special class, the bare scalar class, the derivative-only scalar class, the genuine vector class, and the genuine higher-rank STF branch.
- Status: Proven. The higher-rank STF branch supports a universal self-witness threshold theorem `w_Y >= 3` for genuine parity-even STF primitive families with rank `L >= 3`, while the stronger mixed-pattern theorem fails at audited rank `L = 4`.
- Status: Proven. The negative uniqueness no-go therefore coexists with the closed positive finite-family collapse theorem.

## Known Boundary Escapes

| Boundary | Status | Smallest explicit counterexample | Exact broken layer |
| --- | --- | --- | --- |
| `A5` dropped | Proven | Smooth-flat local monopole model `m_A(Y)=m_0+\alpha e^{-1/Y^2}\Theta(Y)` | Analytic monopole Taylor-jet collapse |
| `A3` dropped | Proven | One-coordinate causal power-law hereditary kernel | Local reduction to a monopole function of instantaneous normal-form coordinates |
| `A4` dropped | Proven | One-state local analytic `chi` model | Y-only monopole reduction, with a separate finite state-augmented salvage theorem surviving |
| `A8` dropped | Proven | Infinite low-weight STF tower | Candidate operator-space finiteness before reduction |

## Supporting Backbone

- Status: Proven. The positive theorem is detailed in [`finite-family-collapse-theorem.md`](finite-family-collapse-theorem.md) and frozen in [`collapse-bridge-status.md`](collapse-bridge-status.md).
- Status: Proven. The negative branch is detailed in [`family-admission-theorem.md`](family-admission-theorem.md).
- Status: Proven. The irreducible family-envelope closure is detailed in [`irreducible-envelope-theorem.md`](irreducible-envelope-theorem.md).
- Status: Proven. The higher-rank STF threshold layer is detailed in [`stf-self-witness-theorem.md`](stf-self-witness-theorem.md).
- Status: Proven. The exact assumption-drop failures are summarized in [`boundary-escape-map.md`](boundary-escape-map.md).
