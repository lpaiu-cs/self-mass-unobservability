# Failure Ledger

- Status: Proven. Inside the stated theorem domain `A1`-`A8`, the positive finite-family collapse theorem is closed.
- Status: Note. This file is no longer a to-do list. It is the sharp boundary-risk register for exact assumption-drop failures and scope escapes.
- Status: Proven. The older milestone-style ledger is preserved in [`archive/failure-ledger-history.md`](archive/failure-ledger-history.md).

## Theorem-Domain Assumptions

- Status: Proven. The closed positive theorem holds only inside the stated free-fall theorem domain recorded in [`assumptions-ledger.md`](assumptions-ledger.md).
- Status: Proven. The assumption-drop rows below are not open tasks. They are the exact places where the closed theorem stops applying or changes form.

## Exact Boundary Escapes

| ID | Status | Boundary drop | Smallest explicit model | Exact theorem layer broken | Replacement bookkeeping |
| --- | --- | --- | --- | --- | --- |
| F-A5 | Proven | `A5` dropped, locality kept | Smooth-flat local monopole model `m_A(Y)=m_0+\alpha e^{-1/Y^2}\Theta(Y)` | Finite analytic Taylor jet in Lemma 55 | Non-Taylor monopole germ data |
| F-A3 | Proven | `A3` dropped, analyticity kept on the instantaneous side | One-coordinate causal power-law hereditary kernel | Local reduction to a monopole function of instantaneous normal-form coordinates, so the local forms of Lemmas 55 and 56 fail | Memory kernel or spectral data |
| F-A4 | Proven | `A4` dropped, locality and analyticity kept | One-state local analytic `chi` model | Y-only monopole reduction, so the Y-only forms of Lemmas 55 and 56 fail | Finite local state-space data `(Y^I, chi^a)` and state-evolution parameters |
| F-A4+ | Proven | `A4` dropped but finitely many local states retained explicitly | Same one-state `chi` model | Original no-state theorem fails, but a second positive branch survives | Finite state-augmented collapse theorem with state data kept separate from Wilson coefficients |
| F-A8 | Proven | Local weight-spectrum finiteness dropped | Infinite low-weight STF tower | Candidate operator-space finiteness before reduction | No finite pre-reduction primitive-family catalog exists below the cutoff |

## Empirical Or Out-Of-Domain Non-Closure

| ID | Status | Scope break | Smallest explicit model | Exact package limit | Replacement path |
| --- | --- | --- | --- | --- | --- |
| F-A7 | Counterexample candidate | Fixed-order cutoff dropped | No explicit smallest all-orders model is promoted in the current repo | The fixed-order operator-finiteness statement no longer applies automatically | A separate all-orders closure theorem would be required |
| F-Scope | Counterexample candidate | Package applied outside the parity-even, nonspinning, local free-fall theorem domain | Nonmetric clock coupling or another out-of-domain sector | The closed theorem package no longer applies as stated | Separate sector-specific theory or empirical bookkeeping |

## Reading Rule

- Status: Proven. The negative uniqueness no-go is not listed here as a failure of the positive theorem, because it is a separate closed theorem branch.
- Status: Proven. Raw survivor-count, rank, and nullity notes are not failure statements unless an authoritative or supporting theorem note promotes them.
