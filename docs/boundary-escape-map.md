# Boundary Escape Map

- Status: Proven. This note records the exact assumption-drop failures of the closed free-fall theorem package.
- Status: Proven. It is not a work log; it is the paper-facing counterexample map.

| Boundary | Status | Smallest explicit counterexample | Exact theorem layer broken | Replacement bookkeeping data |
| --- | --- | --- | --- | --- |
| `A5`: analyticity dropped, locality kept | Proven | `m_A(Y)=m_0+\alpha e^{-1/Y^2}\Theta(Y)` | Lemma 55 fails at the finite analytic Taylor-jet step | Non-Taylor monopole germ data |
| `A3`: locality dropped, analyticity kept on the instantaneous variable side | Proven | One-coordinate causal power-law kernel | Local reduction to `m_A(Y(tau))` fails, so the local forms of Lemmas 55 and 56 no longer apply | Memory kernel or spectral data |
| `A4`: no-state assumption dropped, locality and analyticity kept | Proven | One-state local analytic `chi` model | Y-only reading of Lemmas 55 and 56 fails because the monopole response is not a function of `Y^I` alone | Finite local state-space data `(Y^I, chi^a)` and state-evolution parameters |
| `A8`: local weight-spectrum finiteness dropped | Proven | Infinite low-weight STF tower | Candidate scalar operator-space finiteness fails before reduction | No finite pre-reduction primitive-family catalog exists below the cutoff |

## A4 Salvage Branch

- Status: Proven. The `A4` counterexample does not kill finite collapse altogether.
- Status: Proven. Once finitely many local state variables are kept explicit, a finite state-augmented collapse theorem survives.
- Status: Proven. The augmented bookkeeping remains separate from higher-multipole Wilson coefficients.

## Reading Rule

- Status: Proven. Each row above keeps the primitive-family envelope and fixed-order counting distinct from the exact theorem layer that fails.
- Status: Proven. The counterexamples are assumption-drop sharp within the current repo scope.
- Status: Proven. None of these rows should be read as a reopening of the closed positive theorem inside its stated domain.
