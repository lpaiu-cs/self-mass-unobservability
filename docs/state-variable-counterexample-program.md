# State-Variable Counterexample Program

## Target

- Status: Proven. The target of this boundary-stress program is no longer family-envelope closure, fixed-order local operator finiteness, nonanalytic monopole response, or genuinely hereditary memory.
- Status: Proven. The target is the smallest local analytic finite-state monopole model that keeps the finite primitive-family envelope and fixed-order local operator closure intact but breaks the original no-state theorem.
- Status: Proven. The same program must then ask whether a finite state-augmented collapse theorem survives once finitely many local state variables are kept explicitly.

## Layer Separation

- Status: Proven. The original local no-state theorem is the positive collapse theorem already closed under `A1`-`A8`, including the no-orbital-timescale-state assumption `A4`.
- Status: Proven. The finite-state local extension branch keeps locality `A3` and analyticity `A5` but admits finitely many explicit local state variables `\chi^a`.
- Status: Proven. The hereditary branch is separate: it drops locality `A3` and is now represented sharply by the causal power-law kernel model in [`hereditary-counterexample-program.md`](hereditary-counterexample-program.md).
- Status: Proven. The local nonanalytic branch is also separate: it keeps locality but drops analyticity `A5`, as recorded in [`nonanalytic-counterexample-program.md`](nonanalytic-counterexample-program.md).

## Canonical One-State Model

- Status: Proven. The smallest explicit local finite-state candidate uses one corrected scalar normal-form coordinate `Y(\tau)` and one local internal state `\chi(\tau)`:

```math
m_A(Y,\chi) = m_0 + \alpha Y + \lambda \chi,
\qquad
\dot{\chi} = -\frac{1}{T_h}\chi + \frac{\beta}{T_h}Y.
```

- Status: Proven. This model is local on the worldline and analytic in the variables kept explicitly.
- Status: Proven. It changes no primitive-family content and no local scalar operator catalog. It changes only the monopole bookkeeping layer.

## Boundary Split

- Status: Proven. The adiabatic or slaved control case is not a sharp `A4` escape, because when `\varepsilon_\chi := \Omega_{\rm orb} T_h \ll 1` and the homogeneous mode is absent or decays, `\chi` collapses back into an algebraic function of `Y`.
- Status: Proven. The genuinely dynamical one-state case is the sharp `A4` boundary, because when `\varepsilon_\chi = O(1)` or the homogeneous mode survives, the monopole response depends on independent orbital-timescale state data.

## Current Verdict

- Status: Proven. The smallest explicit local finite-state counterexample to the no-state theorem is the one-state local analytic `\chi` model above.
- Status: Proven. It breaks the original Y-only theorem because the monopole response is no longer a function of the instantaneous scalar normal-form coordinates `Y^I` alone.
- Status: Proven. A second positive branch survives: once finitely many local state variables are kept explicitly, a finite state-augmented collapse theorem survives as a finite-dimensional local state-space statement.
- Status: Proven. The augmented state data remain separate from higher-multipole Wilson coefficients.
