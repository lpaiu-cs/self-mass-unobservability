# Lemma 66: A4 Collapse Failure

## Statement

- Status: Proven. Admit the genuinely dynamical one-state local model of [`64-local-finite-state-construction.md`](64-local-finite-state-construction.md) and keep locality `A3`, analyticity `A5`, finite admitted primitive-family content `A8`, and the fixed-order local operator basis unchanged.
- Status: Proven. Then fixed family content, fixed-order candidate operator finiteness, and finite local normal-form closure remain intact.
- Status: Proven. The first broken step is the original no-state reduction to a monopole response of the form `m_A(Y)` depending only on the instantaneous scalar normal-form coordinates.
- Status: Proven. Therefore Lemma 55 fails in its original Y-only form, and Lemma 56 fails in its original Y-only sensitivity interpretation.

## Proof

1. Status: Proven. The dynamical one-state model introduces no new primitive family and no new local scalar operator class. It changes only the monopole bookkeeping sector by adding the explicit local state variable `\chi`.
2. Status: Proven. Because the local operator basis is unchanged, the finiteness statements of Lemmas 53 and 54 remain intact as local operator-space claims.
3. Status: Proven. Fix a readout time `\tau_0` and choose two local states with the same instantaneous `Y(\tau_0)` but different `\chi(\tau_0)`.
4. Status: Proven. Since

```math
m_A(Y,\chi)=m_0+\alpha Y+\lambda \chi,
```

one has

```math
m_A\!\left(Y(\tau_0),\chi_1(\tau_0)\right)
\neq
m_A\!\left(Y(\tau_0),\chi_2(\tau_0)\right)
```

whenever `\chi_1(\tau_0)\neq \chi_2(\tau_0)` and `\lambda \neq 0`.
5. Status: Proven. Therefore the monopole response is not a single-valued function of `Y(\tau_0)` alone.
6. Status: Proven. The exact failed step of the original theorem is the no-state prerequisite built into Lemma 55: the response is not of the form `m_A(Y)` before any Taylor expansion in `Y` can be attempted.
7. Status: Proven. Lemma 56 then fails secondarily in its original form because the monopole sensitivities are no longer determined by a Y-only jet.

## Boundary

- Status: Proven. This is a genuine `A4` failure, not an `A3` or `A5` failure, because locality and analyticity remain active in the explicit model.
- Status: Proven. The next step is not another failure statement but a salvage question: whether the collapse theorem revives once `\chi` is kept explicitly.
