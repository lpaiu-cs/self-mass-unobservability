# Lemma 65: Adiabatic Versus Dynamical State Boundary

## Statement

- Status: Proven. For the one-state local model

```math
\dot{\chi} = -\frac{1}{T_h}\chi + \frac{\beta}{T_h}Y,
```

the sharp `A4` boundary is not the mere existence of `\chi`.
- Status: Proven. The sharp boundary is whether `\chi` is algebraically or adiabatically eliminable back into a Y-only description, or instead carries independent orbital-timescale data.

## Adiabatic Or Slaved Regime

1. Status: Proven. Write the adiabatic parameter as

```math
\varepsilon_\chi := \Omega_{\rm orb} T_h.
```

2. Status: Proven. If `\varepsilon_\chi \ll 1` and the homogeneous mode is absent or decays, the formal solution gives

```math
\chi
=
\beta \left[
1 - T_h D_\tau + T_h^2 D_\tau^2 - \cdots
\right]Y.
```

3. Status: Proven. At leading order this collapses to the algebraic slaving relation `\chi=\beta Y`, with higher terms suppressed by orbital-time derivatives.
4. Status: Proven. Substituting back into `m_A(Y,\chi)` returns a Y-only local analytic monopole model, plus higher-derivative corrections that stay inside the local EFT.
5. Status: Proven. Therefore the adiabatic or slaved regime is not a sharp `A4` escape.

## Genuinely Dynamical Regime

1. Status: Proven. If `\varepsilon_\chi = O(1)` or the homogeneous mode survives as independent initial data, then `\chi(\tau)` is not an algebraic function of the instantaneous value `Y(\tau)`.
2. Status: Proven. In that regime two solutions can agree on `Y(\tau_0)` at a readout time `\tau_0` but disagree on `\chi(\tau_0)`.
3. Status: Proven. The monopole response `m_A(Y,\chi)` then distinguishes those states.
4. Status: Proven. Therefore the no-state theorem fails sharply in the genuinely dynamical regime.

## Verdict

- Status: Proven. Adiabatic eliminability is not enough to beat the original theorem.
- Status: Proven. A sharp `A4` escape requires a local state variable that survives on the orbital timescale or carries independent initial-state data.
