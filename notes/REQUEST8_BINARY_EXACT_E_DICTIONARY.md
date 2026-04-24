# Request 8: Binary Exact-in-e Dictionary

## Scope

This memo upgrades the binary clock dictionary from an `O(e^1)` rate-harmonic
projection to an exact-in-e leading-order timing projector.

Assumptions remain:

- free-fall sector = GR,
- propagation sector = GR,
- clock sector = decoupled EFT only,
- carrier orbit = Newtonian Kepler binary,
- post-Keplerian order = leading clock `O(c^-2)`,
- no triple systems, source chasing, or external timing runtimes.

## Eccentric-Anomaly Projector

The Newtonian Kepler relation is

```text
M = E - e sin E = n (t - T0),
```

so

```text
dt/dE = (1 - e cos E) / n .
```

Each clock contribution is projected as

```text
d tau_i / dE = (d tau_i / dt) (dt/dE),
Delta_i(E) = integral (d tau_i / dE) dE .
```

The spin-frequency absorption is performed against coordinate time,

```text
t - T0 = (E - e sin E) / n,
```

not against `E` alone. This distinction is why a term that is constant in
`d tau/dE` can still leave an `e sin E` Einstein-delay-basis residual after the
secular coordinate-time trend is removed.

## Exact Contributions

With

```text
r = a (1 - e cos E)
U_ext = G m_B / r
v_A^2 = (m_B/M)^2 G M (2/r - 1/a),
```

the exact leading-order pieces are:

```text
baseline kinetic:
d tau/dE = -G m_B^2 (1 + e cos E) / (2 M a c^2 n)

baseline potential:
d tau/dE = -G m_B / (a c^2 n)

zeta_1 s_A clock:
d tau/dE = -G m_B s_A zeta_1 / (a c^2 n)

zeta_2 s_A^2 clock:
d tau/dE = -G m_B s_A^2 zeta_2 / (a c^2 n)
```

The integrated residuals contain only a coordinate-time secular piece and a
single `sin E` periodic piece.

## Einstein-Delay Basis Comparison

After removing the fitted secular spin-redshift trend, every nonzero periodic
piece is proportional to

```text
sin E .
```

The exact-in-e periodic coefficients are:

```text
baseline kinetic:   -G e m_B^2 / (M a c^2 n)
baseline potential: -G e m_B / (a c^2 n)
zeta_1 term:        -G e m_B s_A zeta_1 / (a c^2 n)
zeta_2 term:        -G e m_B s_A^2 zeta_2 / (a c^2 n)
```

Thus the decoupled clock terms do survive symbolically, but their periodic
timing shape is exactly the Einstein-delay `sin E` basis at this leading order.

## Classification

- Secular pieces: `absorbed` by spin phase or spin frequency.
- Periodic `sin E` pieces: `degenerate` with Einstein delay gamma.
- Independent periodic structure beyond the Einstein-delay basis: none in this
  exact-in-e leading-order binary projector.

## Theorem Status

The exact-in-e result strengthens the gamma-only theorem candidate for the
leading-order binary clock sector: no distinguishable timing basis beyond
secular spin-redshift structure and Einstein-delay `sin E` appears here.

This memo does not claim the final repository theorem. The binary branch still
needs the planned theorem-or-counterexample closure note before the stop rule is
formally invoked.
