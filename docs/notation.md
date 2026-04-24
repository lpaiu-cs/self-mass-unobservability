# Notation

## Bodies And Orbit

- `A`: clock-carrying pulsar body.
- `B`: binary companion.
- `m_A`, `m_B`: body masses.
- `M = m_A + m_B`: total binary mass.
- `a`: relative semi-major axis.
- `e`: orbital eccentricity.
- `E`: eccentric anomaly.
- `r = a (1 - e cos E)`: Newtonian relative separation.
- `n`: mean motion.
- `c`: speed of light.
- `G`: Newtonian gravitational constant.

The MVP uses the Newtonian binary orbit as the carrier solution and expands
only the clock-rate map.

## Clock EFT

```text
d tau_A / dt
= 1 - v_A^2/(2 c^2)
  - alpha_A U_ext / c^2
```

with

```text
alpha_A = 1 + zeta_1 s_A + zeta_2 s_A^2 .
```

The decoupled clock correction is

```text
delta alpha_A = zeta_1 s_A + zeta_2 s_A^2 .
```

For the binary MVP,

```text
U_ext = G m_B / r .
```

The barycentric speed of body `A` is represented by

```text
v_A^2 = (m_B / M)^2 G M (2 / r - 1 / a)
```

at the Newtonian carrier-orbit level.

## Classification Labels

Each projected term is labeled with exactly one status:

- `observable`: produces a timing structure distinguishable from standard
  nuisance or post-Keplerian parameters at the stated order.
- `absorbed`: removable by redefining a standard nuisance parameter such as
  spin phase, spin frequency, or spin-frequency derivative.
- `degenerate`: has the same timing structure as a standard post-Keplerian
  parameter or another fitted model component.
- `null`: absent at the stated truncation or cancels before projection.

No row in the dictionary should remain unclassified.

## Harmonic Language

This repo records the harmonic in the clock-rate expansion and its timing
projection:

- constant clock-rate terms project to secular pulse phase or spin-frequency
  renormalization,
- `cos E` clock-rate terms project to the standard binary Einstein-delay
  harmonic after orbital-time integration,
- any harmonic not present at the declared truncation is marked `null`.
