# Request 8: Binary Clock Dictionary

## Purpose

Build the binary timing-observable dictionary for the decoupled clock sector of
the self-mass-unobservability EFT.

The working question is whether the decoupled clock EFT predicts anything in a
binary pulsar beyond Einstein-delay renormalization and trivial secular
spin-redshift renormalization.

## Assumptions

- free-fall sector = GR
- propagation sector = GR
- clock sector = decoupled EFT only
- binary first; no triple extension in this branch until the binary dictionary
  closes
- no PEP, MLRS, Tempo2, Nutimo, or other external runtime work
- Request 3 scale checks are not evidence for this dictionary

## Clock Model

```text
d tau_A / dt
= 1 - v_A^2/(2 c^2)
  - (1 + zeta_1 s_A + zeta_2 s_A^2) U_ext / c^2
```

For the binary MVP:

```text
U_ext = G m_B / r
r = a (1 - e cos E)
v_A^2 = (m_B / M)^2 G M (2 / r - 1 / a)
```

## Declared Expansion

- expansion variable: binary eccentricity `e` on a Newtonian Kepler ellipse
- post-Keplerian order: leading clock `O(c^-2)`
- default eccentricity order in the first machine table: `O(e^1)`
- periodic/secular split: constant rate terms are secular in pulse phase;
  `cos E` rate terms project to the usual Einstein-delay harmonic
- source status: no external source chasing; all terms originate from the
  clock EFT above

## First Skeleton Classification

The first table is intentionally conservative:

- GR kinetic and GR potential clock terms are baselines, not new clock-sector
  signals.
- `zeta_1 s_A U_ext` and `zeta_2 s_A^2 U_ext` inherit the binary potential
  harmonic at `O(c^-2)`.
- constant pieces are marked `absorbed` by spin phase/frequency structure.
- periodic `cos E` pieces are marked `degenerate` with Einstein-delay gamma at
  this truncation.
- absent independent harmonics are marked `null` rather than promoted to a
  conclusion.

This is not yet a final gamma-only theorem.
