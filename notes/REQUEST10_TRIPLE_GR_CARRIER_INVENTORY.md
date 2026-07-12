# Request 10.3: Triple GR Carrier Inventory / 3-Carrier Bridge

## Goal

Status: Counterexample candidate. Determine whether existing hierarchical
triple GR carrier structure can supply a third positive frequency sample for
the shared-`tau_chi` dynamic-chi transfer test.

The goal is not to find a new clock-sector harmonic. The dynamic-chi branch
uses carrier frequencies as samples of the same transfer law

```text
G(Omega) = c_Y + beta/(1 + i Omega tau_chi),
beta = alpha c_chi.
```

No external timing runtime, source chasing, or empirical estimator is used.

## Imported Triple Carrier Input

Status: Imported from prior work. The leading triple outer-potential
dictionary and basis-rank audit identified existing GR carrier structure:

- inner-only carrier,
- outer-only carrier,
- outer-gradient/dipole combination carrier,
- eccentric sidebands of the outer-dipole carrier.

Status: Proven. The scalar decoupled clock EFT did not create a new timing
shape from those carriers; it reduced to carrier-amplitude rescaling. That
does not remove the carrier frequencies themselves from the GR timing
inventory.

## Minimal Three-Carrier Inventory

Status: Counterexample candidate. The minimal inventory for the dynamic-chi
bridge is

```text
Omega_in,
Omega_out,
Omega_c = |Omega_in - Omega_out|.
```

Here `Omega_c` is the leading GR outer-dipole combination carrier. It is usable
as a transfer sample only if its projection amplitude is nonzero and not an
arbitrary complex nuisance independent of the other carriers.

## Genericity Conditions

Status: Proven. The three-carrier bridge collapses back to a two-carrier or
one-carrier bridge when

- `Omega_in = Omega_out`,
- `Omega_in = 2 Omega_out`,
- `Omega_out = 2 Omega_in`,
- the outer-dipole projection amplitude is zero,
- the projection nuisance is arbitrary per carrier.

The `2:1` exclusions ensure

```text
|Omega_in - Omega_out|
```

is distinct from both `Omega_in` and `Omega_out`.

## Comparator Pressure

Status: Counterexample candidate. With three distinct positive carriers, the
real shared-coefficient derivative comparator is pressured through degree
`N <= 4`.

This follows from the Request 10.1 boundary

```text
K_required = floor((N+1)/2) + 1.
```

For `K = 3`, this covers `N = 1, 2, 3, 4`.

Status: Proven. Three carriers are still insufficient for real degree `N >= 5`.

Status: Counterexample candidate. For the deliberately generous complex
polynomial comparator, three carriers are enough to distinguish degree
`N <= 1` because the conservative boundary is

```text
K_required = N + 2.
```

Status: Proven. Three carriers are insufficient for complex degree `N >= 2`.

## Inventory Table

| carrier set | comparator | verdict | surviving target |
| --- | --- | --- | --- |
| `Omega_in, Omega_out` | real derivative EFT | distinguishes degree 1 and 2 only | low-order pole residual |
| `Omega_in, Omega_out, |Omega_in-Omega_out|` | real degree `N <= 4` | distinguishable | shared-`tau_chi` odd-channel pole residual |
| `Omega_in, Omega_out, |Omega_in-Omega_out|` | real degree `N >= 5` | degenerate at three carriers | need more carrier frequencies or order prior |
| `Omega_in, Omega_out, |Omega_in-Omega_out|` | complex degree `N <= 1` | distinguishable | full rational transfer pole residual |
| `Omega_in, Omega_out, |Omega_in-Omega_out|` | complex degree `N >= 2` | degenerate at three carriers | need more carrier frequencies or order prior |
| eccentric outer-dipole sidebands | same gates | potential extension | higher-order comparator pressure |
| arbitrary projection per carrier | any finite carrier set | collapse | none |

## Interpretation

Status: Counterexample candidate. Request 10.3 strengthens the positive
dynamic-chi branch. Existing triple GR carriers can supply a third frequency
sample without relying on the binary tensor `3E` loophole.

Status: Counterexample candidate. The best current bridge is conditional:
if the outer-dipole combination carrier has a finite-dimensional projection
and the system is nonresonant, the three-carrier inventory reaches real
degree `N <= 4` and complex degree `N <= 1` comparator pressure.

Status: Proven. This is not a final generic observability theorem. Higher
real derivative order, higher complex polynomial order, or arbitrary
per-carrier projection nuisance still close the bridge unless more carrier
frequencies or a derivative-order prior are supplied.

## Machine Outputs

Status: Proven. The inventory table is written to

```text
outputs/tsv/dynamic_chi_triple_gr_carrier_inventory.tsv
outputs/json/dynamic_chi_triple_gr_carrier_inventory.json
```
