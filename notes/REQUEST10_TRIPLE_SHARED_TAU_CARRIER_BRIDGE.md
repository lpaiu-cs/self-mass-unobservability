# Request 10.2b: Triple Shared-Tau Carrier Bridge

## Goal

Status: Counterexample candidate. Connect the Request 10.1 shared-`tau_chi`
multi-frequency theorem candidate to a concrete timing carrier inventory.

The binary tensor/tidal `3E` branch is now conditional after the exact-e +
Shapiro gate. This note therefore asks a cleaner question: can the existing
inner and outer GR clock carriers of a hierarchical triple provide enough
distinct frequencies to test the same dynamic-chi transfer law?

No triple timing runtime, source chasing, or empirical estimator is used here.

## Forcing Model

Status: Conjectural. Model the selected drive invariant for the same clock
body as

```text
F(t)
= F_in cos(Omega_in t + phi_in)
 + F_out cos(Omega_out t + phi_out).
```

Status: Proven. The one-state response remains linear:

```text
tau_chi dot chi + chi = alpha F(t).
```

The measured readout has transfer

```text
G(Omega) = c_Y + beta/(1 + i Omega tau_chi),
beta = alpha c_chi.
```

Thus the two triple carriers sample

```text
G(Omega_in),  G(Omega_out).
```

## Comparator Classes

Status: Proven. Each carrier separately remains degenerate with a
frequency-local `{F, dot F}` fit.

Status: Counterexample candidate. The pair
`Omega_in, Omega_out` is enough to break a real shared-coefficient degree-`1`
or degree-`2` derivative EFT, provided the two frequencies are distinct and
both carrier amplitudes are nonzero.

Status: Proven. The same two carriers are not enough to break all finite
derivative comparators. For real coefficients, degree `N >= 3` needs more
distinct positive frequencies. For the deliberately overpowered complex
polynomial comparator, degree `N >= 1` also needs more than two samples.

## Bridge Table

| carrier inventory | comparator | verdict | surviving target |
| --- | --- | --- | --- |
| one carrier | frequency-local `{F, dot F}` | degenerate | none |
| `Omega_in, Omega_out` | real degree-1 derivative EFT | distinguishable | shared-`tau_chi` odd-channel pole residual |
| `Omega_in, Omega_out` | real degree-2 derivative EFT | distinguishable | shared-`tau_chi` odd-channel pole residual |
| `Omega_in, Omega_out` | real degree-`N >= 3` derivative EFT | degenerate at two carriers | need more carrier frequencies |
| `Omega_in, Omega_out` | complex degree-0 polynomial | distinguishable | full rational transfer residual |
| `Omega_in, Omega_out` | complex degree-`N >= 1` polynomial | degenerate at two carriers | need at least three carriers |
| commensurate or missing carrier | any bridge comparator | collapse | one-frequency degeneracy |
| linear two-frequency drive | sideband basis | null | no sidebands |

## Interpretation

Status: Counterexample candidate. Existing triple inner/outer carriers can be
a meaningful dynamic-chi bridge without requiring a new tensor harmonic,
because the novelty target is the shared rational transfer relation, not a new
timing shape at each individual frequency.

Status: Counterexample candidate. The bridge is strongest against the physical
low-order real derivative comparator: two distinct carriers already test
degree `1` and `2`.

Status: Proven. The bridge is not yet a universal no-go breaker. If the
allowed comparator includes higher-order real derivative coefficients or a
complex degree-1 polynomial, two carriers are insufficient. The next input
would be either an additional surviving carrier frequency or a justified
truncation/power-counting prior on derivative order.

## Collapse Conditions

Status: Proven. The bridge collapses when

- `Omega_in = Omega_out`,
- `F_in = 0` or `F_out = 0`,
- `beta = alpha c_chi = 0`,
- `tau_chi = 0`,
- the whole frequency band lies below the derivative-truncation tolerance.

## Machine Outputs

Status: Proven. The bridge table is written to

```text
outputs/tsv/dynamic_chi_triple_shared_tau_bridge.tsv
outputs/json/dynamic_chi_triple_shared_tau_bridge.json
```
