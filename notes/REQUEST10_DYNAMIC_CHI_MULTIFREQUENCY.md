# Request 10.1: Shared-Tau Multi-Frequency Theorem/Counterexample

## Goal

Status: Counterexample candidate. Determine whether the one-state
dynamic-chi transfer remains distinguishable once the comparator is upgraded
from a single-frequency derivative fit to a shared-coefficient finite
derivative EFT across multiple drive frequencies.

No source chasing, runtime work, or empirical estimator work is used here.

## Model

Status: Proven. The model is

```text
tau_chi dot chi + chi = alpha F(t)
delta(m/m0) = c_Y F(t) + c_chi chi(t).
```

With `z = i Omega` and `beta = alpha c_chi`, the readout transfer is

```text
G(z) = c_Y + beta/(1 + tau_chi z).
```

Status: Proven. A finite derivative comparator of degree `N` has

```text
P_N(z) = sum_{n=0}^{N} d_n z^n.
```

The key rule is that the same coefficients `d_n` must fit every sampled
frequency. If the coefficients are allowed to change by frequency, the
problem collapses back to the one-frequency degeneracy.

## Comparator Classes

Status: Proven. Against an instantaneous static comparator, one finite
frequency already leaves a sine quadrature outside the basis.

Status: Proven. Against a one-frequency `{F, dot F}` comparator, the same
quadrature is degenerate. One derivative coefficient can fit the missing
quadrature at that frequency.

Status: Proven. Against a shared-coefficient complex polynomial comparator,
`N + 1` distinct samples can be fit exactly by interpolation.

Status: Counterexample candidate. Against that same generous complex
comparator, `N + 2` distinct off-pole samples expose the rational pole. The
interpolation residual after fitting `z_0, ..., z_N` and testing `z_*` is

```text
R_N(z_*) =
  beta (-tau_chi)^(N+1)
  prod_{j=0}^{N} (z_* - z_j)
  /
  prod_{j=0}^{N} (1 + tau_chi z_j)(1 + tau_chi z_*).
```

Status: Proven. With physical real derivative coefficients and positive
frequencies, the comparator splits into

```text
P_N(i Omega) = E_M(Omega^2) + i Omega O_L(Omega^2),
M = floor(N/2), L = floor((N-1)/2).
```

Status: Counterexample candidate. The real-coefficient comparator has a
sharper boundary. It can absorb at most `floor((N+1)/2)` positive frequency
samples before the next sample exposes the odd-channel pole

```text
Im G / Omega = - beta tau_chi/(1 + tau_chi^2 Omega^2).
```

## Minimality Table

| N_frequencies | comparator_basis | verdict | surviving_target |
| --- | --- | --- | --- |
| `1` | instantaneous static `F` | distinguishable from static-only | finite-`Omega tau_chi` quadrature |
| `1` | frequency-local `{F, dot F}` | degenerate | none |
| `K <= floor((N+1)/2)` | real shared-coefficient degree-`N` derivative EFT | degenerate | none at finite sampled set |
| `floor((N+1)/2)+1` | real shared-coefficient degree-`N` derivative EFT | distinguishable | shared-`tau_chi` odd-channel pole residual |
| `K <= N+1` | complex shared-coefficient degree-`N` polynomial | degenerate | none at finite sampled set |
| `N+2` | complex shared-coefficient degree-`N` polynomial | distinguishable | full rational transfer pole residual |
| any within `|Omega tau_chi| <= rho` | Taylor derivative EFT through order `N` | operationally degenerate below tolerance | none if tolerance exceeds `|beta| rho^(N+1)` |
| `2` | sideband search basis | null | none |

## Collapse Conditions

Status: Proven. The exact multi-frequency obstruction collapses when

- `beta = alpha c_chi = 0`,
- `tau_chi = 0`,
- sampled frequencies are not distinct,
- the overpowered complex sample set hits the pole `1 + tau_chi z_j = 0`.

Status: Proven. For physical `z_j = i Omega_j` with real positive
`Omega_j` and real positive `tau_chi`, the pole condition is never met.

Status: Proven. The low-frequency operational no-go boundary is

```text
|G(i Omega) - G_N(i Omega)| <= |beta| rho^(N+1)
```

for `|Omega tau_chi| <= rho`, after truncating the derivative expansion at
order `N`.

## Interpretation

Status: Counterexample candidate. The real novelty target is not harmonic
creation or a single phase-lag measurement. It is the shared finite-`tau_chi`
rational transfer relation across enough distinct drive frequencies.

Status: Counterexample candidate. For a physical real-coefficient degree-`N`
derivative comparator, the first exact multi-frequency distinguishability gate
is `floor((N+1)/2)+1` positive frequency samples.

Status: Counterexample candidate. If one grants the comparator complex
coefficients, the conservative gate is `N+2` distinct samples.

## Machine Outputs

Status: Proven. The machine-readable audit is written to

```text
outputs/tsv/dynamic_chi_multifrequency_audit.tsv
outputs/json/dynamic_chi_multifrequency_audit.json
```
