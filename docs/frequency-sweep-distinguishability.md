# Frequency-Sweep Distinguishability

Status: Counterexample candidate. M2 asks whether the relaxation transfer

```text
G(z) = c_Y + beta/(1 + tau_chi z),  beta = alpha c_chi,  z = i Omega
```

can be absorbed by a finite local derivative EFT with shared coefficients,

```text
P_N(z) = sum_{n=0}^{N} d_n z^n.
```

## Generous Complex-Coefficient Boundary

Status: Proven. If the coefficients `d_n` are allowed to be complex and chosen
freely for one finite frequency set, then `N + 1` distinct complex sample
points can be fit exactly by polynomial interpolation.

Status: Proven. Therefore, `K <= N + 1` sampled frequencies do not distinguish
the one-state relaxation model from this deliberately overpowered degree-`N`
polynomial comparator.

Status: Proven. If the same degree-`N` polynomial is required to match the
relaxation transfer at `N + 2` distinct sample points `z_0, ..., z_N, z_*`,
then the residual after fitting the first `N + 1` points is

```text
R_N(z_*) =
  beta (-tau_chi)^(N+1)
  prod_{j=0}^{N} (z_* - z_j)
  /
  prod_{j=0}^{N} (1 + tau_chi z_j)(1 + tau_chi z_*).
```

Status: Proven. For distinct finite frequencies, real positive `tau_chi`, and
`beta != 0`, this residual is nonzero.

Status: Proven. The obstruction remains valid even if the polynomial
coefficients are allowed to be complex, so it also applies to the physically
smaller real-coefficient derivative EFT class.

## Real-Coefficient Derivative EFT

Status: Proven. If the derivative coefficients are real and the sweep uses
positive real frequencies, then

```text
P_N(i Omega) = E_M(Omega^2) + i Omega O_L(Omega^2),
M = floor(N/2),  L = floor((N-1)/2).
```

Status: Proven. The relaxation transfer has

```text
Re G = c_Y + beta/(1 + tau_chi^2 Omega^2),
Im G / Omega = - beta tau_chi/(1 + tau_chi^2 Omega^2).
```

Status: Proven. The odd channel gives the sharper positive-frequency boundary:
a degree-`N` real-coefficient derivative EFT can absorb at most
`floor((N+1)/2)` positive frequency samples by interpolation before the next
sample exposes the pole.

Status: Proven. The residual in the odd channel after fitting the allowed
samples is the interpolation residual of
`-beta tau_chi/(1 + tau_chi^2 u)` in `u=Omega^2`, so it is nonzero for
`beta tau_chi != 0` and distinct positive frequencies.

## Theorem Candidate

Status: Counterexample candidate. With physical real derivative coefficients,
a frequency sweep over at least `floor((N+1)/2)+1` distinct positive drive
frequencies distinguishes the one-state relaxation response from any
shared-coefficient degree-`N` local derivative EFT, provided

```text
beta != 0,  tau_chi != 0,
Omega_j != Omega_k for j != k.
```

Status: Counterexample candidate. Even if one grants the comparator complex
coefficients, `N + 2` distinct drive frequencies distinguish the same response,
provided

```text
beta != 0,  tau_chi != 0,
z_j != z_k for j != k,
1 + tau_chi z_j != 0.
```

Status: Proven. For the physical sweep `z_j = i Omega_j` with real
`Omega_j > 0` and `tau_chi > 0`, the pole condition
`1 + tau_chi z_j = 0` is never met.

## Low-Frequency No-Go Boundary

Status: Proven. The Taylor derivative approximation through order `N` is

```text
G_N(z) = c_Y + beta sum_{n=0}^{N} (-tau_chi z)^n.
```

Status: Proven. Its exact residual is

```text
G(z) - G_N(z) =
  beta (-tau_chi z)^(N+1) / (1 + tau_chi z).
```

Status: Proven. For `z=i Omega` and `|Omega tau_chi| <= rho`, the residual
magnitude obeys

```text
|G(i Omega) - G_N(i Omega)| <= |beta| rho^(N+1).
```

Status: Proven. Thus the relaxation model collapses operationally to a
degree-`N` derivative EFT over a low-frequency band whenever the experimental
tolerance is larger than this truncation bound.

## M2 Classification

Status: Proven. One-frequency quadrature is not decisive against a derivative
EFT comparator.

Status: Proven. `N + 1` frequency points are not decisive against a degree-`N`
complex polynomial comparator because interpolation can absorb them.

Status: Counterexample candidate. The first exact physical distinguishability
target is `floor((N+1)/2)+1` distinct positive frequencies with shared real
derivative coefficients.

Status: Counterexample candidate. The conservative complex-coefficient target
is `N + 2` distinct frequencies.

Status: Counterexample candidate. The observable class is the pole-bearing
frequency-transfer relation, not harmonic creation.

## Link To M3 Dictionary

Status: Counterexample candidate. The first concrete source of positive
frequency samples is the orbital-harmonic template in
[`forcing-observable-dictionary.md`](forcing-observable-dictionary.md).

Status: Proven. For a power-law invariant `(a/r)^p`, the `n` and `2n`
harmonics through `O(e^2)` provide the first two-point positive-frequency test
whenever `e != 0`, `p != 0`, and `p != -3`.

Status: Counterexample candidate. Projection-channel assumptions are audited in
[`projection-channel-audit.md`](projection-channel-audit.md). The frequency
sweep theorem applies directly after deprojection when `Lambda(z)` is known, or
after nuisance accounting when `Lambda(z)` has shared finite parameters.

## Machine Outputs

Status: Proven. The frequency-sweep audit writes the minimality table to

```text
outputs/tsv/dynamic_chi_multifrequency_audit.tsv
outputs/json/dynamic_chi_multifrequency_audit.json
```
