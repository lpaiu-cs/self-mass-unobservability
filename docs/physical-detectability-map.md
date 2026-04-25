# Physical Detectability Map

Status: Counterexample candidate. This note replaces the dimensionless
relative cutoff `eta` by a parameterized SNR criterion. It does not import
real instrument numbers.

## Noise Model

Status: Conjectural. For each extracted harmonic, define

```math
{\rm SNR}_k=\frac{|\hat O_k|}{\sigma_k},
\qquad
\sigma_k=\sigma_0 k^q .
```

The parameters `sigma_0`, `q`, and `SNR_min` are channel inputs, not empirical
claims made by this repository.

Status: Counterexample candidate. A harmonic is physically usable when

```math
{\rm SNR}_k\ge {\rm SNR}_{\min}.
```

## Observable Amplitude

Status: Counterexample candidate. The modeled observable amplitude is

```math
|\hat O_k|
=
Q_0 |A_k(p,e)|
\left|
\frac{\rho^2}{\rho^2-k^2+i\delta k}
\right|
|\Lambda_k|,
```

where `Q_0` is a channel amplitude scale. For the acceleration-like channel
`|\Lambda_k|=1`; for the range-like template the audit uses

```math
|\Lambda_k|=|\kappa_r^2-k^2|^{-1},
```

with `kappa_r=kappa/n`.

## Detectability Criterion

Status: Proven. A physical design is not a novelty claim unless the detectable
set

```math
D=\{k:{\rm SNR}_k\ge{\rm SNR}_{\min}\}
```

satisfies

```math
|D|\ge N+2+K_\Lambda
```

and brackets the resonance:

```math
\exists k_-,k_+\in D,\qquad k_-<\rho<k_+ .
```

Status: Proven. If either condition fails, the correct verdict is an
SNR-underbudget no-go for the stated channel and noise model.

Status: Proven. Because the modeled SNR is linear in `Q_0`, a finite harmonic
audit also has a minimum physical scale:

```math
Q_{0,\min}
=
\min\left\{
Q_0:
|D(Q_0)|\ge N+2+K_\Lambda
\ \hbox{and}\ 
D(Q_0)\ \hbox{brackets}\ \rho
\right\}.
```

If this set is empty for the chosen harmonic cutoff, the design is a physical
detectability no-go for that cutoff.

## Default Map

Status: Counterexample candidate. The symbolic default map uses
`p=2`, `e=0.3`, `rho=3/2`, `delta=0.2`, `H=6`, and
`SNR_min=5`.

Status: Counterexample candidate. The default acceleration map varies
`Q_0` with `sigma_0=10^-3`. The default range map varies `Q_0` with
`sigma_0=10^-6` because the range projection template suppresses amplitudes
by `|\kappa_r^2-k^2|^{-1}` in the chosen dimensionless normalization.

Status: Counterexample candidate. These are proof-of-design examples, not
instrument forecasts.

Status: Counterexample candidate. With `N=1`, the default acceleration channel
has `K_Lambda=1` and requires `4` detectable complex samples. At
`sigma_0=10^-3`, it is SNR-underbudget for `Q_0=0.1` and `Q_0=0.3`, but
becomes physically budget-breaking at `Q_0=1.0`. The deterministic audit finds
`Q_{0,min}=0.909151` for `H=6`.

Status: Counterexample candidate. With `N=1`, the default range channel has
`K_Lambda=2` and requires `5` detectable complex samples. At
`sigma_0=10^-6`, it is SNR-underbudget for `Q_0=0.1`, but becomes physically
budget-breaking at `Q_0=0.3`. The deterministic audit finds
`Q_{0,min}=0.288909` for `H=6`.

## Verdict

Status: Counterexample candidate. The second-order line is now a conditional
observational design candidate: for a chosen channel noise model, one can
classify whether the eccentric harmonic set is physically budget-breaking,
SNR-underbudget, or unable to bracket the resonance.

Status: Proven. If the noise model leaves fewer detectable harmonics than
`N+2+K_Lambda`, the formal exact-in-e harmonic dictionary does not establish a
new observable.
