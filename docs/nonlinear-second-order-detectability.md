# Nonlinear Second-Order Detectability

Status: Counterexample candidate. This note promotes the M12 structural
result into a generated-line detectability audit. It asks whether exact-in-e
forcing supplies enough usable nonlinear generated samples to beat a finite
static nonlinear comparator after projection nuisance is counted.

## Question

Status: Proven. A single generated line does not count as novelty. A static
nonlinear local coefficient can match one generated complex sample.

Status: Counterexample candidate. The nonlinear second-order detectability
question is:

```math
\hbox{Does the usable generated-line set exceed the static nonlinear budget
while carrying the same }D_2(z)\hbox{ seen by the linear lines?}
```

The shared denominator is

```math
D_2(z)=\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2 .
```

For nonlinear internal drive,

```math
\mu_\chi\ddot\chi+\gamma_\chi\dot\chi+\omega_\chi^2(\chi-\alpha F)
=\beta_{F2}F^2,
```

the generated component at harmonic `k n` has the schematic magnitude

```math
|\hat O^{\rm gen}_k|
=
Q_\beta\, |B_k(p,e)|
\left|
\frac{1}{\rho^2-k^2+i\delta k}
\right|
|\Lambda_k|,
```

where

```math
B_k(p,e)=A_k(2p,e)
```

is the exact-in-e cosine coefficient of `F^2` when
`F=(a/r)^p`.

Status: Conjectural. This audit assumes the nonlinear generated component can
be separated from the linear component at the same harmonic by model fit,
calibration, or a joint multi-line diagnostic. If that component is not
separable, exact-in-e harmonics do not become clean generated observables.

## Usable Generated Samples

Status: Counterexample candidate. The relative-cutoff version defines

```math
W_k=
|B_k(p,e)|
\left|
\frac{1}{\rho^2-k^2+i\delta k}
\right|
|\Lambda_k|,
\qquad
U_\eta=\{k:W_k/W_{\max}\ge\eta\}.
```

Status: Counterexample candidate. The SNR version defines

```math
{\rm SNR}_k=
\frac{Q_\beta W_k}{\sigma_0 k^q},
\qquad
U_{\rm SNR}=\{k:{\rm SNR}_k\ge{\rm SNR}_{\min}\}.
```

The parameters `Q_beta`, `sigma_0`, `q`, and `SNR_min` are design inputs, not
empirical claims.

## Budget Criterion

Status: Proven. For a degree-`M` static nonlinear generated-line comparator
and `K_Lambda` shared projection nuisance parameters, the generated-line
sample count must satisfy

```math
|U|\ge M+2+K_\Lambda .
```

Status: Proven. For a resonance-assisted claim, the usable generated samples
must also bracket the internal denominator:

```math
\exists k_-,k_+\in U,
\qquad k_-<\rho<k_+ .
```

Status: Proven. If either condition fails, the result is a
generated-underbudget no-go for the stated cutoff or SNR model.

## Default Audit

Status: Counterexample candidate. The default audit uses `p=2`, `rho=3/2`,
`delta=0.2`, `H_gen=6`, degree `M=1`, and exact-in-e coefficients through
`B_k=A_k(4,e)`.

Status: Counterexample candidate. Under the relative cutoff `eta=10^-3`,
`e=0.1` is generated-underbudget for both acceleration and range channels.
For `e=0.3`, the acceleration channel has `5` usable generated samples against
the `4`-sample requirement, and the range channel has `6` usable generated
samples against the `5`-sample requirement. For `e=0.6`, both channels are
also generated-budget-breaking.

Status: Counterexample candidate. Under the SNR example with
`SNR_min=5`, `sigma_0=10^-3` for acceleration and `sigma_0=10^-6` for range,
the `e=0.3` acceleration channel is underbudget at `Q_beta=0.3` but
budget-breaking at `Q_beta=1.0`. The deterministic minimum scale is
`Q_beta=0.429654`.

Status: Counterexample candidate. In the same SNR example, the `e=0.3` range
channel is underbudget at `Q_beta=0.1` but budget-breaking at `Q_beta=0.3`.
The deterministic minimum scale is `Q_beta=0.118468`.

## Verdict

Status: Counterexample candidate. The nonlinear second-order branch has a
parameterized positive detectability region: for moderate eccentricity,
finite shared projection nuisance, and sufficient generated-line amplitude,
the exact-in-e nonlinear generated samples can exceed the static nonlinear
comparator budget while bracketing the shared quadratic denominator.

Status: Proven. This is still not an instrument forecast. It is a conditional
generated-line design theorem.

Status: Proven. If generated components cannot be separated, if projection
nuisance is arbitrary per frequency, if generated-line amplitudes fall below
the cutoff/SNR threshold, or if the sample count remains below
`M+2+K_Lambda`, the nonlinear second-order branch collapses back to a static
nonlinear mimicry/no-go statement.
