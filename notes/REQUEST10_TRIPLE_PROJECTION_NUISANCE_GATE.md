# Request 10.4: Triple Projection-Nuisance Realism Gate

## Goal

Status: Counterexample candidate. Determine whether the three-carrier
dynamic-chi bridge

```text
Omega_in, Omega_out, |Omega_in-Omega_out|
```

survives once the triple timing projection from carrier response to measured
observable is allowed to carry nuisance parameters.

No timing runtime, source chasing, or empirical estimator is used here. The
task is a symbolic observability gate.

## Imported Result

Status: Imported from prior work. Request 10.3 established that a nonresonant
hierarchical triple can provide three existing GR carrier samples:

```text
Omega_in,
Omega_out,
Omega_c = |Omega_in-Omega_out|.
```

Status: Counterexample candidate. With finite projection nuisance, those three
samples pressure real shared-coefficient derivative comparators through degree
`N <= 4` and the deliberately generous complex comparator through degree
`N <= 1`.

Status: Proven. The same bridge collapses if the projection nuisance is
arbitrary and complex at every carrier.

## Projection Model

Status: Conjectural. Write the measured carrier amplitude as

```text
O_k = Lambda_k(theta) G(z_k) F_k,
z_k = i Omega_k,
G(z) = c_Y + beta/(1 + tau_chi z),
beta = alpha c_chi.
```

The core question is not just whether there are three carriers. It is whether
the three projection factors `Lambda_k(theta)` are tied by a finite shared
parameter vector `theta`, or whether they behave as independent complex
nuisance amplitudes.

## Nuisance Classes

Status: Counterexample candidate. A calibrated common real scale preserves the
bridge:

```text
Lambda_k = Gamma.
```

The scale either deprojects directly or cancels in transfer ratios. It cannot
absorb frequency-dependent phase and amplitude changes carrier by carrier.

Status: Counterexample candidate. A finite real shared projection nuisance,
for example

```text
Lambda_k = g0 + g1 q_k
```

or the range-like form

```text
Lambda_k = Gamma/(kappa^2 + z_k^2),
```

does not collapse the bridge pointwise. It does, however, consume rank and
requires a projection-prior or Jacobian audit before being promoted to an
observable claim.

Status: Counterexample candidate. A finite complex shared projection

```text
Lambda_k = ell0 + ell1 z_k
```

is an overpowered but still finite comparator. It is not equivalent to
independent complex freedom at every carrier, but three carriers are not enough
to make this class a final positive result without a prior, more carriers, or
an explicit projection model.

Status: Proven. An arbitrary per-carrier complex projection collapses the
bridge:

```text
Lambda_k = O_k/(G(z_k) F_k).
```

This choice absorbs the shared transfer law point by point. Under this
comparator, the three-carrier bridge is not runtime-motivated.

## Gate Verdict

| projection class | verdict | runtime relevance | reason |
| --- | --- | --- | --- |
| calibrated common real scale | distinguishable | runtime-motivated | deprojection preserves `G(z_k)` |
| unknown common real scale | distinguishable | runtime-motivated | ratios remove the scale |
| finite real shared geometry | conditional | conditional | finite, not pointwise arbitrary, but rank/prior dependent |
| range-like shared projection | conditional | conditional | survives if `Gamma,kappa` are shared or calibrated |
| finite complex shared projection | conditional | conditional | shared but high-dimensional enough to require more carrier/prior input |
| arbitrary per-carrier complex projection | collapse | not-runtime-motivated | absorbs `G(z_k)` independently |
| zero or singular projection | collapse | not-runtime-motivated | carrier sample is lost |

## Interpretation

Status: Counterexample candidate. Request 10.4 keeps the dynamic-chi branch
alive under calibrated or finite shared projection assumptions. This is a
projection-realism gate, not a measured-observable claim.

Status: Proven. If the realistic triple timing projection is effectively an
arbitrary complex amplitude for each carrier, the shared-`tau_chi` bridge
collapses.

Status: Conjectural. The next decisive issue is the physical size of the
projection nuisance manifold. External runtime becomes scientifically
motivated only if the relevant triple timing projection is finite-dimensional
or independently calibrated enough that it cannot emulate independent complex
amplitudes at the three carriers.

## Machine Outputs

Status: Proven. The projection-nuisance gate table is written to

```text
outputs/tsv/dynamic_chi_triple_projection_nuisance_gate.tsv
outputs/json/dynamic_chi_triple_projection_nuisance_gate.json
```
