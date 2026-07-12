# Request 10.7: External Nutimo Runtime-Worthiness Pilot

## Goal

Status: Counterexample candidate. Define the smallest external `Nutimo`
runtime pilot that can decide whether the dynamic-chi branch is still alive in
the named `J0337` timing implementation.

This request is not a local build task. It is a runtime-worthiness contract:
what must be checked, in what order, and what result stops the branch.

## Imported Boundary

Status: Imported from prior work. Request 10.6 source-inspected the public
`Nutimo` triple timing core and found a conditional positive implementation
class:

- standard integrated three-body core: finite state/geometry projection,
  `conditional`, `runtime-motivated`,
- explicit harmonic special case such as `RN_PL`: per-harmonic amplitude/phase
  nuisance, `collapse`, `not-runtime-motivated`.

Status: Imported from prior work. The target carrier set remains

```text
Omega_in, Omega_out, |Omega_in - Omega_out|.
```

Status: Imported from prior work. The target transfer law remains

```text
G(z) = c_Y + beta/(1 + tau_chi z).
```

## Pilot Stage 1: Configuration Closure

Status: Counterexample candidate. Before any Jacobian or injection test, the
runtime configuration must state:

- active parameter blocks,
- active delay flags,
- target carrier list,
- whether `specialcase` is `RN_PL` or an equivalent harmonic soak-up mode,
- whether any carrier-level amplitude/phase nuisance is enabled.

Status: Proven. If explicit per-harmonic amplitude/phase nuisance is enabled
on the target carriers, the pilot stops immediately. That configuration belongs
to the Request 10.6 collapse class.

Status: Counterexample candidate. If the standard integrated three-body core
is active and harmonic soak-up is disabled, proceed to the named Jacobian rank
gate.

## Pilot Stage 2: Named Jacobian Rank Gate

Status: Counterexample candidate. The first real runtime calculation should be
the finite fitted-parameter residual Jacobian near the baseline solution, not a
posterior run.

The gate is:

```text
rank(projection nuisance on 3 complex carriers) <= 5
and
dynamic_chi_column not in standard finite-parameter span.
```

Status: Counterexample candidate. If this holds, the named implementation is
consistent with the finite shared phase/geometric manifold from Request 10.5.
The dynamic-chi branch remains `conditional` and `runtime-motivated`.

Status: Proven. If the effective projection rank reaches `6/6`, or if the
dynamic-chi test column is absorbed by the named fitted-parameter span, the
branch collapses for this named implementation class.

## Pilot Stage 3: Minimal Synthetic Shared-Tau Injection

Status: Counterexample candidate. Only after the rank gate passes, inject a
minimal synthetic carrier perturbation obeying

```text
delta O_k proportional to Lambda_k F_k G(i Omega_k)
```

for the three target carriers.

Status: Counterexample candidate. The pass condition is not a final
astrophysical constraint. The pass condition is that after refitting the
standard finite nuisance parameters, the remaining carrier relation is still
best described by one shared `tau_chi`.

Status: Proven. If the standard finite parameters or an admitted harmonic
nuisance remove the shared-`tau_chi` relation, the pilot records a conditional
no-go rather than escalating to long runtime.

## Stop Rules

Status: Proven. Stop before runtime if the task becomes dependency repair
rather than a science test.

Status: Proven. Stop after configuration closure if harmonic carrier soak-up
is enabled on the target carriers.

Status: Proven. Stop after the Jacobian rank gate if the projection nuisance
is effectively rank `6/6` or if the dynamic-chi test column lies in the named
finite nuisance span.

Status: Counterexample candidate. Escalate to a real external runtime
experiment only if configuration closure and the rank gate both pass.

## Machine Outputs

Status: Proven. The pilot contract is written to

```text
outputs/tsv/dynamic_chi_nutimo_runtime_worthiness_pilot.tsv
outputs/json/dynamic_chi_nutimo_runtime_worthiness_pilot.json
```
