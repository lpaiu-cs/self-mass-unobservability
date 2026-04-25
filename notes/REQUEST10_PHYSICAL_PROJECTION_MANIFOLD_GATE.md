# Request 10.5: Physical Projection-Manifold Derivation / Runtime-Worthiness Gate

## Goal

Status: Counterexample candidate. Decide whether the triple carrier projection
needed by the dynamic-chi bridge is closer to a finite shared physical
manifold or to arbitrary per-carrier complex nuisance.

This is the runtime-worthiness gate. No external timing runtime, source
chasing, PEP, Nutimo, or empirical estimator is used.

## Imported Boundary

Status: Imported from prior work. Request 10.4 classified projection nuisance
classes:

- calibrated/common real projection: `distinguishable`,
- finite shared real/complex projection: `conditional`,
- arbitrary per-carrier complex projection: `collapse`.

Status: Counterexample candidate. The remaining question is where the actual
hierarchical-triple timing projection sits in that class table.

## Carrier Coordinates

Status: Conjectural. Model the three measured carrier projection factors as
complex amplitudes:

```text
Lambda_in,
Lambda_out,
Lambda_c,
Omega_c = |Omega_in - Omega_out|.
```

The carrier vector has six real coordinates:

```text
Re Lambda_in, Im Lambda_in,
Re Lambda_out, Im Lambda_out,
Re Lambda_c, Im Lambda_c.
```

Status: Proven. If these six coordinates are independent nuisance parameters,
the projection can absorb the shared transfer relation point by point. This is
the `arbitrary per-carrier complex` collapse case.

## Minimal Physical Phase-Locked Manifold

Status: Counterexample candidate. The leading outer-dipole combination carrier
is not naturally phase-free. A minimal phase-locked projection has

```text
Lambda_in  = A_in  exp(i phi_in),
Lambda_out = A_out exp(i phi_out),
Lambda_c   = A_c   exp(i(phi_in - phi_out)).
```

Status: Proven. This manifold has five real parameters

```text
A_in, A_out, A_c, phi_in, phi_out
```

inside the six-real-dimensional carrier vector. Its generic Jacobian rank is
`5`, not `6`.

Status: Proven. The phase-lock constraint is

```text
arg Lambda_c - arg Lambda_in + arg Lambda_out = 0 mod pi.
```

Equivalently, the symbolic residual

```text
sin[(phi_in - phi_out) - phi_in + phi_out]
```

is exactly zero.

## Amplitude-Linked Variant

Status: Counterexample candidate. A stronger finite-manifold model links the
combination amplitude to the monopole amplitudes:

```text
Lambda_c = H A_in A_out exp(i(phi_in - phi_out + sigma)).
```

Status: Counterexample candidate. With calibrated `H` and floating `sigma`,
the generic Jacobian rank is still `5`; if both `H` and `sigma` are calibrated,
the rank drops further. This strengthens the statement that the physical
projection is finite and shared, not pointwise arbitrary.

## Gate Table

| projection manifold | generic rank | verdict | runtime worthiness |
| --- | ---: | --- | --- |
| calibrated/common real projection | 0 projection-rank contribution | distinguishable | runtime-motivated |
| phase-locked outer-dipole projection | 5 of 6 | conditional | runtime-motivated |
| amplitude-linked outer-dipole projection | 5 of 6 in the tested floating-orientation model | conditional | runtime-motivated |
| finite shared complex projection | at most finite template rank | conditional | conditional |
| independent complex carrier projection | 6 of 6 | collapse | not-runtime-motivated |
| independently floated combination phase | 6 of 6 | collapse | not-runtime-motivated |

## Interpretation

Status: Counterexample candidate. The physical projection-manifold audit moves
the dynamic-chi branch one step closer to runtime-worthiness. A phase-locked
outer-dipole carrier is a finite shared manifold, not arbitrary per-carrier
complex freedom.

Status: Proven. If a timing model floats the combination-carrier amplitude and
phase independently from the inner and outer carrier phases, the bridge
collapses exactly like arbitrary complex projection.

Status: Conjectural. The next runtime decision depends on whether the named
triple timing implementation enforces the phase/geometric links or effectively
fits independent complex amplitudes. Only the first case motivates external
runtime work as a scientific follow-up.

## Machine Outputs

Status: Proven. The physical projection-manifold gate is written to

```text
outputs/tsv/dynamic_chi_triple_projection_manifold_gate.tsv
outputs/json/dynamic_chi_triple_projection_manifold_gate.json
```
