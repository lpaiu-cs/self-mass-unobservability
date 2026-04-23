# Request 5 Phase B: Build Feasibility

`REQUEST5_J0337_PHASEB_PUBLIC_INPUTS.md` established that public `TOA`, `par`,
and `Nutimo` code inputs for `PSR J0337+1715` are real and locally mirrored.

That leaves the next bounded question:

> can the current workspace actually build/run the public `Nutimo` release
> without turning `Request 5` into an infrastructure project?

This note closes that question for the current host.

## Artifacts

- feasibility script:
  [request5_j0337_phaseB_build_feasibility.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request5_j0337_phaseB_build_feasibility.py)
- machine summary:
  [request5_j0337_phaseB_build_feasibility_summary.json](/Users/lpaiu/vs/lab/self-mass-unobservability/request5_j0337_phaseB_build_feasibility_summary.json)
- visual summary:
  [request5_j0337_phaseB_build_feasibility_summary.svg](/Users/lpaiu/vs/lab/self-mass-unobservability/request5_j0337_phaseB_build_feasibility_summary.svg)
- diagnostic table:
  [request5_j0337_phaseB_build_feasibility_diagnostics.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request5_j0337_phaseB_build_feasibility_diagnostics.tsv)

## Host Result

The current host does **not** close the `Nutimo` runtime stack.

### First failure: default compiler frontend

`makefile-original` uses `g++`.
On this host that resolves to Apple clang, and the default build dies
immediately on:

```text
clang++: error: unsupported option '-fopenmp'
```

So the public release is not even using a locally compatible default compiler
setup out of the box.

### Second failure: missing scientific stack

When the compile probe is repeated with `CXX=g++-14`, the OpenMP frontend issue
goes away, but compilation then fails on:

```text
fatal error: boost/numeric/odeint.hpp: No such file or directory
```

That is consistent with the dependency declarations in the `Nutimo` README and
install script, which expect:

- `Boost`
- `Tempo2`
- `Minuit`
- `Acor`
- Cython/Python extension support

The current host probe finds:

- `tempo2`: missing
- `root-config`: missing
- `cython`: missing
- `mpic++`: missing
- `boost/numeric/odeint.hpp`: not found

## Interpretation

This is a materially different situation from `Request 4`.

For `Request 4`, the problem was that the public residual-bearing promotion
layer simply was not exposed.

For `Request 5`, the public path is genuinely open:

- public `TOA` files exist,
- public `Nutimo` code exists,
- public results bundles exist.

The blocker is therefore not public availability.
It is local dependency provisioning.

## Correct Boundary

Because this project has already used stop rules to avoid scope creep, the
right conclusion here is:

1. `J0337 Phase B` is a **real public-input branch**,
2. but it is **not yet a locally runnable branch** in the current workspace,
3. and the next escalation should **not** be an open-ended local dependency
   rebuild unless that infrastructure work becomes an explicit project goal,
4. so for the present project state, full `Phase B` execution should be treated
   as requiring a pre-provisioned research environment or an external-analysis
   route.

## Status Line

The most accurate one-line status is now:

> `J0337 Phase B` is publicly reproducible at the data/code level, but the
> current host does not close the `Nutimo` dependency/runtime stack, so the
> local workspace stops short of a full three-body TOA refit.
