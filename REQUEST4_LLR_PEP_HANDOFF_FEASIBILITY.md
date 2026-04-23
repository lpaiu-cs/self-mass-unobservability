# Request 4: PEP Handoff Feasibility

After the frozen `MLRS` branch, `PEP` was ranked as the strongest public-code
hand-off candidate for the final weak-field `LLR` posterior.

This note closes the next bounded question:

> can `PEP` actually be acquired and built in the current workspace without
> turning `Request 4` into a porting project?

## Artifacts

- feasibility script:
  [request4_llr_pep_handoff_feasibility.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_pep_handoff_feasibility.py)
- machine summary:
  [request4_llr_pep_handoff_summary.json](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_pep_handoff_summary.json)
- visual summary:
  [request4_llr_pep_handoff_summary.svg](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_pep_handoff_summary.svg)
- diagnostic table:
  [request4_llr_pep_handoff_diagnostics.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_pep_handoff_diagnostics.tsv)

## Acquisition Result

Public acquisition succeeded.

The bounded scout cloned:

- `pep_core` from `https://gitlab.com/jbattat/pep_core.git`
- `pep_doc` from `https://gitlab.com/jbattat/pep_doc.git`

Current pinned commits in this workspace:

- `pep_core`: `570410b`
- `pep_doc`: `146d799`

The repository contents are also consistent with a real `LLR`-capable code
family rather than a false positive:

- `README.md` says the package builds with `make`,
- `bigtest/` contains the `bigtest.obsllr` input,
- `peputil/makefile` explicitly uses `bigtest.obsllr` in the validation path.

So the external-estimator scout was directionally correct:
`PEP` is the right public code family to test next.

## Current Host Verdict

On the current machine, `PEP` is **not** locally buildable with the available
toolchain.

Observed host/toolchain:

- host: `Darwin arm64`
- compiler: `GNU Fortran (Homebrew GCC 14.2.0_1) 14.2.0`
- only detected `gfortran` path: `/opt/homebrew/bin/gfortran`

The blocker is not a missing-file problem.
It is a numeric-kind / toolchain-compatibility problem.

## Why It Fails

### 1. `REAL*10` is pervasive

The bounded scan finds:

- `166` files with `REAL*10`
- `207` source hits

This is not a one-file fix.

### 2. The current compiler rejects `REAL*10` directly

A trivial Fortran probe

```fortran
      real*10 x
      x = 1.0_10
      end
```

already fails with:

```text
Error: Old-style type declaration REAL*10 not supported
Error: Invalid real kind 10
```

So the first failure is a host/compiler capability mismatch, not a project
setup mistake.

### 3. Some sub-builds also expect `-m128bit-long-double`

Both

- `verify/makefile`
- `peputil/makefile`

hard-code `-m128bit-long-double`.

On the current host compiler this flag is rejected outright:

```text
gfortran: error: unrecognized command-line option '-m128bit-long-double'
```

So even beyond the first `REAL*10` failure, part of the package assumes a
different long-double environment than the current arm64 toolchain.

### 4. Minimal build repro fails immediately

The bounded reproducible failure is:

```text
make -C pep -f makefile a1ut1.o
```

This already dies at the first `REAL*10` source.

That is a better boundary than trying to push through the full root `make`,
because it proves the incompatibility without turning the scout into a long
compile log.

## What This Means For Request 4

This is **not** evidence that the `PEP` hand-off choice was wrong.

It means:

- `PEP` remains the strongest public-code estimator family found,
- but the current arm64 `gfortran` environment is not a viable local build
  target for it,
- and rewriting `REAL*10` usage locally would become a porting project.

That last point matters.
`Request 4` was explicitly kept bounded after the `MLRS` stop rule.
Turning `PEP` into a source-port effort would violate the same scope rule.

## Correct Next Step

The next bounded escalation is now very narrow:

1. try only a separate `x86_64` / legacy-compatible build environment for
   `PEP`,
2. do **not** start a local source-port of `PEP`,
3. if that bounded environment probe also fails or is unavailable, hand the
   final weak-field posterior to a collaboration / external-analysis-center
   route such as `JPL` or `POLAC/ELPN`.

## Status Line

The most accurate one-line status is now:

> `PEP` is confirmed as the right public external estimator family for
> `Request 4`, but it is not locally buildable on the current arm64 toolchain;
> the only justified next step is a bounded legacy/x86 environment probe, not a
> local porting effort.
