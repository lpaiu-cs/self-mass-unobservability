# Request 4: PEP Environment Probe

`REQUEST4_LLR_PEP_HANDOFF_FEASIBILITY.md` already established that:

- `PEP` is the right public external estimator family,
- but the current `arm64` `gfortran` cannot build it because of pervasive
  `REAL*10` usage and `-m128bit-long-double` assumptions.

That left exactly one justified local escalation:

> is there already a bounded `x86_64` / legacy-compatible environment on this
> machine that would let `PEP` build without turning into a provisioning or
> porting project?

This note closes that question.

## Artifacts

- probe script:
  [request4_llr_pep_environment_probe.py](request4_llr_pep_environment_probe.py)
- machine summary:
  [request4_llr_pep_environment_summary.json](request4_llr_pep_environment_summary.json)
- visual summary:
  [request4_llr_pep_environment_summary.svg](request4_llr_pep_environment_summary.svg)
- diagnostic table:
  [request4_llr_pep_environment_diagnostics.tsv](request4_llr_pep_environment_diagnostics.tsv)

## Result

The bounded local environment probe is negative.

What the host has:

- `Rosetta` is available

What the host does **not** already have:

- no `/usr/local` `Homebrew` prefix
- no `/usr/local/bin/gfortran*`
- no `/usr/local/bin/gcc*`
- no detected `docker`, `podman`, `colima`, `lima`, `nerdctl`, or `orb`
  runtime

So there is no already-provisioned local path to:

- launch an `x86_64` toolchain,
- run a containerized legacy build,
- or perform a bounded off-architecture `PEP` compile.

## Consequence

This matters because `Request 4` was already under an explicit bounded-project
rule.

At this point, installing a new `x86_64` toolchain stack, reviving an old build
environment, or standing up a new VM/container runtime would no longer be a
short hand-off probe.
It would be new infrastructure work.

That means the correct project consequence is:

1. stop local `PEP` escalation here,
2. keep the `PEP` scout artifacts as proof that the external estimator family
   is correct,
3. treat actual `PEP` execution as requiring a separate pre-provisioned
   environment or outside analyst / collaboration support,
4. keep `Request 4`'s final weak-field posterior on the external-estimator side,
   not in this local workspace.

## Status Line

The most accurate one-line status is now:

> `PEP` is still the right external estimator family for `Request 4`, but this
> host does not already provide the bounded legacy/x86 environment needed to run
> it, so local escalation stops here and the final weak-field posterior must be
> handed off to an external estimator environment or analysis route.
