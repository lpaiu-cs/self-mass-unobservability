# Request 4: Real-Data LLR Status

Request 4 is intentionally **not** the same thing as Request 3.

- Request 3 is a toy Sun-Earth-Moon sensitivity and injection-recovery study.
- Request 4 is a **real normal-point refit** on top of an existing LLR estimation stack.

That distinction is the whole reason the original Section 6 was repaired.

## What Is Already Done

The workspace already contains the pieces that belong **before** a real-data fit:

- Request 1: self-subtraction does not create a new COM force.
- Request 2: literal internal-structure self-unobservability is inconsistent.
- Request 3: mock LLR sensitivity and injection-recovery.
- Request 5 Phase A: strong-field prior translation for the tied model.
- Request 6: clock-sector toy joint fit.

So the EFT side is staged. What is still missing for Request 4 is the **production LLR analysis layer**.

## Role In The Project After Freezing Request 6

With Request 6 now effectively frozen as a support/local-audit branch, Request 4
becomes one of the two main routes to the tied-vs-decoupled verdict.

That project role should be read carefully:

- Request 4 is the weak-field / Solar-System side of the final verdict,
- Request 5 is the strong-field / pulsar free-fall side,
- Request 6 is now support material for the clock sector rather than the place
  where the final model choice is decided.

So Request 4 is no longer just "one more pending item." It is part of the
remaining mainline program.

## Why Request 4 Was Left Pending

The current repository has:

- no downloaded public LLR normal points,
- no validated normal-point ingest and quality-control chain,
- no established LLR force-model / partial-derivative estimator,
- no troposphere, station bias, reflector bias, and lunar-tidal estimation stack.

Without those, a "real-data posterior" would just be the toy integrator with public files glued on top, which is not a credible LLR result.

## External Inputs Required

The public data side is available, but it must be attached to a real estimator.

- The ILRS data page states that laser ranging data are distributed as full-rate data and normal points, and that LLR archives are maintained for Apollo and Lunokhod reflectors.
- The ILRS normal-point format page documents the data-record layout needed for ingest.
- The APOLLO normal-point page states that APOLLO normal points are public and also mirrored through the NASA CDDIS archive.

Those are the right starting points for Request 4, but they are still just data and formats. They are not the actual inference pipeline.

## What Counts As a Real Request 4 Implementation

A defensible Request 4 needs all of the following:

1. pick one concrete public release of normal points,
2. reproduce a GR baseline fit on that release,
3. add the EFT remapping

```math
\delta_{\rm SEP} = \sigma_1 (s_E-s_M) + \sigma_2 (s_E^2-s_M^2),
```

4. fit `sigma_1`, `sigma_2` together with at least:
   - Earth-Moon initial state,
   - station and reflector bias,
   - troposphere or Earth-orientation surrogate,
   - lunar tidal parameter,
   - SRP and thermal surrogates.

Only then is it meaningful to compare the recovered `sigma_1` scale with published Nordtvedt-style LLR bounds.

## Current Practical Status

So the honest state is:

- Request 4 is **specified**.
- Request 4 is **not numerically executed** in this workspace yet.
- The blocker is not the EFT. The blocker is the missing real LLR estimator and curated data ingest.

## Immediate Next Step

If Request 4 is the next priority, the next concrete action is:

1. choose the public data source to target first,
2. download a fixed release of those normal points into the workspace,
3. either attach to an existing LLR codebase or build a narrow ingest-plus-linearized-fit scaffold around that release.

Anything short of that is still Request 3 territory, not Request 4.

Given the present project state, that is also the correct strategic next move
if the goal is the tied-vs-decoupled verdict rather than further local clock
audits.
