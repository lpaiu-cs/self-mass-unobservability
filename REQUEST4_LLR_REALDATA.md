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

Two cautions should stay attached to that statement:

- the current APOLLO-only branch is a single-station scaffold, not yet a
  canonical weak-field verdict path,
- and the APOLLO public page itself says that the legacy text archive is an old
  format while the newer CRD stream is mirrored through CDDIS.

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
- A pinned public APOLLO normal-point release and ingest scaffold now exist in
  the workspace; see `REQUEST4_LLR_APOLLO_INGEST.md`.
- That ingest already exposed one archive inconsistency: the public `group_b`
  file currently contains `514` records while the APOLLO page describes `506`.
- A narrow baseline nuisance/design layer now also exists; see
  `REQUEST4_LLR_BASELINE_SCAFFOLD.md`.
- A first APOLLO-only baseline surrogate fit also now exists; see
  `REQUEST4_LLR_BASELINE_FIT.md`.
- That surrogate improves a geocentric `DE421` nominal residual by about
  `12x` in RMS, but still leaves residuals at `O(10^5-10^6 m)`.
- A CRD pivot scout now also exists; see `REQUEST4_LLR_CRD_PIVOT_SCOUT.md`.
- That scout shows the canonical CDDIS CRD root is real but currently gated by
  `Earthdata Login`, while the public ILRS support tarballs are directly
  accessible and can already be pinned and parsed in this workspace.
- So Request 4 is still **not a real weak-field parameter fit** in this
  workspace.
- The blocker is now very concrete: the self-built APOLLO-only surrogate does
  not close cleanly enough, which points toward an existing LLR estimator stack
  or the CRD/ILRS canonical path rather than indefinite growth of the bespoke
  branch. The pivot scout sharpens that further: authenticated CRD access
  and/or estimator integration is now the actual bottleneck.

## Immediate Next Step

If Request 4 remains the next priority, the next concrete action is now:

1. stop growing the APOLLO-only surrogate as if it were on track to become a
   final weak-field estimator by itself,
2. attach the already-pinned data branch to an existing LLR codebase or
   estimation stack, or else pivot to the CRD/ILRS canonical path,
3. only on top of that more credible estimator layer add the EFT remapping and
   nuisance estimation for `sigma_1`, `sigma_2`.

Anything short of that is still ingest territory, not the real Request 4 fit.

Given the present project state, that is also the correct strategic next move
if the goal is the tied-vs-decoupled verdict rather than further local clock
audits.
