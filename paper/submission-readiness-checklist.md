# Submission Readiness Checklist

This checklist closes the repository as a paper-ready dynamic-loophole
package. It does not claim that the prose is submission-final.

## Claim Consistency

- [x] Strongest honest claim is stated consistently: one-state is control/no-go;
  nonlinear second-order internal visibility is the live parameterized
  counterexample candidate.
- [x] The package states that the result is a parameterized design theorem,
  not an instrument forecast.
- [x] The package states that it is not a real-data detection claim.
- [x] The package states that arbitrary per-frequency projection nuisance
  collapses identifiability.
- [x] The package states that sideband existence alone is not a dynamic
  signature.
- [x] The package states that robustness-map fractions are not population
  measures.

## Artifact Traceability

- [x] Every major result is mapped to a repository artifact in
  [`claim-stack.md`](claim-stack.md).
- [x] The chronological result route is mapped in
  [`result-ledger.md`](result-ledger.md).
- [x] Non-claims are isolated in [`non-claims.md`](non-claims.md).
- [x] Reviewer risks are mapped to repository answers in
  [`reviewer-risk-notes.md`](reviewer-risk-notes.md).

## Draft Package

- [x] Clean submission-style draft exists:
  [`dynamic-loophole-submission-draft.md`](dynamic-loophole-submission-draft.md).
- [x] Existing `paper/manuscript.md` is not overwritten.
- [x] Existing status-labeled internal draft
  `paper/dynamic-loophole-manuscript.md` is not overwritten.
- [x] The clean draft removes internal `Status:` labels from prose.
- [x] The clean draft preserves the caveats in normal paper language.

## Figure And Table Package

- [x] Minimal four-figure package is specified in
  [`figure-inventory.md`](figure-inventory.md).
- [x] Minimal three-table package is specified in
  [`table-inventory.md`](table-inventory.md).
- [x] Every figure includes a one-sentence non-claim.
- [x] Every table includes a non-overclaim boundary.

## Scope Guard

- [x] No new loophole model is opened.
- [x] No new static family audit is opened.
- [x] No new empirical claim is added.
- [x] No instrument forecast is added.
- [x] No LLR, MLRS, PEP, Nutimo, runtime, build, clock-sector, or
  pulsar-timing branch is reopened.

## Required Verification

The required verification commands for this packaging pass are:

```powershell
python symbolic\analytic_phase_boundary.py
python symbolic\nonlinear_robustness_map.py
python symbolic\checks\test_symbolic.py
python -m compileall symbolic
git diff --check
```

- [x] Verification commands have been rerun after the final packaging edits.

## Human Polishing Queue

- [ ] Convert the clean draft into target journal style.
- [ ] Produce Figure 1 through Figure 4 from the mapped artifacts.
- [ ] Produce Table 1 through Table 3 from the mapped inventories.
- [ ] Decide whether a short system-anchored support note is needed before
  external circulation.
