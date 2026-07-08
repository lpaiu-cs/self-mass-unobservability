# Reproducibility Note

- Status: Note. This note records the release-stage validation path for the frozen theorem package.
- Status: Proven. It separates reliable smoke checks from the slower monolithic symbolic runner.

## Reliable Checks

- Status: Proven. The following checks pass reliably in the current Windows worktree and are the recommended release validation path:
  - `python symbolic\checks\test_symbolic_smoke.py`
  - `python symbolic\family_envelope_census.py`
  - `python symbolic\irrep_family_census.py`
  - `python symbolic\weight_spectrum_demo.py`
  - `python symbolic\nonanalytic_jet_demo.py`
  - `python symbolic\hereditary_kernel_demo.py`
  - `python symbolic\stateful_counterexample_demo.py`
  - `python -m compileall symbolic`
- Status: Proven. These checks cover the closed family-envelope statement, the A8 sharpening, the explicit A5/A3/A4 boundary escapes, and a syntax sweep of the symbolic package.

## Monolithic Runner Status

- Status: Proven. The monolithic regression entry point `python symbolic\checks\test_symbolic.py` timed out in the current environment during the publication-freeze pass.
- Status: Proven. That timeout is classified here as a reproducibility or infrastructure issue, not as an unresolved mathematical contradiction.
- Status: Proven. The timeout does not by itself overturn any closed theorem statement, because the package-level smoke subset above still exercises the current authoritative closure statements directly.

## Recommended Smoke-Test Subset

- Status: Proven. The current recommended reader-facing smoke path is:
  1. `python symbolic\checks\test_symbolic_smoke.py`
  2. `python symbolic\weight_spectrum_demo.py`
  3. `python symbolic\nonanalytic_jet_demo.py`
  4. `python symbolic\hereditary_kernel_demo.py`
  5. `python symbolic\stateful_counterexample_demo.py`
  6. `python -m compileall symbolic`
- Status: Proven. This subset is intentionally lightweight and is the preferred reproducibility path for future readers unless the monolithic runner is separately hardened.

## Reading Rule

- Status: Proven. A smoke-check pass means the frozen package claims and sharp boundary examples still materialize coherently in the current worktree.
- Status: Proven. A monolithic-runner timeout should be treated as a tooling issue unless it is accompanied by a contradictory symbolic output or a failing package-level smoke check.
