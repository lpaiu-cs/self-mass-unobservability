# Reproducibility Note

- Status: Note. This note records the release-stage validation path for the frozen theorem package.
- Status: Note. It separates reliable smoke checks from the slower monolithic symbolic runner.

## Prerequisites

- Status: Note. These checks require a Python environment with `numpy` and `sympy` installed (see the README Quick Start). Use the interpreter where `import sympy` succeeds: on this worktree that is `python3`, while inside an activated project `.venv` it may be `python`. The commands below use `python3`; substitute your environment's interpreter if it differs.

## Reliable Checks

- Status: Note. The following checks pass reliably in the current Windows worktree and are the recommended release validation path:
  - `python3 symbolic\checks\test_symbolic_smoke.py`
  - `python3 symbolic\family_envelope_census.py`
  - `python3 symbolic\irrep_family_census.py`
  - `python3 symbolic\weight_spectrum_demo.py`
  - `python3 symbolic\nonanalytic_jet_demo.py`
  - `python3 symbolic\hereditary_kernel_demo.py`
  - `python3 symbolic\stateful_counterexample_demo.py`
  - `python3 -m compileall symbolic`
- Status: Note. These checks cover the closed family-envelope statement, the A8 sharpening, the explicit A5/A3/A4 boundary escapes, and a syntax sweep of the symbolic package.

## Monolithic Runner Status

- Status: Note. The monolithic regression entry point `python3 symbolic\checks\test_symbolic.py` timed out in the current environment during the publication-freeze pass.
- Status: Note. That timeout is classified here as a reproducibility or infrastructure issue, not as an unresolved mathematical contradiction.
- Status: Note. The timeout does not by itself overturn any closed theorem statement, because the package-level smoke subset above still exercises the current authoritative closure statements directly.

## Recommended Smoke-Test Subset

- Status: Note. The current recommended reader-facing smoke path is:
  1. `python3 symbolic\checks\test_symbolic_smoke.py`
  2. `python3 symbolic\weight_spectrum_demo.py`
  3. `python3 symbolic\nonanalytic_jet_demo.py`
  4. `python3 symbolic\hereditary_kernel_demo.py`
  5. `python3 symbolic\stateful_counterexample_demo.py`
  6. `python3 -m compileall symbolic`
- Status: Note. This subset is intentionally lightweight and is the preferred reproducibility path for future readers unless the monolithic runner is separately hardened.

## Reading Rule

- Status: Note. A smoke-check pass means the frozen package claims and sharp boundary examples still materialize coherently in the current worktree.
- Status: Note. A monolithic-runner timeout should be treated as a tooling issue unless it is accompanied by a contradictory symbolic output or a failing package-level smoke check.
