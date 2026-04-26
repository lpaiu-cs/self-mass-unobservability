# collapse-theorem

Status: Proven. This repository is now a closed theorem and counterexample package for the free-fall sector of the self-mass-unobservability program.

## Authoritative Entry Points

- Status: Proven. Start with [`docs/theorem-package.md`](docs/theorem-package.md) for the closed theorem domain, the main positive theorem, the negative uniqueness no-go, and the sharp boundary escapes.
- Status: Proven. Use [`docs/boundary-escape-map.md`](docs/boundary-escape-map.md) for the exact assumption-drop counterexamples and their replacement bookkeeping.
- Status: Proven. Use [`docs/paper-outline.md`](docs/paper-outline.md) for the paper-facing structure of the current repo.
- Status: Proven. Use [`docs/repo-status-matrix.md`](docs/repo-status-matrix.md) to distinguish authoritative notes, supporting theorem notes, historical notes, and raw bookkeeping.

## Verification

- Status: Proven. The symbolic regression entry point is `python symbolic\checks\test_symbolic.py`.
- Status: Proven. The A8 sharpness check is `python symbolic\weight_spectrum_demo.py`.
- Status: Proven. A lightweight syntax sweep is `python -m compileall symbolic`.
