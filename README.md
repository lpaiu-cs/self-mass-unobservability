# self-mass-unobservability

Status: Note. This repository carries two closed tracks: the frozen static theorem-and-counterexample package (Paper A), and the dynamic-chi measurement program that instruments the theorem's `A4` escape (Paper B).
Status: Note. The static mathematical content is frozen unless a direct contradiction is found during cleanup or later review.

## Authoritative Entry Points

- Start with [`docs/theorem-package.md`](docs/theorem-package.md) for the closed theorem domain, the main positive theorem, the negative uniqueness no-go, and the sharp boundary escapes.
- Use [`docs/boundary-escape-map.md`](docs/boundary-escape-map.md) for the exact assumption-drop counterexamples and their replacement bookkeeping.
- Use [`docs/paper-outline.md`](docs/paper-outline.md) for the paper-facing structure of the theorem track, and [`paper/paper-A-collapse-theorem.md`](paper/paper-A-collapse-theorem.md) for the current Paper A draft.
- Use [`docs/release-note.md`](docs/release-note.md) for the compact publication and handoff summary.
- Use [`docs/reproducibility-note.md`](docs/reproducibility-note.md) for the recommended smoke tests, and [`verification/`](verification/README.md) for the independent re-derivation of the family survivor counts (twice-corrected: the rank-4 undercount fix and the 2026-07-12 gradient-kinematics correction; corrected table `E,B,S,V,T,Q,U,Z = 5,16,30,15,17,23,17,21`).

## Dynamic-Chi Program (Requests 10.x, Paper B)

The dynamic track instruments the theorem's `F-A4+` state-augmented salvage — the finite state pair `(beta, tau_chi)` of an orbital-timescale internal state — on the pulsar triple PSR J0337+1715, and quotes the program's first real-data upper limits.

- Entry point: [`paper/paper-B-dynamic-sep-limit.md`](paper/paper-B-dynamic-sep-limit.md) (draft manuscript; LaTeX via `paper/build_paper_b.py`).
- Headline: `|delta Delta| < 1.68e-9` (95% statistical CL x `K_dyn = 10` systematic anchor; `tau_chi = 2 d`, worst drive phase) over `tau_chi in [2, 500] d`, per `request10_external/sep_dynamic/sep_phase_marg_10_8e.json`; Fisher-only floor `2.79e-10` and the full anchor bracket quoted alongside; clock-sector companion limits from the same data in the 10.7 chain.
- Request chain and pre-registrations: `notes/REQUEST10_*.md` (10.1 counting theorems through the 10.8f review response; every stage pre-registered and committed before its data look).
- Return artifacts, gate record, and deterministic reproduction scripts: [`request10_external/`](request10_external/README.md).
- The dynamic track's early working ledger is preserved at [`docs/failure-ledger-dynamic-chi.md`](docs/failure-ledger-dynamic-chi.md); the theorem track's boundary-risk register remains [`docs/failure-ledger.md`](docs/failure-ledger.md).

## Empirical Program (Requests 1-7)

The repository also carries the data-side program that motivated the theorem track. It is organized by request:

- `Request 1`: center-of-mass decoupling through quadrupole order.
- `Request 2`: internal-structure consistency check showing literal self-gravity removal is not a viable stellar-structure model.
- `Request 3`: mock LLR injection-recovery study for the weak-field free-fall sector.
- `Request 4`: real-data LLR ingest / estimator hand-off program. The data-side and software-interface scouts are present, but this repo does not yet contain a final weak-field posterior.
- `Request 5`: strong-field `J0337` program. Phase A is implemented; Phase B currently stops at public-input and build-feasibility gates.
- `Request 6`: local clock-sector audit with leave-one-out mass reconstruction and follow-up covariance probes.
- `Request 7`: provisional joint consistency scaffold combining the living strong-field and clock branches.

## Repository Layout

The repository is organized by request so the root stays clean and each stage
keeps its memo, script, notebooks, and generated outputs together.

- `request1/` ... `request7/`: per-request scripts, memos, notebooks, and checked-in outputs.
- `section6/`: repaired orchestration note for the staged program.
- [`data/`](data/README.md): mirrored public inputs, vendored external code, and bounded hand-off workspaces used by Requests 4 and 5.
- [`paper/`](paper/README.md): manuscript source and LaTeX draft for the interim paper.

## Quick Start

Use Python `3.10+` (tested in this workspace with Python `3.12`).

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Optional notebook environment:

```bash
pip install -r requirements-notebooks.txt
```

## Common Entry Points

You can run the main self-contained stages either directly or through `make`. If your default `python3` is older than `3.10`, override it explicitly, for example `make PYTHON=python3.12 request7`.

```bash
make request1
make request2
make request3
make request5-phaseA
make request6
make request7
make paper-tex
```

Equivalent direct commands:

```bash
(cd request1 && python3 request1_com_decoupling.py)
(cd request2 && python3 request2_internal_structure.py)
(cd request3 && python3 request3_llr_mock.py)
(cd request5 && python3 request5_j0337_phaseA.py)
(cd request6 && python3 request6_clock_sector.py)
(cd request7 && python3 request7_joint_consistency_scaffold.py)
```

Most scripts write outputs into their own request directory by default when run
from that directory or via `make`. Request 4 and Request 5 Phase B scouts
additionally depend on the mirrored assets under `data/`, and some of those
probes assume host tools such as `make`, `gfortran`, or macOS `textutil`.

## Request Index

- [`request1/REQUEST1_COM_DECOUPLING.md`](request1/REQUEST1_COM_DECOUPLING.md)
- [`request2/REQUEST2_INTERNAL_STRUCTURE.md`](request2/REQUEST2_INTERNAL_STRUCTURE.md)
- [`request3/REQUEST3_LLR_MOCK.md`](request3/REQUEST3_LLR_MOCK.md)
- [`request4/REQUEST4_LLR_REALDATA.md`](request4/REQUEST4_LLR_REALDATA.md)
- [`request5/REQUEST5_J0337_PHASEA.md`](request5/REQUEST5_J0337_PHASEA.md)
- [`request5/REQUEST5_J0337_PHASEB_PUBLIC_INPUTS.md`](request5/REQUEST5_J0337_PHASEB_PUBLIC_INPUTS.md)
- [`request5/REQUEST5_J0337_PHASEB_BUILD_FEASIBILITY.md`](request5/REQUEST5_J0337_PHASEB_BUILD_FEASIBILITY.md)
- [`request6/REQUEST6_CLOCK_SECTOR.md`](request6/REQUEST6_CLOCK_SECTOR.md)
- [`request6/REQUEST6_B1913_COVARIANCE.md`](request6/REQUEST6_B1913_COVARIANCE.md)
- [`request6/REQUEST6_LOW_SIDE_EXTENSIONS.md`](request6/REQUEST6_LOW_SIDE_EXTENSIONS.md)
- [`request6/REQUEST6_LOW_SIDE_COVARIANCE_PROXY.md`](request6/REQUEST6_LOW_SIDE_COVARIANCE_PROXY.md)
- [`request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md`](request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md)
- [`section6/SECTION6_REPAIRED.md`](section6/SECTION6_REPAIRED.md)

## Data And External Code

The workspace uses mirrored public releases and vendored external code snapshots
for the bounded Request 4 and Request 5 scouts, but the public repository keeps
only the manifests, memos, summaries, and patches needed to reconstruct those
steps. See [`data/README.md`](data/README.md) for the directory-level policy.

## Notes For Public Release

- Cache files and local machine artifacts are ignored in `.gitignore`.
- Generated artifacts and vendored data/code are marked in `.gitattributes` so GitHub language stats stay readable.
- `LICENSE` (MIT) and `CITATION.cff` are now included with release-stage citation metadata.
