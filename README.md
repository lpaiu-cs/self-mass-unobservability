# self-mass-unobservability

Research workspace for staged calculations around self-mass unobservability and the tied-vs-decoupled EFT question. The repository combines analytic checks, mock-data studies, public-input scouts, and provisional joint consistency scaffolds.

## Current Scope

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
- This pass does not choose a publication license or citation metadata. Add `LICENSE` and `CITATION.cff` before the final public release if you want GitHub to expose them explicitly.
