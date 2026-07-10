# Request 10.7a: External Nutimo Pilot — Return Artifacts

Status: Imported from prior work. This directory contains the **return
artifacts of the external collaborator contract** defined in
`notes/REQUEST10_EXTERNAL_NUTIMO_HANDOFF_PACKET.md` on branch
`lpaiu/minimal-nonlinear-sideband`. Per that packet, the Nutimo runtime run is
performed *outside* this repository; only the requested artifacts are returned
here.

Status: Proven. The runtime environment was provisioned and executed in WSL
Ubuntu-22.04 (outside the repo tree) on 2026-07-10, with explicit user
approval, from the public release of record:

- Code+data: Zenodo <https://zenodo.org/records/13899771> (Voisin et al.,
  A&A, doi:10.1051/0004-6361/202452100) — `nutimo.tar.bz2`, `Data.tar.bz2`,
  `Analysis.tar.bz2`.
- Baseline model: `parfile-planetGR-max-bestfit` (standard integrated
  three-body core + integrated 4th body / circum-ternary planet,
  `specialcase` empty ⇒ no RN_PL harmonic soak-up on the target carriers).

Note on scope: `AGENTS.md` (main branch) forbids reopening Nutimo/build work
*in this repo*; this delivery respects that boundary — the build and run live
in the external WSL environment, and this directory only carries the packet's
requested return artifacts on a side branch.

## Artifact inventory (as requested by the packet)

| file | packet item |
| --- | --- |
| `configuration_manifest.json` | configuration closure + target carrier inventory |
| `baseline_fit_summary.json` | preflight environment boundary + baseline fit |
| `finite_jacobian.npy` (+ `finite_jacobian_meta.json`) | ntoa×nfit residual derivative matrix |
| `carrier_projection_rank.json` | six-real-coordinate three-carrier projection rank + singular values |
| `dynamic_chi_test_column.tsv` (+ `dynamic_chi_span_test.json`) | synthetic shared-τ_χ residual columns + span membership |
| `synthetic_injection_recovery.json` | minimal shared-τ_χ injection, linearized refit recovery |
| `decision_summary.md` | one-page verdict with Request 10.7 labels |

Reproduction scripts (`jacobian_runner.py`, `prep_steps.py`,
`make_manifest.py`, `analysis_stage23.py`) are included under `scripts/`.
