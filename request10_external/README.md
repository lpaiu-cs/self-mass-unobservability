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

## Request 10.7b: joint-fit real-data upper limit

The conditional promotion of the 10.7a verdict was executed as the
pre-registered experiment `../notes/REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md`
(pre-registration committed before the first real-data fit). Additional
artifacts:

| file | content |
| --- | --- |
| `baseline_planetGR.npz` | real-data input of record (res/errs/toas of the 10.7a baseline run) |
| `joint_fit_upper_limit[_v2|_v2trunc].json` | pipeline outputs: v1 Jacobian / corrected v2 (QUOTED conservative) / truncated (CANDIDATE) |
| `beta_limit_curve[_v2|_v2trunc].tsv` | `tau_chi, beta_hat, sigma_beta, u95` over quoted [10,327] d + diagnostic grids |
| `finite_jacobian_v2.npy` (+`_meta.json`, `jac_v2/`) | corrected Jacobian: 7 planet columns re-derived at perturbative steps |
| `carrier_projection_rank_v2.json` | 10.7a Stage-2 gates re-adjudicated with v2 (unchanged: rank <= 5/6) |
| `jointfit_linearization_check.json`, `jointfit_lindiag.json` | 10.7b stop-rule firing + failure diagnosis (planet secant columns) |
| `jointfit_gateLA.json`, `gateLA_params.json` | 10.7c gates: A (absorption reality) PASS 0.003%, L (displacement linearity) FAIL 0.10 vs 0.05 (demoted to diagnostic by 10.7d) |
| `joint_fit_upper_limit_rn.json`, `beta_limit_curve_rn.tsv` | **QUOTED headline limit (10.7e)**: truncated + red-marginalized, window [1,327] d, min u95 = 0.202 us @ 1 d |
| `joint_fit_upper_limit_phys[_rn].json`, `beta_limit_curve_phys[_rn].tsv` | 10.7e E2 physical anchor: dimensionless beta_phys < 0.36-0.59 (Request-8 dictionary drive) |
| `pole_amplitude_bounds.json`, `residual_periodogram.json` | 10.7e E1 per-carrier pole bounds (inner/dif ~50 ns); low-frequency excess 12.6x justifying E3 |
| `joint_fit_upper_limit_10_7d.json`, `beta_limit_curve_10_7d.tsv` | 10.7d white-noise version (superseded by the pre-registered 10% rule at ratio 1.1001) |
| `jointfit_gateD.json`, `gateD_params.json`, `gateD_columns.npz` | 10.7d promotion gates: D1 (live null) 5/5 PASS, D2 (estimator calibration at u95) 5/5 PASS |
| `jointfit_dtheta52*.npy`, `xb52.npy`, `xc52.npy` | validation displacements and signal columns |
| `joint_fit_summary.md` | 10.7b/c/d verdict with Request 10.7 labels and the quoted limit |

Reproduction: `scripts/joint_fit_upper_limit.py` (deterministic, artifacts-only;
argv: jacobian, suffix, SV cut), `scripts/build_jac_v2.py`,
`scripts/rank_gate_v2.py`, and the WSL-side scripts described in
`../notes/REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md` and
`../notes/REQUEST10_EXTERNAL_JOINT_FIT_AMENDMENT_10_7C.md`.

## Request 10.8/10.8b: dynamic free-fall SEP channel (`sep_dynamic/`)

Design and gates: `../notes/REQUEST10_8_DYNAMIC_SEP_POLARIZATION_CHANNEL.md`,
`../notes/REQUEST10_8B_REALDATA_DYNAMIC_SEP_FIT.md`.

| file | content |
| --- | --- |
| `sep_dynamic/col_SEP_D.npz`, `sep_F0_report.json` | static SEP response column (F0; lever arm ~3.4e10 us/Delta, linear below the turn-wrap) |
| `sep_dynamic/col_R5_*.npz`, `sep_R5_report.json` | fixed-absorber hypothesis test (refuted; anchor downgraded to published 1.8e-6) |
| `sep_dynamic/patches/sepdyn.patch` | integrator patch: SEP_D(t) on pulsar pairs in rhs_GR_nbody (A=0 gate bit-exact) |
| `sep_dynamic/sep_ramp_probe.npz`, `sep_t2_results_v4.json` | time-convention calibration (integral-kernel 1/2 signature; TIMEDAYS units) |
| `sep_dynamic/sep_dynamic_columns.npz` | 6 exact dynamic response columns ({in,out,dif} x {cos,sin}) |
| `sep_dynamic/sep_feasibility_gates.json` | F2 (survival to 43.8%) and F3 (anchored 6.1e-8) — both PASS |
| `sep_dynamic/sep_joint_fit_10_8b.json`, `sep_beta_limit_curve_10_8b.tsv` | **QUOTED 10.8b limit: beta_ff < 1.1e-7 (95%, tau = 2 d), window [2,500] d** |
| `sep_dynamic/sep_gateG.json` | live gates: G1, G2a, G2b all PASS (6/6) |

Reproduction: `scripts/sep_static_sensitivity.py`,
`scripts/sep_feasibility_gates.py`, `scripts/sep_joint_fit_10_8b.py`
(all deterministic, artifacts-only), plus the WSL-side scripts per the notes.
