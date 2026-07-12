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

## Request 10.8-10.8f: dynamic free-fall SEP channel (`sep_dynamic/`)

Design and gates: `../notes/REQUEST10_8_DYNAMIC_SEP_POLARIZATION_CHANNEL.md`,
`../notes/REQUEST10_8B_REALDATA_DYNAMIC_SEP_FIT.md`, then
`..._10_8C/D/E` and the 10.8f review response
(`../notes/REQUEST10_8F_REVIEW_RESPONSE.md`: an independent review caught a
causal/anticausal quadrature swap; every data-touching stage was re-run on
the corrected causal templates — the rows below are the corrected record).

| file | content |
| --- | --- |
| `sep_dynamic/col_SEP_D.npz`, `sep_F0_report.json` | static SEP response column (F0; lever arm ~3.4e10 us/Delta, linear below the turn-wrap) |
| `sep_dynamic/col_R5_*.npz`, `sep_R5_report.json` | fixed-absorber hypothesis test (refuted; anchor downgraded to published 1.8e-6) |
| `sep_dynamic/patches/sepdyn.patch` | integrator patch: SEP_D(t) on pulsar pairs in rhs_GR_nbody (A=0 gate bit-exact) |
| `sep_dynamic/sep_ramp_probe.npz`, `sep_t2_results_v4.json` | time-convention calibration (integral-kernel 1/2 signature; TIMEDAYS units) |
| `sep_dynamic/sep_dynamic_columns.npz` | 6 exact dynamic response columns ({in,out,dif} x {cos,sin}) |
| `sep_dynamic/sep_feasibility_gates.json` | F2 (per-column survival 8-27%; causal shared-tau peak 21.1%) and F3 (anchored 6.87e-8) — both PASS |
| `sep_dynamic/sep_joint_fit_10_8b.json`, `sep_beta_limit_curve_10_8b.tsv` | zero-phase fit: no detection (Z = 1.328, p = 0.39; anticausal control quiet); K10 curve min 1.36e-9 (2 d) |
| `sep_dynamic/anchor_resolution_10_8c.json`, `sep_beta_limit_curve_10_8c.tsv` | K anchor diagnosis (rho spectrum; K_dyn = 10, rule branch 1) + re-anchored curve; generator `scripts/anchor_resolution_10_8c.py` |
| `sep_dynamic/turnV_columns.npz`, `turn_search_10_8d.json` | turn-fixing measurement + lattice turn-search: PASS on the tested lattice (7 viable cells; worst beta dev 1.01e-10 vs tol 2.07e-10) |
| `sep_dynamic/sep_phase_marg_10_8e.json`, `sep_limit_curve_10_8e.tsv` | **HEADLINE: \|dDelta\| < 1.68e-9 (95% statistical CL x K_dyn = 10, worst phase, tau = 2 d)**; Fisher floor 2.79e-10; full-rank bracket 3.53e-9; K934 1.57e-7; no detection (E1 p = 0.26, control equally elevated) |
| `sep_dynamic/sep_gateG2.json`, `sep_gateG2wp.json`, `gateG2_inputs.npz`, `gateG2_params.json` | **live gates of record (C4-complete harness): 9/9 PASS at tau = 2/18/52 d + 3/3 PASS on the worst-phase headline pair (t_off = 104.35 d)** — G1 offsets <= 0.039 sigma_F, G2a <= 0.096 sigma_F, G2b <= 0.13% |
| `sep_dynamic/sep_gateG.json`, `sep_gateG_adjudication.json` | superseded first gate run (G1/G2a FAIL) + two-layer diagnosis: the offsets were a harness span defect (C4 guard direction missing) — span-only delta reproduces them on res0 with no live refit; guard-kept span reproduces z_lin to ~1e-16 |

Reproduction: `scripts/sep_common.py` (shared estimator core; seed
20260710), `scripts/sep_static_sensitivity.py`,
`scripts/sep_feasibility_gates.py`, `scripts/sep_joint_fit_10_8b.py`,
`scripts/anchor_resolution_10_8c.py`, `scripts/turn_search_10_8d.py`,
`scripts/sep_phase_marg_10_8e.py`, `scripts/make_gateG2_inputs.py`,
`scripts/gateG_adjudicate.py`
(all deterministic, artifacts-only), plus the WSL runtime scripts archived
verbatim with an artifact map in `scripts/wsl/README.md`.
