# WSL runtime scripts of record (10.8f C11)

These scripts ran INSIDE the WSL Nutimo runtime (`~/work/nutimo_pilot/run_planetGR`,
see `notes/` pre-registrations and the memory note on the pilot environment);
they are the producers of every live-measurement artifact committed under
`request10_external/` and `request10_external/sep_dynamic/`. They are archived
here verbatim for reproducibility — they do NOT run from the repository tree.

Launch pattern (from Windows; scripts are staged through `tr -d '\r'` because
the repo may carry CRLF):

```bash
export LD_LIBRARY_PATH="$HOME/work/nutimo_pilot/nutimo/src"
export TEMPO2="$HOME/work/nutimo_pilot/install/third_party/tempo2"
export PYTHONPATH="$HOME/work/nutimo_pilot/nutimo/src"
cd "$HOME/work/nutimo_pilot/run_planetGR"
# copy the script's declared repo inputs into the run dir first, then:
tr -d '\r' < wsl_<name>.py > /tmp/wsl_<name>.py
python3 /tmp/wsl_<name>.py 2>&1 | tee /tmp/<name>.log
# copy the script's declared outputs back into request10_external/
```

Long jobs were launched detached (`setsid nohup ... &`). External inputs not
in this repository: the runtime tree itself (Nutimo build + tempo2 +
DE430 ephemeris), the TOA/parfile pair
(`0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim`,
`parfile-planetGR-max-bestfit`), and the published-fit MCMC chain
`Analysis/Synthesis/MCMC-planetGR.npz` (used only for step sizes and the
10.8c posterior widths; the values consumed are cached in
`steps_provenance.json` and `scripts/anchor_resolution_10_8c.py`).
The dynamic-drive integrator patch is committed at
`sep_dynamic/patches/sepdyn.patch`.

| script | request | committed artifact(s) of record | status |
|---|---|---|---|
| `wsl_rederive_planet_jac.py` | 10.7b | `finite_jacobian_v2.npy`, `finite_jacobian_v2_meta.json` | of record |
| `wsl_lincheck.py` | 10.7b | `jointfit_linearization_check.json` | of record |
| `wsl_lindiag.py` | 10.7b | `jointfit_lindiag.json` | of record |
| `wsl_gateLA.py` | 10.7c | `jointfit_gateLA.json` (inputs: `gateLA_params.json`, `xb52/xc52.npy`, `jointfit_dtheta52_trunc.npy`) | of record |
| `wsl_gateD.py` | 10.7d | `jointfit_gateD.json` (inputs: `gateD_columns.npz`, `gateD_params.json`) | of record |
| `wsl_sepF0.py` | 10.8 F0 | `sep_dynamic/col_SEP_D.npz`, `sep_dynamic/sep_F0_report.json` | of record |
| `wsl_sepR5.py` | 10.8 R5 | `sep_dynamic/col_R5_*.npz`, `sep_dynamic/sep_R5_report.json` | of record |
| `wsl_t2cols.py` | 10.8 T2 | `sep_dynamic/sep_t2_results.json`, probe columns | superseded (patched-baseline integrity + first drive probes) |
| `wsl_t2cols_v2.py` | 10.8 T2 | `sep_dynamic/sep_t2_results_v2.json` | superseded (spectral-scan calibration attempt; stopped per its own rule) |
| `wsl_t2cols_v3.py` | 10.8 T2 | none | superseded iteration of v2; no artifact of record |
| `wsl_t2cols_v4.py` | 10.8 T2 | `sep_dynamic/sep_t2_results_v4.json`, `sep_dynamic/sep_dynamic_columns.npz`, `sep_dynamic/sep_ramp_probe.npz` | **calibration + production columns of record** (quasi-static ramp: integrator time in TIMEDAYS; integral-kernel 1/2 signature 0.502/0.0219) |
| `wsl_turnV.py` | 10.8d V | `sep_dynamic/turnV_columns.npz`, `sep_dynamic/turnV_meta.json` | of record |
| `wsl_gateG.py` | 10.8b G / 10.8f R6 | `sep_dynamic/sep_gateG.json` (inputs: `sep_dynamic/gateG_inputs.npz`, `sep_dynamic/gateG_params.json`) | of record — CAUSAL-template re-run; the pre-C1 anticausal run is in git history (52153b3) |

Adjudication of the 10.8f R6 gate outcome (G1/G2a FAIL, G2b PASS, common-mode
null-offset diagnosis): `scripts/gateG_adjudicate.py` →
`sep_dynamic/sep_gateG_adjudication.json`.
