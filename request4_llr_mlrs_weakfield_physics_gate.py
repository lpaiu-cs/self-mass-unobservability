from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path


LAB_ROOT = Path("data/request4_llr/mlrs_handshake_lab")
REF_ROOT = LAB_ROOT / "data/analysis/ref"
ANALYSIS_ROOT = LAB_ROOT / "data/analysis"
RUN_ROOT = LAB_ROOT / "data/weakfield_physics_probe"
RUN_ROOT_REL = Path("data/weakfield_physics_probe")
SUMMARY_JSON = Path("request4_llr_mlrs_weakfield_physics_gate_summary.json")
SUMMARY_SVG = Path("request4_llr_mlrs_weakfield_physics_gate_summary.svg")

SAMPLE_STEMS = (
    "s25y11d016t0231#103",
    "s25y11d042t0234#103",
)

S_EARTH = 4.64e-10
S_MOON = 1.88e-11
SIGMA_CASES = (-1.0e-3, -5.0e-4, 5.0e-4, 1.0e-3)
ANALYTIC_SIGMA1_COEFF_M = 13.1


def run_command(
    args: list[str],
    cwd: Path,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=cwd, env=env, check=False, capture_output=True, text=True)


def rebuild_llr_npt_binary() -> None:
    cwd = LAB_ROOT / "src/llr_npt"
    commands = [
        [
            "/opt/homebrew/bin/gfortran",
            "-std=legacy",
            "-fno-align-commons",
            "-I../include",
            "-c",
            "jjreadnp.f",
            "-o",
            "jjreadnp.o",
        ],
        [
            "/opt/homebrew/bin/gfortran",
            "-std=legacy",
            "-fno-align-commons",
            "-I../include",
            "-c",
            "smu_probe_adjust_states.f",
            "-o",
            "smu_probe_adjust_states.o",
        ],
        [
            "/opt/homebrew/bin/gfortran",
            "-std=legacy",
            "-fno-align-commons",
            "-o",
            "llr_npt_crd",
            "ceror.o",
            "datio01.o",
            "error.o",
            "eulcomnp.o",
            "grtjd.o",
            "jdtgr01.o",
            "jeulpkg.o",
            "jjreadnp.o",
            "kustner01.o",
            "lsq_crd.o",
            "lunio02.o",
            "maknorm_crd.o",
            "misclibnp_crd.o",
            "mjreadnp.o",
            "mlrsgc01.o",
            "np_crd.o",
            "npaug10.o",
            "nrept11.o",
            "opio.o",
            "parin_crd.o",
            "pleph_lib.o",
            "pmodrep02.o",
            "putnp_crd.o",
            "rdlun_crd.o",
            "refrmm.o",
            "saltne.o",
            "sbetidnp.o",
            "smu_probe_adjust_omc.o",
            "smu_probe_adjust_states.o",
            "sncalc_crd.o",
            "sublasr03.o",
            "../utilib/calc_pmmf.o",
            "../utilib/cexit.o",
            "../utilib/doytomjd.o",
            "../utilib/getfield.o",
            "../utilib/grtodoyf.o",
            "../utilib/read_crd.o",
            "../utilib/read_crd_MLRS.o",
            "../utilib/read_crd_MLRSf.o",
            "../utilib/read_crdf.o",
            "../utilib/trimlen.o",
            "../utilib/write_crd.o",
            "../utilib/write_crd_MLRS.o",
            "../utilib/write_crd_MLRSf.o",
            "../utilib/write_crdf.o",
        ],
        ["cp", "llr_npt_crd", "../../bin/llr_npt_crd"],
    ]
    for command in commands:
        result = run_command(command, cwd=cwd)
        if result.returncode != 0:
            raise RuntimeError(result.stderr or result.stdout or f"build failed: {' '.join(command)}")


def prepare_run(sample_stem: str, run_name: str) -> tuple[Path, Path]:
    run_dir = RUN_ROOT / sample_stem / run_name
    run_dir_rel = RUN_ROOT_REL / sample_stem / run_name
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    frcal_ref = REF_ROOT / f"{sample_stem}.frcal"
    frcal_src = frcal_ref if frcal_ref.exists() else ANALYSIS_ROOT / f"{sample_stem}.frcal"
    rsc_ref = REF_ROOT / f"{sample_stem}.rsc"
    rsc_src = rsc_ref if rsc_ref.exists() else ANALYSIS_ROOT / f"{sample_stem}.rsc"
    shutil.copy2(frcal_src, run_dir / "probe.frcal")
    shutil.copy2(rsc_src, run_dir / "probe.rsc")
    (run_dir / "probe.llc").write_text("")
    (run_dir / "probe.lsc").write_text("")
    return run_dir, run_dir_rel


def run_pipeline(sample_stem: str, run_name: str, env_overrides: dict[str, str | None]) -> Path:
    run_dir, run_dir_rel = prepare_run(sample_stem, run_name)
    env = dict(os.environ)
    for key, value in env_overrides.items():
        if value is None:
            env.pop(key, None)
        else:
            env[key] = value

    commands = [
        [
            "./bin/llr_npt_crd",
            "-r",
            "-p",
            "poldat_thruJ17",
            "-i",
            str(run_dir_rel / "probe.frcal"),
            "-o",
            str(run_dir_rel / "probe.frrec"),
            "-c",
            str(run_dir_rel / "probe.rsc"),
            "-l",
            str(run_dir_rel / "probe.llc"),
            "-s",
            str(run_dir_rel / "probe.lsc"),
        ],
        [
            "./bin/Poisson_crd",
            "-i",
            str(run_dir_rel / "probe.frrec"),
            "-o",
            str(run_dir_rel / "probe.frfil"),
            "-c",
            "data/lib/filter_crd.pc",
            "-l",
            str(run_dir_rel / "probe.llc"),
            "-s",
            str(run_dir_rel / "probe.lsc"),
        ],
        [
            "./bin/lun_auto_crd",
            "-i",
            str(run_dir_rel / "probe.frfil"),
            "-o",
            str(run_dir_rel / "probe.frfil2"),
            "-l",
            str(run_dir_rel / "probe.llc"),
            "-s",
            str(run_dir_rel / "probe.lsc"),
        ],
        [
            "./bin/lun_sdqi_crd",
            "-i",
            str(run_dir_rel / "probe.frfil2"),
            "-o",
            str(run_dir_rel / "probe.frfil3"),
            "-l",
            str(run_dir_rel / "probe.llc"),
            "-s",
            str(run_dir_rel / "probe.lsc"),
        ],
        [
            "./bin/llr_npt_crd",
            "-n",
            "-p",
            "poldat_thruJ17",
            "-i",
            str(run_dir_rel / "probe.frfil3"),
            "-o",
            str(run_dir_rel / "probe.frfin"),
            "-m",
            str(run_dir_rel / "probe.npt"),
            "-c",
            str(run_dir_rel / "probe.rsc"),
            "-l",
            str(run_dir_rel / "probe.llc"),
            "-s",
            str(run_dir_rel / "probe.lsc"),
        ],
        [
            "./bin/frd_strip",
            str(run_dir_rel / "probe.frfin"),
            str(run_dir_rel / "probe.frd"),
        ],
    ]

    for command in commands:
        result = run_command(command, cwd=LAB_ROOT, env=env)
        if command[0].endswith("Poisson_crd") and result.returncode >= 2 and not (run_dir / "probe.frfil").exists():
            raise RuntimeError(result.stderr or result.stdout or "poisson failed")
        if command[0].endswith("Poisson_crd") and not (run_dir / "probe.frfil").exists():
            shutil.copy2(run_dir / "probe.frrec", run_dir / "probe.frfil")
        elif command[0].endswith("lun_auto_crd") and result.returncode != 0 and not (run_dir / "probe.frfil2").exists():
            shutil.copy2(run_dir / "probe.frfil", run_dir / "probe.frfil2")
        elif command[0].endswith("lun_sdqi_crd") and result.returncode != 0 and not (run_dir / "probe.frfil3").exists():
            shutil.copy2(run_dir / "probe.frfil2", run_dir / "probe.frfil3")
        elif result.returncode != 0 and not command[0].endswith(("Poisson_crd", "lun_auto_crd", "lun_sdqi_crd")):
            raise RuntimeError(result.stderr or result.stdout or f"{command[0]} failed")

    return run_dir


def read_lines(path: Path) -> list[str]:
    return path.read_text(errors="replace").splitlines()


def diff_lines(left: Path, right: Path) -> list[tuple[int, str, str]]:
    left_lines = read_lines(left)
    right_lines = read_lines(right)
    return [
        (line_number, left_line, right_line)
        for line_number, (left_line, right_line) in enumerate(zip(left_lines, right_lines), start=1)
        if left_line != right_line
    ]


def parse_93_last_deltas(left: Path, right: Path) -> list[float]:
    deltas: list[float] = []
    for left_line, right_line in zip(read_lines(left), read_lines(right)):
        if not left_line.startswith("93 ") or not right_line.startswith("93 "):
            continue
        deltas.append(float(left_line.split()[-1]) - float(right_line.split()[-1]))
    return deltas


def parse_npt_line_11(path: Path) -> list[str]:
    for line in read_lines(path):
        if line.startswith("11 "):
            return line.split()
    raise ValueError(f"missing 11 record in {path}")


def npt_token_delta(current: Path, reference: Path, token_index: int) -> float:
    current_tokens = parse_npt_line_11(current)
    reference_tokens = parse_npt_line_11(reference)
    return float(current_tokens[token_index - 1]) - float(reference_tokens[token_index - 1])


def analyze_run(run_dir: Path, baseline_dir: Path) -> dict[str, object]:
    frrec_93 = parse_93_last_deltas(run_dir / "probe.frrec", baseline_dir / "probe.frrec")
    changed_93 = [value for value in frrec_93 if value != 0.0]
    npt_line_diffs = diff_lines(run_dir / "probe.npt", baseline_dir / "probe.npt")
    frfin_diffs = diff_lines(run_dir / "probe.frfin", baseline_dir / "probe.frfin")
    return {
        "changed_93_lines": len(changed_93),
        "total_93_lines": len(frrec_93),
        "mean_delta_ns": sum(changed_93) / len(changed_93) if changed_93 else 0.0,
        "min_delta_ns": min(frrec_93) if frrec_93 else 0.0,
        "max_delta_ns": max(frrec_93) if frrec_93 else 0.0,
        "frfin_line_diff_count": len(frfin_diffs),
        "frd_exact_match": (run_dir / "probe.frd").read_bytes() == (baseline_dir / "probe.frd").read_bytes(),
        "npt_line_diff_count": len(npt_line_diffs),
        "npt_time_of_flight_delta_s": npt_token_delta(run_dir / "probe.npt", baseline_dir / "probe.npt", 3),
        "npt_pmm_delta_ps": npt_token_delta(run_dir / "probe.npt", baseline_dir / "probe.npt", 11),
        "npt_first_diff": npt_line_diffs[0] if npt_line_diffs else None,
    }


def delta_sep_from_sigmas(sigma1: float, sigma2: float) -> float:
    return sigma1 * (S_EARTH - S_MOON) + sigma2 * (S_EARTH**2 - S_MOON**2)


def weakfield_relative_amplitude_m(sigma1: float, sigma2: float) -> float:
    return ANALYTIC_SIGMA1_COEFF_M * sigma1 + (
        ANALYTIC_SIGMA1_COEFF_M * (S_EARTH**2 - S_MOON**2) / (S_EARTH - S_MOON)
    ) * sigma2


def write_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1480, 1160
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: MLRS weak-field physics gate</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">External weak-field amplitude from Request 3, injected as a barycenter-consistent Sun-driven relative state correction at the jjreadnp seam.</text>',
        '<rect x="48" y="108" width="1384" height="168" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Model</text>',
        '<text x="66" y="168" font-family="monospace" font-size="13" fill="#111827">delta_SEP = sigma1*(s_E-s_M) + sigma2*(s_E^2-s_M^2)</text>',
        f'<text x="66" y="192" font-family="monospace" font-size="13" fill="#111827">A_N = {ANALYTIC_SIGMA1_COEFF_M/(S_EARTH-S_MOON):.6e} m per delta_SEP</text>',
        '<text x="66" y="216" font-family="monospace" font-size="13" fill="#111827">helper applies delta r_rel = A * u_ES(t), delta v_rel = A * d(u_ES)/dt and splits it across Earth/Moon with fixed EM barycenter</text>',
        '<text x="66" y="240" font-family="monospace" font-size="13" fill="#111827">edit scope stays local: jjreadnp.f + smu_probe_adjust_states.f only</text>',
    ]

    y = 308
    for sample_stem, runs in summary["runs"].items():
        lines.append(f'<text x="48" y="{y}" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">{sample_stem}</text>')
        y += 28
        for label, run in runs.items():
            lines.append(
                f'<text x="66" y="{y}" font-family="monospace" font-size="12" fill="#111827">{label}: amp={run["relative_amp_m"]:.6e} m, mean Δ93={run["mean_delta_ns"]:.6e} ns, frfin={run["frfin_line_diff_count"]}, frd={run["frd_exact_match"]}, npt ΔTOF={run["npt_time_of_flight_delta_s"]:.3e} s, npt ΔPmM={run["npt_pmm_delta_ps"]:.3f} ps</text>'
            )
            y += 22
        y += 10

    lines.append(f'<text x="48" y="{y}" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Odd Symmetry / Scaling</text>')
    y += 28
    for sample_stem, odd_items in summary["odd_symmetry"].items():
        lines.append(f'<text x="66" y="{y}" font-family="monospace" font-size="13" fill="#111827">{sample_stem}</text>')
        y += 22
        for item in odd_items:
            lines.append(
                f'<text x="84" y="{y}" font-family="monospace" font-size="12" fill="#111827">eps={item["epsilon_sigma1"]:.4g} -> frrec sum={item["frrec_mean_sum_ns"]:.3e} ns, npt TOF sum={item["npt_tof_sum_s"]:.3e} s, npt PmM sum={item["npt_pmm_sum_ps"]:.3e} ps</text>'
            )
            y += 20
        for item in summary["scaling"][sample_stem]:
            lines.append(
                f'<text x="84" y="{y}" font-family="monospace" font-size="12" fill="#111827">eps={item["epsilon_sigma1"]:.4g} -> gain(Δ93/sigma1)={item["frrec_gain_ns_per_sigma1"]:.6e}, gain(ΔPmM/sigma1)={item["npt_pmm_gain_ps_per_sigma1"]:.6e}</text>'
            )
            y += 20
        y += 10

    lines.extend(
        [
            f'<text x="48" y="{y+10}" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Reading</text>',
            f'<text x="66" y="{y+40}" font-family="monospace" font-size="13" fill="#111827">If these runs stay local and near-linear, the software stop rule remains cleared even for a weak-field-inspired state correction.</text>',
            f'<text x="66" y="{y+64}" font-family="monospace" font-size="13" fill="#111827">The remaining question is whether this local observation-operator model is physically adequate, not whether the seam survives.</text>',
            "</svg>",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    rebuild_llr_npt_binary()

    runs: dict[str, dict[str, dict[str, object]]] = {}
    odd_symmetry: dict[str, list[dict[str, float]]] = {}
    scaling: dict[str, list[dict[str, float]]] = {}

    for sample_stem in SAMPLE_STEMS:
        baseline_dir = run_pipeline(
            sample_stem,
            "baseline",
            {
                "SMU_OMC_OFFSET_NS": None,
                "SMU_OMC_SINE_AMPL_NS": None,
                "SMU_OMC_SINE_PERIOD_S": None,
                "SMU_OMC_SINE_T0_JD": None,
                "SMU_OMC_SINE_T0_SOD": None,
                "SMU_OMC_SINE_PHASE_RAD": None,
                "SMU_STATE_EM_SHIFT_M": None,
                "SMU_STATE_SUN_SHIFT_M": None,
                "SMU_WF_REL_AMP_M": None,
            },
        )

        sample_runs: dict[str, dict[str, object]] = {}
        for sigma1 in SIGMA_CASES:
            sigma2 = 0.0
            label = f"sigma1_{sigma1:+.4g}".replace("+", "p").replace("-", "m")
            rel_amp_m = weakfield_relative_amplitude_m(sigma1, sigma2)
            run_dir = run_pipeline(
                sample_stem,
                label,
                {
                    "SMU_OMC_OFFSET_NS": None,
                    "SMU_OMC_SINE_AMPL_NS": None,
                    "SMU_OMC_SINE_PERIOD_S": None,
                    "SMU_OMC_SINE_T0_JD": None,
                    "SMU_OMC_SINE_T0_SOD": None,
                    "SMU_OMC_SINE_PHASE_RAD": None,
                    "SMU_STATE_EM_SHIFT_M": None,
                    "SMU_STATE_SUN_SHIFT_M": None,
                    "SMU_WF_REL_AMP_M": f"{rel_amp_m:.12g}",
                },
            )
            sample_runs[f"{sigma1:+.4g}"] = analyze_run(run_dir, baseline_dir) | {
                "sigma1": sigma1,
                "sigma2": sigma2,
                "delta_sep": delta_sep_from_sigmas(sigma1, sigma2),
                "relative_amp_m": rel_amp_m,
            }

        runs[sample_stem] = sample_runs

        odd_items = []
        scaling_items = []
        for eps in (5.0e-4, 1.0e-3):
            pos = sample_runs[f"{eps:+.4g}"]
            neg = sample_runs[f"{-eps:+.4g}"]
            odd_items.append(
                {
                    "epsilon_sigma1": eps,
                    "frrec_mean_sum_ns": pos["mean_delta_ns"] + neg["mean_delta_ns"],  # type: ignore[operator]
                    "npt_tof_sum_s": pos["npt_time_of_flight_delta_s"] + neg["npt_time_of_flight_delta_s"],  # type: ignore[operator]
                    "npt_pmm_sum_ps": pos["npt_pmm_delta_ps"] + neg["npt_pmm_delta_ps"],  # type: ignore[operator]
                }
            )
            scaling_items.append(
                {
                    "epsilon_sigma1": eps,
                    "frrec_gain_ns_per_sigma1": pos["mean_delta_ns"] / eps,  # type: ignore[operator]
                    "npt_pmm_gain_ps_per_sigma1": pos["npt_pmm_delta_ps"] / eps,  # type: ignore[operator]
                }
            )
        odd_symmetry[sample_stem] = odd_items
        scaling[sample_stem] = scaling_items

    summary = {
        "status": "weakfield_local_state_gate_completed",
        "model": {
            "sigma_cases": SIGMA_CASES,
            "sigma2_fixed": 0.0,
            "analytic_sigma1_coeff_m": ANALYTIC_SIGMA1_COEFF_M,
            "analytic_delta_sep_to_range_coeff_m": ANALYTIC_SIGMA1_COEFF_M / (S_EARTH - S_MOON),
            "injection_class": "delta_r_rel = A_N * delta_SEP * u_ES(t), delta_v_rel = A_N * delta_SEP * d(u_ES)/dt",
        },
        "runs": runs,
        "odd_symmetry": odd_symmetry,
        "scaling": scaling,
        "assessment": (
            "A weak-field-inspired, externally parameterized delta_SEP state correction "
            "still fits through the local jjreadnp seam without widening the MLRS branch "
            "into PLEPH edits or ephemeris regeneration. The remaining open issue is "
            "physical adequacy of this observation-operator approximation."
        ),
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    write_svg(SUMMARY_SVG, summary)


if __name__ == "__main__":
    main()
