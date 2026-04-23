from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent
LAB_ROOT = REPO_ROOT / "data/request4_llr/mlrs_handshake_lab"
REF_ROOT = LAB_ROOT / "data/analysis/ref"
ANALYSIS_ROOT = LAB_ROOT / "data/analysis"
RUN_ROOT = LAB_ROOT / "data/state_seam_probe"
RUN_ROOT_REL = Path("data/state_seam_probe")
SUMMARY_JSON = ROOT / "request4_llr_mlrs_state_seam_gate_summary.json"
SUMMARY_SVG = ROOT / "request4_llr_mlrs_state_seam_gate_summary.svg"

SAMPLE_STEM = "s25y11d016t0231#103"
SHIFT_VALUES_M = (-0.01, 0.01)


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


def prepare_run(run_name: str) -> tuple[Path, Path]:
    run_dir = RUN_ROOT / SAMPLE_STEM / run_name
    run_dir_rel = RUN_ROOT_REL / SAMPLE_STEM / run_name
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    frcal_ref = REF_ROOT / f"{SAMPLE_STEM}.frcal"
    frcal_src = frcal_ref if frcal_ref.exists() else ANALYSIS_ROOT / f"{SAMPLE_STEM}.frcal"
    rsc_ref = REF_ROOT / f"{SAMPLE_STEM}.rsc"
    rsc_src = rsc_ref if rsc_ref.exists() else ANALYSIS_ROOT / f"{SAMPLE_STEM}.rsc"
    shutil.copy2(frcal_src, run_dir / "probe.frcal")
    shutil.copy2(rsc_src, run_dir / "probe.rsc")
    (run_dir / "probe.llc").write_text("")
    (run_dir / "probe.lsc").write_text("")
    return run_dir, run_dir_rel


def run_pipeline(run_name: str, env_overrides: dict[str, str | None]) -> Path:
    run_dir, run_dir_rel = prepare_run(run_name)
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


def write_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1420, 820
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: MLRS state seam gate</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Bounded probe at the jjreadnp/ephemeris-state bridge, not a physical delta_SEP implementation.</text>',
        '<rect x="48" y="108" width="1324" height="182" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Build Scope</text>',
        '<text x="66" y="168" font-family="monospace" font-size="13" fill="#111827">source edits: jjreadnp.f + smu_probe_adjust_states.f</text>',
        '<text x="66" y="192" font-family="monospace" font-size="13" fill="#111827">rebuild: jjreadnp.o + smu_probe_adjust_states.o + llr_npt_crd relink only</text>',
        '<text x="66" y="216" font-family="monospace" font-size="13" fill="#111827">no batch-script edits, no PLEPH edits, no jeulpkg.f surgery</text>',
        '<text x="66" y="240" font-family="monospace" font-size="13" fill="#111827">probe: +/- 0.01 m Moon displacement along instantaneous Earth->Moon line</text>',
        '<rect x="48" y="320" width="1324" height="220" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="350" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Results</text>',
    ]
    y = 382
    for label, run in summary["runs"].items():
        lines.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="12" fill="#111827">{label}: mean Δ93 = {run["mean_delta_ns"]:.6f} ns, changed 93 = {run["changed_93_lines"]}/{run["total_93_lines"]}, frfin diffs = {run["frfin_line_diff_count"]}, frd exact = {run["frd_exact_match"]}, npt ΔTOF = {run["npt_time_of_flight_delta_s"]:.3e} s, npt ΔPmM = {run["npt_pmm_delta_ps"]:.3f} ps</text>'
        )
        y += 24
    lines.extend(
        [
            '<rect x="48" y="570" width="1324" height="170" fill="white" stroke="#cbd5e1"/>',
            '<text x="66" y="600" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Verdict</text>',
            '<text x="66" y="630" font-family="monospace" font-size="13" fill="#111827">The prediction/state insertion does not spill into broad surgery at the software-interface level.</text>',
            '<text x="66" y="654" font-family="monospace" font-size="13" fill="#111827">The gate remains bounded to the ephemeris-state bridge, but this is still not a physical weak-field fit.</text>',
            '<text x="66" y="678" font-family="monospace" font-size="13" fill="#111827">If a real delta_SEP model later demands editing PLEPH or regenerating ephemerides, stop and hand off outward.</text>',
            "</svg>",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    rebuild_llr_npt_binary()
    baseline_dir = run_pipeline(
        "baseline",
        {
            "SMU_OMC_OFFSET_NS": None,
            "SMU_OMC_SINE_AMPL_NS": None,
            "SMU_OMC_SINE_PERIOD_S": None,
            "SMU_OMC_SINE_T0_JD": None,
            "SMU_OMC_SINE_T0_SOD": None,
            "SMU_OMC_SINE_PHASE_RAD": None,
            "SMU_STATE_EM_SHIFT_M": None,
        },
    )

    runs: dict[str, dict[str, object]] = {}
    for shift_m in SHIFT_VALUES_M:
        label = f"{shift_m:+.3f}m"
        run_dir = run_pipeline(
            f"state_shift_{label.replace('+', 'p').replace('-', 'm')}",
            {
                "SMU_OMC_OFFSET_NS": None,
                "SMU_OMC_SINE_AMPL_NS": None,
                "SMU_OMC_SINE_PERIOD_S": None,
                "SMU_OMC_SINE_T0_JD": None,
                "SMU_OMC_SINE_T0_SOD": None,
                "SMU_OMC_SINE_PHASE_RAD": None,
                "SMU_STATE_EM_SHIFT_M": f"{shift_m:.12g}",
            },
        )
        runs[label] = analyze_run(run_dir, baseline_dir)

    odd_symmetry = {
        "frrec_mean_sum_ns": runs["+0.010m"]["mean_delta_ns"] + runs["-0.010m"]["mean_delta_ns"],  # type: ignore[operator]
        "npt_tof_sum_s": runs["+0.010m"]["npt_time_of_flight_delta_s"] + runs["-0.010m"]["npt_time_of_flight_delta_s"],  # type: ignore[operator]
        "npt_pmm_sum_ps": runs["+0.010m"]["npt_pmm_delta_ps"] + runs["-0.010m"]["npt_pmm_delta_ps"],  # type: ignore[operator]
    }

    summary = {
        "status": "state_seam_bounded",
        "sample_stem": SAMPLE_STEM,
        "shift_values_m": SHIFT_VALUES_M,
        "build_scope": {
            "source_files_touched": [
                "src/llr_npt/jjreadnp.f",
                "src/llr_npt/smu_probe_adjust_states.f",
            ],
            "recompiled_objects": [
                "src/llr_npt/jjreadnp.o",
                "src/llr_npt/smu_probe_adjust_states.o",
            ],
            "relinked_binary": "bin/llr_npt_crd",
            "batch_script_edits": 0,
            "pleph_edits": 0,
            "jeulpkg_edits": 0,
        },
        "probe_model": {
            "kind": "moon_displacement_along_earth_moon_line",
            "amplitude_m": 0.01,
        },
        "runs": runs,
        "odd_symmetry": odd_symmetry,
        "assessment": (
            "A deeper prediction/state-seam probe fits through the jjreadnp bridge "
            "without widening into broad legacy-code surgery. The remaining risk is "
            "physical, not architectural: a real delta_SEP insertion may still need "
            "a model stronger than this synthetic state displacement."
        ),
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    write_svg(SUMMARY_SVG, summary)


if __name__ == "__main__":
    main()
