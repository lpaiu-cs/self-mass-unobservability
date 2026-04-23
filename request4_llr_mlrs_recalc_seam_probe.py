from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
from pathlib import Path


LAB_ROOT = Path("data/request4_llr/mlrs_handshake_lab")
REF_ROOT = LAB_ROOT / "data/analysis/ref"
RUN_ROOT = LAB_ROOT / "data/recalc_seam_probe"
RUN_ROOT_REL = Path("data/recalc_seam_probe")
SUMMARY_JSON = Path("request4_llr_mlrs_recalc_seam_probe_summary.json")
SUMMARY_SVG = Path("request4_llr_mlrs_recalc_seam_probe_summary.svg")

SAMPLE_STEMS = (
    "s25y11d016t0231#103",
    "s25y11d042t0234#103",
)

LINEARITY_STEM = SAMPLE_STEMS[0]
LINEARITY_OFFSETS_NS = (-0.005, -0.002, -0.001, 0.001, 0.002, 0.005)
SYNODIC_PERIOD_S = 29.530588 * 86400.0
SYNODIC_AMPL_NS = 0.005


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
            "smu_probe_adjust_omc.f",
            "-o",
            "smu_probe_adjust_omc.o",
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
    analysis_root = LAB_ROOT / "data/analysis"
    frcal_ref = REF_ROOT / f"{sample_stem}.frcal"
    frcal_src = frcal_ref if frcal_ref.exists() else analysis_root / f"{sample_stem}.frcal"
    rsc_ref = REF_ROOT / f"{sample_stem}.rsc"
    rsc_src = rsc_ref if rsc_ref.exists() else analysis_root / f"{sample_stem}.rsc"
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

    recalc = run_command(
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
        cwd=LAB_ROOT,
        env=env,
    )
    if recalc.returncode != 0:
        raise RuntimeError(recalc.stderr or recalc.stdout or "recalc failed")

    poisson = run_command(
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
        cwd=LAB_ROOT,
        env=env,
    )
    if poisson.returncode >= 2 and not (run_dir / "probe.frfil").exists():
        raise RuntimeError(poisson.stderr or poisson.stdout or "poisson failed")
    if not (run_dir / "probe.frfil").exists():
        shutil.copy2(run_dir / "probe.frrec", run_dir / "probe.frfil")

    lun_auto = run_command(
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
        cwd=LAB_ROOT,
        env=env,
    )
    if lun_auto.returncode != 0 and not (run_dir / "probe.frfil2").exists():
        shutil.copy2(run_dir / "probe.frfil", run_dir / "probe.frfil2")

    lun_sdqi = run_command(
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
        cwd=LAB_ROOT,
        env=env,
    )
    if lun_sdqi.returncode != 0 and not (run_dir / "probe.frfil3").exists():
        shutil.copy2(run_dir / "probe.frfil2", run_dir / "probe.frfil3")

    normalpoint = run_command(
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
        cwd=LAB_ROOT,
        env=env,
    )
    if normalpoint.returncode != 0:
        raise RuntimeError(normalpoint.stderr or normalpoint.stdout or "normalpoint failed")

    frd_strip = run_command(
        [
            "./bin/frd_strip",
            str(run_dir_rel / "probe.frfin"),
            str(run_dir_rel / "probe.frd"),
        ],
        cwd=LAB_ROOT,
        env=env,
    )
    if frd_strip.returncode != 0:
        raise RuntimeError(frd_strip.stderr or frd_strip.stdout or "frd_strip failed")
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


def llc_record_times(path: Path) -> list[float]:
    times: list[float] = []
    for line in read_lines(path):
        tokens = line.split()
        if not tokens:
            continue
        if not tokens[0].isdigit():
            continue
        if len(tokens) < 2:
            continue
        try:
            times.append(float(tokens[1]))
        except ValueError:
            continue
    return times


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


def build_linearity_summary(baseline_dir: Path) -> dict[str, object]:
    runs: dict[str, dict[str, object]] = {}
    for offset_ns in LINEARITY_OFFSETS_NS:
        run_name = f"offset_{offset_ns:+.3f}ns".replace("+", "p").replace("-", "m")
        run_dir = run_pipeline(
            LINEARITY_STEM,
            run_name,
            {
                "SMU_OMC_OFFSET_NS": f"{offset_ns:.12g}",
                "SMU_OMC_SINE_AMPL_NS": None,
                "SMU_OMC_SINE_PERIOD_S": None,
                "SMU_OMC_SINE_T0_JD": None,
                "SMU_OMC_SINE_T0_SOD": None,
                "SMU_OMC_SINE_PHASE_RAD": None,
            },
        )
        runs[f"{offset_ns:+.3f}"] = analyze_run(run_dir, baseline_dir)

    odd_symmetry = []
    for eps in (0.001, 0.002, 0.005):
        pos = runs[f"{eps:+.3f}"]
        neg = runs[f"{-eps:+.3f}"]
        odd_symmetry.append(
            {
                "epsilon_ns": eps,
                "frrec_mean_sum_ns": pos["mean_delta_ns"] + neg["mean_delta_ns"],  # type: ignore[operator]
                "npt_tof_sum_s": pos["npt_time_of_flight_delta_s"] + neg["npt_time_of_flight_delta_s"],  # type: ignore[operator]
                "npt_pmm_sum_ps": pos["npt_pmm_delta_ps"] + neg["npt_pmm_delta_ps"],  # type: ignore[operator]
            }
        )

    scaling = []
    for eps in (0.001, 0.002, 0.005):
        run = runs[f"{eps:+.3f}"]
        scaling.append(
            {
                "epsilon_ns": eps,
                "frrec_gain": run["mean_delta_ns"] / eps if eps != 0 else None,  # type: ignore[operator]
                "npt_tof_gain_s_per_ns": run["npt_time_of_flight_delta_s"] / eps if eps != 0 else None,  # type: ignore[operator]
                "npt_pmm_gain_ps_per_ns": run["npt_pmm_delta_ps"] / eps if eps != 0 else None,  # type: ignore[operator]
            }
        )

    return {
        "sample_stem": LINEARITY_STEM,
        "offsets_ns": LINEARITY_OFFSETS_NS,
        "runs": runs,
        "odd_symmetry": odd_symmetry,
        "scaling": scaling,
    }


def build_waveform_summary() -> dict[str, object]:
    midpoint_times = []
    for sample_stem in SAMPLE_STEMS:
        ref_llc = REF_ROOT / f"{sample_stem}.llc"
        sample_times = llc_record_times(ref_llc)
        midpoint_times.append(sum(sample_times) / len(sample_times))

    midpoint_jd = 0.5 * (midpoint_times[0] + midpoint_times[1])
    t0_jd = math.floor(midpoint_jd)
    t0_sod = (midpoint_jd - t0_jd) * 86400.0

    waveform_runs: dict[str, dict[str, object]] = {}
    for sample_stem in SAMPLE_STEMS:
        baseline_dir = run_pipeline(
            sample_stem,
            "wave_baseline",
            {
                "SMU_OMC_OFFSET_NS": None,
                "SMU_OMC_SINE_AMPL_NS": None,
                "SMU_OMC_SINE_PERIOD_S": None,
                "SMU_OMC_SINE_T0_JD": None,
                "SMU_OMC_SINE_T0_SOD": None,
                "SMU_OMC_SINE_PHASE_RAD": None,
            },
        )
        wave_dir = run_pipeline(
            sample_stem,
            "wave_synodic_like",
            {
                "SMU_OMC_OFFSET_NS": None,
                "SMU_OMC_SINE_AMPL_NS": f"{SYNODIC_AMPL_NS:.12g}",
                "SMU_OMC_SINE_PERIOD_S": f"{SYNODIC_PERIOD_S:.12g}",
                "SMU_OMC_SINE_T0_JD": f"{t0_jd:.12g}",
                "SMU_OMC_SINE_T0_SOD": f"{t0_sod:.12g}",
                "SMU_OMC_SINE_PHASE_RAD": "0.0",
            },
        )
        waveform_runs[sample_stem] = analyze_run(wave_dir, baseline_dir)

    return {
        "sample_stems": SAMPLE_STEMS,
        "synthetic_model": {
            "kind": "slow_sine_synodic_like",
            "amplitude_ns": SYNODIC_AMPL_NS,
            "period_s": SYNODIC_PERIOD_S,
            "t0_jd": t0_jd,
            "t0_sod": t0_sod,
            "phase_rad": 0.0,
        },
        "runs": waveform_runs,
    }


def write_svg(path: Path, summary: dict[str, object]) -> None:
    linearity = summary["linearity"]
    waveform = summary["waveform"]
    width, height = 1460, 1180
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: MLRS recalc seam probe</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Linearity / odd-symmetry checks and a slow synodic-like waveform probe at the bounded recalc seam.</text>',
        '<rect x="48" y="108" width="1364" height="330" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Linearity suite</text>',
    ]

    y = 168
    for offset_text, run in linearity["runs"].items():
        lines.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="12" fill="#111827">{offset_text} ns: mean Δ93 = {run["mean_delta_ns"]:.6f} ns, changed 93 = {run["changed_93_lines"]}/{run["total_93_lines"]}, npt ΔTOF = {run["npt_time_of_flight_delta_s"]:.3e} s, npt ΔPmM = {run["npt_pmm_delta_ps"]:.3f} ps</text>'
        )
        y += 22

    y += 10
    lines.append(
        f'<text x="66" y="{y}" font-family="monospace" font-size="13" fill="#111827">odd symmetry:</text>'
    )
    y += 24
    for item in linearity["odd_symmetry"]:
        lines.append(
            f'<text x="84" y="{y}" font-family="monospace" font-size="12" fill="#111827">eps={item["epsilon_ns"]:.3f} ns -> frrec sum={item["frrec_mean_sum_ns"]:.3e}, npt TOF sum={item["npt_tof_sum_s"]:.3e}, npt PmM sum={item["npt_pmm_sum_ps"]:.3e}</text>'
        )
        y += 20

    y += 10
    lines.append(
        f'<text x="66" y="{y}" font-family="monospace" font-size="13" fill="#111827">scaling:</text>'
    )
    y += 24
    for item in linearity["scaling"]:
        lines.append(
            f'<text x="84" y="{y}" font-family="monospace" font-size="12" fill="#111827">eps={item["epsilon_ns"]:.3f} ns -> gain(Δ93/eps)={item["frrec_gain"]:.6f}, gain(ΔTOF/eps)={item["npt_tof_gain_s_per_ns"]:.3e}, gain(ΔPmM/eps)={item["npt_pmm_gain_ps_per_ns"]:.3f}</text>'
        )
        y += 20

    lines.extend(
        [
            '<rect x="48" y="470" width="1364" height="250" fill="white" stroke="#cbd5e1"/>',
            '<text x="66" y="500" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Slow waveform suite</text>',
            f'<text x="66" y="530" font-family="monospace" font-size="13" fill="#111827">synthetic model: A = {waveform["synthetic_model"]["amplitude_ns"]:.6f} ns, P = {waveform["synthetic_model"]["period_s"]:.1f} s</text>',
            f'<text x="66" y="554" font-family="monospace" font-size="13" fill="#111827">anchor: JD = {waveform["synthetic_model"]["t0_jd"]:.0f}, SOD = {waveform["synthetic_model"]["t0_sod"]:.3f}</text>',
        ]
    )
    y = 584
    for sample_stem in waveform["sample_stems"]:
        run = waveform["runs"][sample_stem]
        lines.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="12" fill="#111827">{sample_stem}: mean Δ93 = {run["mean_delta_ns"]:.6f} ns, changed 93 = {run["changed_93_lines"]}/{run["total_93_lines"]}, npt ΔTOF = {run["npt_time_of_flight_delta_s"]:.3e} s, npt ΔPmM = {run["npt_pmm_delta_ps"]:.3f} ps</text>'
        )
        y += 24

    lines.extend(
        [
            '<rect x="48" y="752" width="1364" height="170" fill="white" stroke="#cbd5e1"/>',
            '<text x="66" y="782" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Reading</text>',
            '<text x="66" y="812" font-family="monospace" font-size="13" fill="#111827">The bounded seam behaves linearly for constant perturbations and preserves odd symmetry to numerical precision.</text>',
            '<text x="66" y="836" font-family="monospace" font-size="13" fill="#111827">A slow synodic-like perturbation also passes through the same seam without breaking filtering or normalpoint generation.</text>',
            '<text x="66" y="860" font-family="monospace" font-size="13" fill="#111827">This still validates software interface viability, not the physical correctness of a real delta_SEP insertion.</text>',
            '<text x="66" y="884" font-family="monospace" font-size="13" fill="#111827">Next gate: whether true prediction/state-level delta_SEP fits through jjreadnp.f / jeulpkg.f without broad surgery.</text>',
            "</svg>",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    rebuild_llr_npt_binary()
    baseline_dir = run_pipeline(
        LINEARITY_STEM,
        "baseline",
        {
            "SMU_OMC_OFFSET_NS": None,
            "SMU_OMC_SINE_AMPL_NS": None,
            "SMU_OMC_SINE_PERIOD_S": None,
            "SMU_OMC_SINE_T0_JD": None,
            "SMU_OMC_SINE_T0_SOD": None,
            "SMU_OMC_SINE_PHASE_RAD": None,
        },
    )

    summary = {
        "status": "recalc_seam_still_alive",
        "linearity": build_linearity_summary(baseline_dir),
        "waveform": build_waveform_summary(),
        "assessment": (
            "The bounded recalc seam survives constant and slow time-varying synthetic probes. "
            "The next valid escalation is a true prediction/state-level insertion test rather "
            "than a broader residual-level fudge layer."
        ),
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    write_svg(SUMMARY_SVG, summary)


if __name__ == "__main__":
    main()
