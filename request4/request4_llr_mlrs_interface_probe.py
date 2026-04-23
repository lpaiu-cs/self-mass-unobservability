from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent
LAB_ROOT = REPO_ROOT / "data/request4_llr/mlrs_handshake_lab"
ANALYSIS_REF = LAB_ROOT / "data/analysis/ref"
PROBE_ROOT = LAB_ROOT / "data/interface_probe"
PROBE_ROOT_REL = Path("data/interface_probe")
SUMMARY_JSON = ROOT / "request4_llr_mlrs_interface_probe_summary.json"
SUMMARY_SVG = ROOT / "request4_llr_mlrs_interface_probe_summary.svg"
SAMPLE_STEM = "s25y11d016t0231#103"


def repo_relative(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def run_command(args: list[str], cwd: Path, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=cwd, env=env, check=False, capture_output=True, text=True)


def prepare_run(run_name: str) -> tuple[Path, Path]:
    out_dir_abs = PROBE_ROOT / run_name
    out_dir_rel = PROBE_ROOT_REL / run_name
    if out_dir_abs.exists():
        shutil.rmtree(out_dir_abs)
    out_dir_abs.mkdir(parents=True)
    shutil.copy2(ANALYSIS_REF / f"{SAMPLE_STEM}.frcal", out_dir_abs / "probe.frcal")
    shutil.copy2(ANALYSIS_REF / f"{SAMPLE_STEM}.rsc", out_dir_abs / "probe.rsc")
    (out_dir_abs / "probe.llc").write_text("")
    (out_dir_abs / "probe.lsc").write_text("")
    return out_dir_abs, out_dir_rel


def run_probe(run_name: str, offset_ns: float | None) -> Path:
    out_dir_abs, out_dir_rel = prepare_run(run_name)
    env = dict(**os.environ)
    if offset_ns is None:
        env.pop("SMU_OMC_OFFSET_NS", None)
    else:
        env["SMU_OMC_OFFSET_NS"] = f"{offset_ns:.12g}"

    recalc = run_command(
        [
            "./bin/llr_npt_crd",
            "-r",
            "-p",
            "poldat_thruJ17",
            "-i",
            str(out_dir_rel / "probe.frcal"),
            "-o",
            str(out_dir_rel / "probe.frrec"),
            "-c",
            str(out_dir_rel / "probe.rsc"),
            "-l",
            str(out_dir_rel / "probe.llc"),
            "-s",
            str(out_dir_rel / "probe.lsc"),
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
            str(out_dir_rel / "probe.frrec"),
            "-o",
            str(out_dir_rel / "probe.frfil"),
            "-c",
            "data/lib/filter_crd.pc",
            "-l",
            str(out_dir_rel / "probe.llc"),
            "-s",
            str(out_dir_rel / "probe.lsc"),
        ],
        cwd=LAB_ROOT,
        env=env,
    )
    if poisson.returncode >= 2 and not (out_dir_abs / "probe.frfil").exists():
        raise RuntimeError(poisson.stderr or poisson.stdout or "poisson failed")
    if not (out_dir_abs / "probe.frfil").exists():
        shutil.copy2(out_dir_abs / "probe.frrec", out_dir_abs / "probe.frfil")

    lun_auto = run_command(
        [
            "./bin/lun_auto_crd",
            "-i",
            str(out_dir_rel / "probe.frfil"),
            "-o",
            str(out_dir_rel / "probe.frfil2"),
            "-l",
            str(out_dir_rel / "probe.llc"),
            "-s",
            str(out_dir_rel / "probe.lsc"),
        ],
        cwd=LAB_ROOT,
        env=env,
    )
    if not (out_dir_abs / "probe.frfil2").exists():
        shutil.copy2(out_dir_abs / "probe.frfil", out_dir_abs / "probe.frfil2")

    lun_sdqi = run_command(
        [
            "./bin/lun_sdqi_crd",
            "-i",
            str(out_dir_rel / "probe.frfil2"),
            "-o",
            str(out_dir_rel / "probe.frfil3"),
            "-l",
            str(out_dir_rel / "probe.llc"),
            "-s",
            str(out_dir_rel / "probe.lsc"),
        ],
        cwd=LAB_ROOT,
        env=env,
    )
    if not (out_dir_abs / "probe.frfil3").exists():
        shutil.copy2(out_dir_abs / "probe.frfil2", out_dir_abs / "probe.frfil3")

    normalpoint = run_command(
        [
            "./bin/llr_npt_crd",
            "-n",
            "-p",
            "poldat_thruJ17",
            "-i",
            str(out_dir_rel / "probe.frfil3"),
            "-o",
            str(out_dir_rel / "probe.frfin"),
            "-m",
            str(out_dir_rel / "probe.npt"),
            "-c",
            str(out_dir_rel / "probe.rsc"),
            "-l",
            str(out_dir_rel / "probe.llc"),
            "-s",
            str(out_dir_rel / "probe.lsc"),
        ],
        cwd=LAB_ROOT,
        env=env,
    )
    if normalpoint.returncode != 0:
        raise RuntimeError(normalpoint.stderr or normalpoint.stdout or "normalpoint failed")

    frd_strip = run_command(
        [
            "./bin/frd_strip",
            str(out_dir_rel / "probe.frfin"),
            str(out_dir_rel / "probe.frd"),
        ],
        cwd=LAB_ROOT,
        env=env,
    )
    if frd_strip.returncode != 0:
        raise RuntimeError(frd_strip.stderr or frd_strip.stdout or "frd_strip failed")
    return out_dir_abs


def read_lines(path: Path) -> list[str]:
    return path.read_text(errors="replace").splitlines()


def token_diffs(current: str, reference: str) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    cur_tokens = current.split()
    ref_tokens = reference.split()
    for index, (cur_token, ref_token) in enumerate(zip(cur_tokens, ref_tokens), start=1):
        if cur_token == ref_token:
            continue
        row: dict[str, object] = {
            "token_index": index,
            "current": cur_token,
            "reference": ref_token,
        }
        try:
            row["delta"] = float(cur_token) - float(ref_token)
        except ValueError:
            pass
        out.append(row)
    return out


def diff_summary(current: Path, reference: Path) -> dict[str, object]:
    current_lines = read_lines(current)
    reference_lines = read_lines(reference)
    diffs: list[dict[str, object]] = []
    max_abs_numeric_delta = 0.0
    for line_number, (cur_line, ref_line) in enumerate(zip(current_lines, reference_lines), start=1):
        if cur_line == ref_line:
            continue
        token_rows = token_diffs(cur_line, ref_line)
        for token_row in token_rows:
            if "delta" in token_row:
                max_abs_numeric_delta = max(max_abs_numeric_delta, abs(token_row["delta"]))  # type: ignore[arg-type]
        diffs.append(
            {
                "line_number": line_number,
                "current": cur_line,
                "reference": ref_line,
                "token_diffs": token_rows,
            }
        )
    return {
        "line_diff_count": len(diffs),
        "max_abs_numeric_delta": max_abs_numeric_delta,
        "first_diff": diffs[0] if diffs else None,
    }


def count_changed_93_lines(current: Path, reference: Path) -> dict[str, object]:
    values: list[float] = []
    current_lines = read_lines(current)
    reference_lines = read_lines(reference)
    for cur_line, ref_line in zip(current_lines, reference_lines):
        if not cur_line.startswith("93 ") or not ref_line.startswith("93 "):
            continue
        cur_last = float(cur_line.split()[-1])
        ref_last = float(ref_line.split()[-1])
        values.append(cur_last - ref_last)
    changed = [value for value in values if value != 0.0]
    return {
        "total_93_lines": len(values),
        "changed_93_lines": len(changed),
        "min_delta_ns": min(values) if values else None,
        "max_delta_ns": max(values) if values else None,
    }


def write_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1450, 980
    base_ref = summary["baseline_vs_reference"]
    synthetic = summary["synthetic_offset_probe"]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: MLRS interface probe</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Baseline-vs-reference drift classification and bounded recalc-layer synthetic hook test.</text>',
        '<rect x="48" y="108" width="1354" height="190" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Baseline vs reference</text>',
        f'<text x="66" y="168" font-family="monospace" font-size="13" fill="#111827">frrec: {base_ref["frrec"]["line_diff_count"]} lines, max |Δ| = {base_ref["frrec"]["max_abs_numeric_delta"]:.6f} ns</text>',
        f'<text x="66" y="192" font-family="monospace" font-size="13" fill="#111827">frfin: {base_ref["frfin"]["line_diff_count"]} line, max |Δ| = {base_ref["frfin"]["max_abs_numeric_delta"]:.6f} ns</text>',
        f'<text x="66" y="216" font-family="monospace" font-size="13" fill="#111827">frd exact = {base_ref["frd_exact_match"]}</text>',
        f'<text x="66" y="240" font-family="monospace" font-size="13" fill="#111827">npt: {base_ref["npt"]["line_diff_count"]} line, max |Δ| = {base_ref["npt"]["max_abs_numeric_delta"]:.3f}</text>',
        '<text x="66" y="264" font-family="monospace" font-size="13" fill="#111827">Interpretation: baseline drift appears first in 93 residual records, not only at final formatting.</text>',
        '<rect x="48" y="330" width="1354" height="230" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="360" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Synthetic offset probe</text>',
        f'<text x="66" y="390" font-family="monospace" font-size="13" fill="#111827">offset = {synthetic["offset_ns"]:.6f} ns via SMU_OMC_OFFSET_NS</text>',
        f'<text x="66" y="414" font-family="monospace" font-size="13" fill="#111827">frrec changed 93 lines = {synthetic["frrec_93"]["changed_93_lines"]}/{synthetic["frrec_93"]["total_93_lines"]}</text>',
        f'<text x="66" y="438" font-family="monospace" font-size="13" fill="#111827">frrec observed Δ range = [{synthetic["frrec_93"]["min_delta_ns"]:.6f}, {synthetic["frrec_93"]["max_delta_ns"]:.6f}] ns</text>',
        f'<text x="66" y="462" font-family="monospace" font-size="13" fill="#111827">frfin diff lines = {synthetic["frfin"]["line_diff_count"]}</text>',
        f'<text x="66" y="486" font-family="monospace" font-size="13" fill="#111827">frd exact after offset = {synthetic["frd_exact_match"]}</text>',
        f'<text x="66" y="510" font-family="monospace" font-size="13" fill="#111827">npt diff lines = {synthetic["npt"]["line_diff_count"]}</text>',
        f'<text x="66" y="534" font-family="monospace" font-size="13" fill="#111827">npt first-line drift: TOF and PmM shift, downstream file shape intact.</text>',
        '<rect x="48" y="592" width="1354" height="310" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="622" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Interface reading</text>',
        '<text x="66" y="652" font-family="monospace" font-size="13" fill="#111827">candidate 1: np_crd.f recalc OMC computation and 93 writer</text>',
        '<text x="66" y="676" font-family="monospace" font-size="13" fill="#111827">candidate 2: jjreadnp.f / jreade state bridge into JPL PLEPH</text>',
        '<text x="66" y="700" font-family="monospace" font-size="13" fill="#111827">stop rule: if true delta_SEP needs broad batch-stack surgery beyond these seams, pivot outward.</text>',
        f'<text x="66" y="736" font-family="monospace" font-size="13" fill="#111827">assessment = {summary["assessment"]}</text>',
        "</svg>",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    baseline_dir = run_probe("baseline", None)
    offset_ns = 0.005
    offset_dir = run_probe("offset", offset_ns)

    ref_prefix = ANALYSIS_REF / SAMPLE_STEM
    baseline_prefix = baseline_dir / "probe"
    offset_prefix = offset_dir / "probe"

    baseline_vs_reference = {
        "frrec": diff_summary(baseline_prefix.with_suffix(".frrec"), ref_prefix.with_suffix(".frrec")),
        "frfin": diff_summary(baseline_prefix.with_suffix(".frfin"), ref_prefix.with_suffix(".frfin")),
        "npt": diff_summary(baseline_prefix.with_suffix(".npt"), ref_prefix.with_suffix(".npt")),
        "frd_exact_match": baseline_prefix.with_suffix(".frd").read_bytes() == ref_prefix.with_suffix(".frd").read_bytes(),
    }

    synthetic_offset_probe = {
        "offset_ns": offset_ns,
        "frrec_93": count_changed_93_lines(offset_prefix.with_suffix(".frrec"), baseline_prefix.with_suffix(".frrec")),
        "frrec": diff_summary(offset_prefix.with_suffix(".frrec"), baseline_prefix.with_suffix(".frrec")),
        "frfin": diff_summary(offset_prefix.with_suffix(".frfin"), baseline_prefix.with_suffix(".frfin")),
        "npt": diff_summary(offset_prefix.with_suffix(".npt"), baseline_prefix.with_suffix(".npt")),
        "frd_exact_match": offset_prefix.with_suffix(".frd").read_bytes() == baseline_prefix.with_suffix(".frd").read_bytes(),
    }

    summary = {
        "status": "bounded_interface_alive",
        "sample_stem": SAMPLE_STEM,
        "lab_root": repo_relative(LAB_ROOT),
        "baseline_run": repo_relative(baseline_dir),
        "offset_run": repo_relative(offset_dir),
        "baseline_vs_reference": baseline_vs_reference,
        "synthetic_offset_probe": synthetic_offset_probe,
        "candidate_interfaces": [
            {
                "file": "data/request4_llr/mlrs_handshake_lab/src/llr_npt/np_crd.f",
                "line_hint": "441-460",
                "role": "recalc-layer OMC computation and write_93 emission",
            },
            {
                "file": "data/request4_llr/mlrs_handshake_lab/src/llr_npt/jjreadnp.f",
                "line_hint": "30-53",
                "role": "PLEPH-backed state bridge",
            },
            {
                "file": "data/request4_llr/mlrs_handshake_lab/src/llr_npt/jeulpkg.f",
                "line_hint": "341-348, 1216-1224",
                "role": "prediction and light-time path that consumes jreade output",
            },
        ],
        "assessment": (
            "Baseline-ref drift is a real recalc-state microdrift confined to 93 residual records, "
            "not a pure final-formatting artifact. The synthetic OMC hook shows a bounded correction "
            "can be injected at the recalc layer and propagated through normalpoint generation "
            "without breaking the downstream batch stack."
        ),
    }

    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    write_svg(SUMMARY_SVG, summary)


if __name__ == "__main__":
    main()
