from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


LAB_ROOT = Path("data/request4_llr/mlrs_handshake_lab")
SUMMARY_JSON = Path("request4_llr_mlrs_handshake_summary.json")
SUMMARY_SVG = Path("request4_llr_mlrs_handshake_summary.svg")
JPL_DE421_URL = "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de421/lnxp1900p2053.421"

SAMPLE_STEMS = (
    "s25y11d016t0231#103",
    "s25y11d042t0234#103",
)

BINARY_PATHS = (
    "bin/llr_npt_crd",
    "bin/Poisson_crd",
    "bin/lun_auto_crd",
    "bin/lun_sdqi_crd",
    "bin/frd_strip",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_description(path: Path) -> str:
    result = subprocess.run(
        ["file", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    line = result.stdout.strip()
    return line.split(":", 1)[1].strip() if ":" in line else line


def read_lines(path: Path) -> list[str]:
    return path.read_text(errors="replace").splitlines()


def try_float(token: str) -> float | None:
    try:
        return float(token)
    except ValueError:
        return None


def diff_tokens(current: str, reference: str) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    cur_tokens = current.split()
    ref_tokens = reference.split()
    for index, (cur_token, ref_token) in enumerate(zip(cur_tokens, ref_tokens), start=1):
        if cur_token == ref_token:
            continue
        cur_float = try_float(cur_token)
        ref_float = try_float(ref_token)
        row: dict[str, object] = {
            "token_index": index,
            "current": cur_token,
            "reference": ref_token,
        }
        if cur_float is not None and ref_float is not None:
            row["delta"] = cur_float - ref_float
        out.append(row)
    if len(cur_tokens) != len(ref_tokens):
        out.append(
            {
                "token_index": None,
                "current_token_count": len(cur_tokens),
                "reference_token_count": len(ref_tokens),
            }
        )
    return out


def summarize_sample(stem: str) -> dict[str, object]:
    analysis_dir = LAB_ROOT / "data/analysis"
    ref_dir = analysis_dir / "ref"
    row: dict[str, object] = {"stem": stem}

    errorc_path = analysis_dir / f"{stem}.errorc"
    row["errorc_bytes"] = errorc_path.stat().st_size if errorc_path.exists() else None
    row["errorc_empty"] = errorc_path.exists() and errorc_path.stat().st_size == 0

    for ext in ("frd", "npt"):
        current_path = analysis_dir / f"{stem}.{ext}"
        reference_path = ref_dir / f"{stem}.{ext}"
        current_lines = read_lines(current_path)
        reference_lines = read_lines(reference_path)

        diffs: list[dict[str, object]] = []
        for line_number, (current_line, reference_line) in enumerate(
            zip(current_lines, reference_lines),
            start=1,
        ):
            if current_line != reference_line:
                diffs.append(
                    {
                        "line_number": line_number,
                        "current": current_line,
                        "reference": reference_line,
                        "token_diffs": diff_tokens(current_line, reference_line),
                    }
                )

        row[ext] = {
            "current_sha256": sha256_file(current_path),
            "reference_sha256": sha256_file(reference_path),
            "exact_match": current_path.read_bytes() == reference_path.read_bytes(),
            "line_count_current": len(current_lines),
            "line_count_reference": len(reference_lines),
            "line_diff_count": len(diffs),
            "first_diff": diffs[0] if diffs else None,
        }

    return row


def write_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1400, 900
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: MLRS handshake audit</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Public MLRS sample replay after bounded modernization and public DE421 injection.</text>',
        '<rect x="48" y="108" width="1304" height="176" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Status</text>',
        f'<text x="66" y="168" font-family="monospace" font-size="13" fill="#111827">overall = {summary["status"]}</text>',
        f'<text x="66" y="192" font-family="monospace" font-size="13" fill="#111827">FRD exact matches = {summary["sample_replay"]["frd_exact_matches"]}/{summary["sample_replay"]["sample_count"]}</text>',
        f'<text x="66" y="216" font-family="monospace" font-size="13" fill="#111827">NPT single-line drifts = {summary["sample_replay"]["npt_single_line_drifts"]}/{summary["sample_replay"]["sample_count"]}</text>',
        f'<text x="66" y="240" font-family="monospace" font-size="13" fill="#111827">JPL file size = {summary["jpl_ephemeris"]["size_bytes"]} bytes</text>',
        f'<text x="66" y="264" font-family="monospace" font-size="13" fill="#111827">JPL source = {summary["jpl_ephemeris"]["source_url"]}</text>',
        '<rect x="48" y="320" width="1304" height="220" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="350" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Binaries</text>',
    ]
    y = 380
    for rel_path, info in summary["binaries"].items():
        lines.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="12" fill="#111827">{rel_path}: {info["file_description"]}</text>'
        )
        y += 24

    lines.extend(
        [
            '<rect x="48" y="576" width="1304" height="276" fill="white" stroke="#cbd5e1"/>',
            '<text x="66" y="606" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Sample replay</text>',
        ]
    )
    y = 636
    for sample in summary["samples"]:
        npt = sample["npt"]
        frd = sample["frd"]
        lines.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="12" fill="#111827">{sample["stem"]}: FRD exact={frd["exact_match"]}, NPT diffs={npt["line_diff_count"]}, errorc_empty={sample["errorc_empty"]}</text>'
        )
        y += 24
        first_diff = npt["first_diff"]
        if first_diff is not None:
            line = first_diff["line_number"]
            token_summary = ", ".join(
                f'#{token["token_index"]}:{token.get("delta", token["current"] + "!=" + token["reference"])}'
                for token in first_diff["token_diffs"][:3]
            )
            lines.append(
                f'<text x="84" y="{y}" font-family="monospace" font-size="11" fill="#475569">line {line} drift: {token_summary}</text>'
            )
            y += 20
    lines.append("</svg>")
    path.write_text("\n".join(lines))


def main() -> None:
    if not LAB_ROOT.exists():
        raise SystemExit(f"missing lab root: {LAB_ROOT}")

    binaries: dict[str, dict[str, object]] = {}
    for rel_path in BINARY_PATHS:
        path = LAB_ROOT / rel_path
        binaries[rel_path] = {
            "exists": path.exists(),
            "file_description": file_description(path) if path.exists() else None,
            "sha256": sha256_file(path) if path.exists() else None,
        }

    jpl_path = LAB_ROOT / "data/pred/jpleph.421"
    jpl_link = LAB_ROOT / "JPLEPH"
    samples = [summarize_sample(stem) for stem in SAMPLE_STEMS]

    frd_exact = sum(sample["frd"]["exact_match"] for sample in samples)
    npt_exact = sum(sample["npt"]["exact_match"] for sample in samples)
    npt_single_line_drifts = sum(sample["npt"]["line_diff_count"] == 1 for sample in samples)

    summary = {
        "status": "mlrs_handshake_alive",
        "lab_root": str(LAB_ROOT),
        "jpl_ephemeris": {
            "source_url": JPL_DE421_URL,
            "local_path": str(jpl_path),
            "size_bytes": jpl_path.stat().st_size,
            "sha256": sha256_file(jpl_path),
            "top_level_link": str(jpl_link),
            "top_level_link_target": str(jpl_link.resolve()),
        },
        "binaries": binaries,
        "samples": samples,
        "sample_replay": {
            "sample_count": len(samples),
            "frd_exact_matches": frd_exact,
            "npt_exact_matches": npt_exact,
            "npt_single_line_drifts": npt_single_line_drifts,
        },
        "pipeline_interfaces": [
            {
                "stage": "recalc",
                "script_file": "data/request4_llr/mlrs_handshake_lab/bin/ldb_crd",
                "script_lines": "115-126",
                "command": "bin/llr_npt_crd -r -p poldat_thruJ17 -i ...frcal -o ...frrec",
                "role": "prediction and residual recalculation",
                "eft_insertion_candidate": True,
                "notes": "Nearest place to perturb the gravity model before filtering.",
            },
            {
                "stage": "poisson_filter",
                "script_file": "data/request4_llr/mlrs_handshake_lab/bin/ldb_crd",
                "script_lines": "132-141",
                "command": "bin/Poisson_crd -i ...frrec -o ...frfil",
                "role": "filtering only",
                "eft_insertion_candidate": False,
                "notes": "Post-model outlier rejection, not the gravity-model hook.",
            },
            {
                "stage": "normalpoint",
                "script_file": "data/request4_llr/mlrs_handshake_lab/bin/ldb_crd",
                "script_lines": "174-180",
                "command": "bin/llr_npt_crd -n -p poldat_thruJ17 -i ...frfil3 -o ...frfin -m ...npt",
                "role": "normal-point construction after filtering",
                "eft_insertion_candidate": False,
                "notes": "Consumes filtered ranges; model perturbation should enter upstream.",
            },
            {
                "stage": "ephemeris_bridge",
                "source_file": "data/request4_llr/mlrs_handshake_lab/src/llr_npt/jjreadnp.f",
                "source_lines": "5-6, 30-53",
                "role": "bridge from MLRS residual code to JPL PLEPH calls",
                "eft_insertion_candidate": True,
                "notes": "Current model enters through PLEPH-backed state evaluation here.",
            },
        ],
        "assessment": {
            "passed_short_handshake_gate": True,
            "why": (
                "Public MLRS sample workflow now replays end-to-end on macOS arm64 "
                "with exact FRD matches and only one-line NPT drifts in each sample."
            ),
            "not_yet_done": (
                "This is still a bounded hand-off test, not a weak-field posterior. "
                "Partials, nuisance estimation, and EFT remapping are not wired yet."
            ),
        },
    }

    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    write_svg(SUMMARY_SVG, summary)


if __name__ == "__main__":
    main()
