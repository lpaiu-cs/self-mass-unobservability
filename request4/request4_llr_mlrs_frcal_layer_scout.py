from __future__ import annotations

import json
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent
MANUAL_DOC = REPO_ROOT / "data/request4_llr/mlrs_handshake_lab/MLRS Lunar Prediction and Normal Point Manual-v1.0.doc"
LDB_CRD = REPO_ROOT / "data/request4_llr/mlrs_handshake_lab/bin/ldb_crd"
PROMOTION_SUMMARY = ROOT / "request4_llr_mlrs_promotion_audit_summary.json"
REFERENCE_FRREC = REPO_ROOT / "data/request4_llr/mlrs_handshake_lab/data/analysis/s25y11d016t0231#103.frrec"
PROMOTED_FRFIL = {
    "apollo15_202503_wetl_promo": REPO_ROOT / "data/request4_llr/mlrs_handshake_lab/data/analysis/apollo15_202503_wetl_promo.frfil",
    "apollo15_202511_apol_promo": REPO_ROOT / "data/request4_llr/mlrs_handshake_lab/data/analysis/apollo15_202511_apol_promo.frfil",
    "apollo15_202502_grsm_promo": REPO_ROOT / "data/request4_llr/mlrs_handshake_lab/data/analysis/apollo15_202502_grsm_promo.frfil",
}
SUMMARY_JSON = ROOT / "request4_llr_mlrs_frcal_layer_summary.json"
SUMMARY_SVG = ROOT / "request4_llr_mlrs_frcal_layer_summary.svg"


@dataclass(frozen=True)
class LineHit:
    line_number: int
    text: str


def convert_manual_to_text() -> list[str]:
    result = subprocess.run(
        ["textutil", "-convert", "txt", "-stdout", str(MANUAL_DOC)],
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout.splitlines()


def first_match(lines: list[str], needle: str) -> LineHit:
    for index, line in enumerate(lines, start=1):
        if needle.lower() in line.lower():
            return LineHit(line_number=index, text=line.strip())
    raise ValueError(f"needle not found: {needle}")


def first_match_file(path: Path, needle: str) -> LineHit:
    for index, line in enumerate(path.read_text().splitlines(), start=1):
        if needle.lower() in line.lower():
            return LineHit(line_number=index, text=line.strip())
    raise ValueError(f"needle not found in {path}: {needle}")


def count_tags(path: Path, tags: tuple[str, ...]) -> dict[str, int]:
    counts = {tag: 0 for tag in tags}
    with path.open() as handle:
        for raw_line in handle:
            tag = raw_line[:2]
            if tag in counts:
                counts[tag] += 1
    return counts


def write_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1440, 960
    evidence = summary["evidence"]
    manual = evidence["local_manual"]
    bundle = evidence["local_bundle"]
    public = evidence["public_promotion"]
    ref_counts = summary["reference_counts"]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: MLRS frcal/93 layer scout</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Bounded search for a public or semi-public raw-CRD to residual-bearing frcal promotion layer.</text>',
        '<rect x="48" y="108" width="1344" height="136" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Verdict</text>',
        f'<text x="66" y="170" font-family="monospace" font-size="13" fill="#111827">missing layer status = {summary["missing_layer_status"]}</text>',
        f'<text x="66" y="194" font-family="monospace" font-size="13" fill="#111827">stop rule triggered = {summary["stop_rule_triggered"]}</text>',
        f'<text x="66" y="218" font-family="monospace" font-size="13" fill="#111827">request4 status = {summary["request4_status"]}</text>',
        '<rect x="48" y="276" width="1344" height="238" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="306" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Evidence chain</text>',
        f'<text x="66" y="338" font-family="monospace" font-size="13" fill="#111827">manual says Poisson input = {manual["poisson_input"]["text"]}</text>',
        f'<text x="66" y="362" font-family="monospace" font-size="13" fill="#111827">manual says bundled tests use = {manual["frcal_input"]["text"]}</text>',
        f'<text x="66" y="386" font-family="monospace" font-size="13" fill="#111827">bundle script says = {bundle["calibration_removed"]["text"]}</text>',
        f'<text x="66" y="410" font-family="monospace" font-size="13" fill="#111827">bundle script input line = {bundle["frcal_input"]["text"]}</text>',
        f'<text x="66" y="434" font-family="monospace" font-size="13" fill="#111827">public promotion recalc pass = {public["recalc_passed_cases"]} / {public["total_cases"]}</text>',
        f'<text x="66" y="458" font-family="monospace" font-size="13" fill="#111827">public promotion poisson pass = {public["poisson_passed_cases"]} / {public["total_cases"]}</text>',
        f'<text x="66" y="482" font-family="monospace" font-size="13" fill="#111827">public promotion 93-bearing cases = {public["cases_with_93_stream"]}</text>',
        '<rect x="48" y="546" width="1344" height="318" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="576" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Reference vs promoted record structure</text>',
        f'<text x="66" y="608" font-family="monospace" font-size="13" fill="#111827">reference frrec counts = 10:{ref_counts["10"]} 12:{ref_counts["12"]} 30:{ref_counts["30"]} 92:{ref_counts["92"]} 93:{ref_counts["93"]}</text>',
    ]
    y = 638
    for label, counts in summary["promoted_counts"].items():
        lines.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="13" fill="#111827">{label} frfil counts = 10:{counts["10"]} 12:{counts["12"]} 20:{counts["20"]} 30:{counts["30"]} 40:{counts["40"]} 93:{counts["93"]}</text>'
        )
        y += 24
    y += 12
    lines.extend(
        [
            f'<text x="66" y="{y}" font-family="sans-serif" font-size="15" font-weight="700" fill="#111827">Interpretation</text>',
            f'<text x="66" y="{y + 28}" font-family="sans-serif" font-size="13" fill="#334155">Official CRD docs describe CRD as the native format before station-side calibration/filtering,</text>',
            f'<text x="66" y="{y + 50}" font-family="sans-serif" font-size="13" fill="#334155">and note that 9x station-defined records are normally stripped before transmittal.</text>',
            f'<text x="66" y="{y + 72}" font-family="sans-serif" font-size="13" fill="#334155">The public MLRS bundle expects fully calibrated .frcal input but also says MLRS-specific calibration code was removed.</text>',
            f'<text x="66" y="{y + 94}" font-family="sans-serif" font-size="13" fill="#334155">So the missing residual-bearing frcal/93 layer is not publicly exposed in the current bounded MLRS path.</text>',
        ]
    )
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    manual_lines = convert_manual_to_text()
    local_manual = {
        "poisson_input": asdict(first_match(manual_lines, "Standard, range calibrated, CRD full rate data file")),
        "frcal_input": asdict(first_match(manual_lines, "These are fully calibrated full rate files.")),
        "writes_93": asdict(first_match(manual_lines, "writes them to the CRD file in special MLRS-only 93 records")),
    }

    local_bundle = {
        "calibration_removed": asdict(first_match_file(LDB_CRD, "MLRS-specific calibration code removed")),
        "frcal_input": asdict(first_match_file(LDB_CRD, '".frcal"')),
    }

    promotion = json.loads(PROMOTION_SUMMARY.read_text())
    promotion_results = promotion["promotion_results"]
    total_cases = len(promotion_results)
    recalc_passed_cases = sum(1 for row in promotion_results if row["recalc_passed"])
    poisson_passed_cases = sum(1 for row in promotion_results if row["poisson_passed"])
    cases_with_93_stream = [
        row["filename"] for row in promotion_results if int(row["frfil_has_93"]) > 0
    ]
    invalid_payload_cases = [
        row["filename"] for row in promotion_results if row["error_signature"] == "invalid_fullrate_payload"
    ]
    public_promotion = {
        "total_cases": total_cases,
        "recalc_passed_cases": recalc_passed_cases,
        "poisson_passed_cases": poisson_passed_cases,
        "cases_with_93_stream": cases_with_93_stream,
        "invalid_payload_cases": invalid_payload_cases,
    }

    reference_counts = count_tags(REFERENCE_FRREC, ("10", "12", "30", "92", "93"))
    promoted_counts = {
        label: count_tags(path, ("10", "12", "20", "30", "40", "93"))
        for label, path in PROMOTED_FRFIL.items()
    }

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "official_public_sources": {
            "software_index": "https://ilrs.gsfc.nasa.gov/technology/software/index.html",
            "crd_overview": "https://ilrs.gsfc.nasa.gov/data_and_products/formats/crd.html",
            "crd_v1_01_manual": "https://ilrs.gsfc.nasa.gov/docs/crd_v1.01.pdf",
            "cddis_archive_access": "https://www.earthdata.nasa.gov/centers/cddis-daac/archive-access",
        },
        "evidence": {
            "local_manual": local_manual,
            "local_bundle": local_bundle,
            "public_promotion": public_promotion,
        },
        "reference_counts": reference_counts,
        "promoted_counts": promoted_counts,
        "missing_layer_status": "not_publicly_exposed_in_bounded_search",
        "stop_rule_triggered": True,
        "request4_status": "hand_off_final_weak_field_posterior_to_external_estimator",
        "verdict": (
            "No public or semi-public raw-CRD to residual-bearing frcal/93 promotion layer "
            "was found in the bounded Request 4 MLRS search."
        ),
        "interpretation": [
            "Official CRD documentation describes MLRS as converting acquisition data to CRD and then doing station-side calibration, filtering, and normal pointing there.",
            "The CRD format manual explicitly says 91/92/93 user-defined station records are normally stripped before transmittal.",
            "The public MLRS bundle expects fully calibrated .frcal input and its ldb_crd script says MLRS-specific calibration code was removed.",
            "Public monthly CRD promotion reaches recalc but still produces no 93-bearing residual stream in the promoted .frfil cases.",
        ],
    }

    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    write_svg(SUMMARY_SVG, summary)


if __name__ == "__main__":
    main()
