from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SUMMARY_JSON = ROOT / "request4_llr_external_estimator_summary.json"
SUMMARY_SVG = ROOT / "request4_llr_external_estimator_summary.svg"
CANDIDATES_TSV = ROOT / "request4_llr_external_estimator_candidates.tsv"


@dataclass(frozen=True)
class Candidate:
    name: str
    owner: str
    kind: str
    llr_fit_engine: str
    public_code_status: str
    public_data_status: str
    integration_effort: str
    recommendation_tier: str
    rationale: str
    primary_sources: tuple[str, ...]


def write_tsv(path: Path, rows: list[Candidate]) -> None:
    fieldnames = list(asdict(rows[0]).keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            data = asdict(row)
            data["primary_sources"] = " | ".join(row.primary_sources)
            writer.writerow(data)


def write_svg(path: Path, summary: dict[str, object], rows: list[Candidate]) -> None:
    width, height = 1500, 980
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: external LLR estimator scout</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Post-MLRS hand-off ranking for the final weak-field posterior branch.</text>',
        '<rect x="48" y="108" width="1404" height="150" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Decision</text>',
        f'<text x="66" y="170" font-family="monospace" font-size="13" fill="#111827">primary hand-off target = {summary["primary_hand_off_target"]}</text>',
        f'<text x="66" y="194" font-family="monospace" font-size="13" fill="#111827">secondary path = {summary["secondary_hand_off_path"]}</text>',
        f'<text x="66" y="218" font-family="monospace" font-size="13" fill="#111827">request4 verdict = {summary["request4_verdict"]}</text>',
        '<rect x="48" y="292" width="1404" height="612" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="322" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Candidates</text>',
    ]
    y = 354
    for row in rows:
        lines.extend(
            [
                f'<text x="66" y="{y}" font-family="sans-serif" font-size="15" font-weight="700" fill="#111827">{row.recommendation_tier} · {row.name}</text>',
                f'<text x="66" y="{y + 22}" font-family="monospace" font-size="12" fill="#111827">{row.owner} | {row.kind} | llr_fit={row.llr_fit_engine} | public_code={row.public_code_status} | effort={row.integration_effort}</text>',
                f'<text x="66" y="{y + 44}" font-family="sans-serif" font-size="12" fill="#334155">{row.rationale}</text>',
            ]
        )
        y += 78
    lines.extend(
        [
            '<rect x="48" y="930" width="1404" height="18" fill="none"/>',
            '</svg>',
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    rows = [
        Candidate(
            name="PEP (Planetary Ephemeris Program)",
            owner="CfA / ASCL / PEP paper",
            kind="full astrometric / ephemeris fit engine with LLR support",
            llr_fit_engine="yes",
            public_code_status="public-code candidate; paper says source/docs/tests are publicly available via GitLab",
            public_data_status="works with external ephemerides and heterogeneous astrometric data",
            integration_effort="medium-high",
            recommendation_tier="A",
            rationale=(
                "Best external hand-off target found: it is described as an open-source general-purpose astrometric analysis program, "
                "its literature explicitly states lunar laser ranging support down to mm-class one-way precision, and ILRS literature shows it has already been used on real LLR data."
            ),
            primary_sources=(
                "https://arxiv.org/abs/2103.16745",
                "https://ascl.net/code/search/orbit",
                "https://ilrs.cddis.eosdis.nasa.gov/lw19/docs/2014/Papers/3148_Martini_paper.pdf",
            ),
        ),
        Candidate(
            name="JPL Lunar Analysis Center / DE pipeline",
            owner="JPL / ILRS",
            kind="mature internal LLR analysis center workflow",
            llr_fit_engine="yes",
            public_code_status="ephemerides and reader toolkits public; full fit engine not identified as public",
            public_data_status="public DE ephemerides and SPICE readers available",
            integration_effort="high or collaboration-only",
            recommendation_tier="B",
            rationale=(
                "Clearly mature and high-performing for LLR, but the public-facing assets are the DE ephemerides and access toolkits rather than a released end-to-end estimator. "
                "This is a credible collaboration / external-analysis route, not a drop-in public code hand-off."
            ),
            primary_sources=(
                "https://ssd.jpl.nasa.gov/planets/eph_export.html",
                "https://ilrs.gsfc.nasa.gov/docs/2020/ilrsreport_2016_section7.pdf",
            ),
        ),
        Candidate(
            name="POLAC / ELPN",
            owner="Paris Observatory Lunar Analysis Center",
            kind="mature lunar analysis-center ephemeris and fit workflow",
            llr_fit_engine="yes",
            public_code_status="public analysis-center capability; no public code path identified in bounded scout",
            public_data_status="public predictions and analysis-center outputs described",
            integration_effort="high or collaboration-only",
            recommendation_tier="B",
            rationale=(
                "Strong science maturity and explicit use for LLR fundamental-physics fits, including variational equations and SME-style extensions, "
                "but the bounded scout did not identify a public code release suitable as a local estimator drop-in."
            ),
            primary_sources=(
                "https://ilrs.gsfc.nasa.gov/docs/2020/ilrsreport_2016_section7.pdf",
                "https://ilrs.gsfc.nasa.gov/lw22/papers/S10/S10-06_Bourgoin_Paper.pdf",
            ),
        ),
        Candidate(
            name="INPOP + CALCEPH",
            owner="IMCCE",
            kind="public ephemerides + open-source ephemeris access layer",
            llr_fit_engine="not in public code found here",
            public_code_status="reader libraries public; fit engine not identified as public",
            public_data_status="ephemeris files public; Moon residual computation offered to users",
            integration_effort="medium for ephemeris side, high for full posterior",
            recommendation_tier="C",
            rationale=(
                "Useful if the project needs a public ephemeris / residual-computation substrate, but the public assets found are the ephemeris products and CALCEPH readers, not a released full LLR estimator."
            ),
            primary_sources=(
                "https://www.imcce.fr/inpop/",
                "https://calceph.imcce.fr/docs/4.0.1/html/c/calceph.intro.html",
                "https://www.imcce.fr/inpop/calceph",
            ),
        ),
        Candidate(
            name="EPM / IAA RAS",
            owner="IAA RAS",
            kind="mature ephemeris and web residual service",
            llr_fit_engine="yes internally",
            public_code_status="officially not distributed",
            public_data_status="ephemerides and web O-C services public",
            integration_effort="high or service-only",
            recommendation_tier="C",
            rationale=(
                "Scientifically mature and explicitly active in LLR residual analysis, but the official site states that the scientific software engine is not distributed, "
                "so this is not a local code hand-off candidate."
            ),
            primary_sources=(
                "https://iaaras.ru/en/dept/ephemeris/online/references/",
                "https://ilrs.gsfc.nasa.gov/docs/2020/ilrsreport_2016_section7.pdf",
            ),
        ),
    ]

    write_tsv(CANDIDATES_TSV, rows)

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "request4_verdict": "freeze MLRS as bounded observation-operator lab; hand off final weak-field posterior externally",
        "primary_hand_off_target": "PEP (if code acquisition / build is feasible)",
        "secondary_hand_off_path": "collaboration or external-analysis-center route via JPL / POLAC-ELPN",
        "public_ephemeris_support_path": "INPOP + CALCEPH or JPL DE + SPICE are support layers, not full fit-engine replacements",
        "candidates": [asdict(row) for row in rows],
        "interpretation": [
            "PEP is the strongest public-code candidate found for an actual LLR-capable estimator hand-off.",
            "JPL and POLAC are scientifically mature but were not identified here as public drop-in estimator releases.",
            "INPOP/CALCEPH and JPL DE/SPICE are valuable ephemeris/residual-access layers but not sufficient by themselves to replace a full weak-field estimator fit engine.",
            "IAA RAS EPM is mature but the official site states the engine is not distributed.",
        ],
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    write_svg(SUMMARY_SVG, summary, rows)


if __name__ == "__main__":
    main()
