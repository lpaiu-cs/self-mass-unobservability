from __future__ import annotations

import csv
import json
import shutil
import tarfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT / "data" / "request5_j0337" / "nutimo_public_release_2026-04-23"
Y2020 = DATA_ROOT / "2020"
Y2025 = DATA_ROOT / "2025"

SUMMARY_JSON = ROOT / "request5_j0337_phaseB_public_inputs_summary.json"
SUMMARY_SVG = ROOT / "request5_j0337_phaseB_public_inputs_summary.svg"
MANIFEST_TSV = ROOT / "request5_j0337_phaseB_public_inputs_manifest.tsv"


@dataclass(frozen=True)
class Asset:
    release: str
    doi: str
    asset: str
    local_status: str
    local_path: str
    size_bytes: int
    role: str


def repo_relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def tar_members(path: Path) -> list[str]:
    with tarfile.open(path, "r:bz2") as handle:
        return sorted(member.name for member in handle.getmembers() if member.isfile())


def write_tsv(path: Path, rows: list[Asset]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1500, 980
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 5: J0337 Phase B public inputs</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Bounded scout of public TOA, Nutimo code, and immediate local dependency gaps.</text>',
        '<rect x="48" y="108" width="1404" height="190" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Decision</text>',
        f'<text x="66" y="170" font-family="monospace" font-size="13" fill="#111827">verdict = {summary["verdict"]}</text>',
        f'<text x="66" y="194" font-family="monospace" font-size="13" fill="#111827">2020 release mirrored = {summary["release_2020_mirrored"]} | 2025 release mirrored = {summary["release_2025_mirrored"]}</text>',
        f'<text x="66" y="218" font-family="monospace" font-size="13" fill="#111827">tim files found = {summary["tim_files_found"]} | par files found = {summary["par_files_found"]}</text>',
        f'<text x="66" y="242" font-family="monospace" font-size="13" fill="#111827">nutimo source files listed = {summary["nutimo_source_files_found"]}</text>',
        f'<text x="66" y="266" font-family="monospace" font-size="13" fill="#111827">missing local deps = {summary["missing_local_dependencies"]}</text>',
        '<rect x="48" y="332" width="1404" height="594" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="362" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Key files</text>',
    ]
    y = 396
    for label, values in (
        ("2020 tim", summary["release_2020_tim_files"]),
        ("2020 par", summary["release_2020_par_files"]),
        ("2025 tim", summary["release_2025_tim_files"]),
        ("2025 code", summary["release_2025_code_highlights"]),
    ):
        lines.append(
            f'<text x="66" y="{y}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">{label}</text>'
        )
        y += 20
        for item in values[:8]:
            lines.append(
                f'<text x="84" y="{y}" font-family="monospace" font-size="12" fill="#334155">{item}</text>'
            )
            y += 18
        y += 18
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    assets = [
        Asset(
            release="2020",
            doi="10.5281/zenodo.3778978",
            asset="README",
            local_status="mirrored",
            local_path=repo_relative(Y2020 / "README"),
            size_bytes=(Y2020 / "README").stat().st_size,
            role="release readme",
        ),
        Asset(
            release="2020",
            doi="10.5281/zenodo.3778978",
            asset="nutimo.tar.bz2",
            local_status="mirrored",
            local_path=repo_relative(Y2020 / "nutimo.tar.bz2"),
            size_bytes=(Y2020 / "nutimo.tar.bz2").stat().st_size,
            role="timing-model source code",
        ),
        Asset(
            release="2020",
            doi="10.5281/zenodo.3778978",
            asset="Data_and_Results.tar.bz2",
            local_status="mirrored",
            local_path=repo_relative(Y2020 / "Data_and_Results.tar.bz2"),
            size_bytes=(Y2020 / "Data_and_Results.tar.bz2").stat().st_size,
            role="TOA, parfile, MCMC results bundle",
        ),
        Asset(
            release="2025",
            doi="10.5281/zenodo.13899771",
            asset="Readme.md",
            local_status="mirrored",
            local_path=repo_relative(Y2025 / "Readme.md"),
            size_bytes=(Y2025 / "Readme.md").stat().st_size,
            role="release readme",
        ),
        Asset(
            release="2025",
            doi="10.5281/zenodo.13899771",
            asset="nutimo.tar.bz2",
            local_status="mirrored",
            local_path=repo_relative(Y2025 / "nutimo.tar.bz2"),
            size_bytes=(Y2025 / "nutimo.tar.bz2").stat().st_size,
            role="timing-model source code",
        ),
        Asset(
            release="2025",
            doi="10.5281/zenodo.13899771",
            asset="Data.tar.bz2",
            local_status="mirrored",
            local_path=repo_relative(Y2025 / "Data.tar.bz2"),
            size_bytes=(Y2025 / "Data.tar.bz2").stat().st_size,
            role="TOA bundle used in paper",
        ),
        Asset(
            release="2025",
            doi="10.5281/zenodo.13899771",
            asset="Analysis.tar.bz2",
            local_status="remote-only",
            local_path="",
            size_bytes=247273012,
            role="full analysis products; left remote in bounded scout",
        ),
    ]
    write_tsv(MANIFEST_TSV, assets)

    members_2020 = tar_members(Y2020 / "Data_and_Results.tar.bz2")
    members_2025_code = tar_members(Y2025 / "nutimo.tar.bz2")
    members_2025_data = tar_members(Y2025 / "Data.tar.bz2")

    tim_2020 = [name for name in members_2020 if name.endswith(".tim") or ".tim-" in name]
    par_2020 = [name for name in members_2020 if ".par" in name]
    tim_2025 = [name for name in members_2025_data if name.endswith(".tim")]
    par_2025 = [name for name in members_2025_data if ".par" in name]

    source_2025 = [name for name in members_2025_code if name.startswith("nutimo/src/")]
    code_highlights = [
        name
        for name in members_2025_code
        if name.endswith(("install_script.sh", "README.md", "makefile-original", ".pyx", ".cpp"))
    ]

    local_commands = {
        name: shutil.which(name)
        for name in ("tempo2", "root-config", "g++", "gcc", "python3", "cython")
    }
    missing_local_dependencies = [
        name
        for name in ("tempo2", "root-config", "cython")
        if not local_commands.get(name)
    ]

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": (
            "Phase B public inputs are open and locally mirrored: public TOA files, parfiles, "
            "and Nutimo source are available from Zenodo. The next gate is build/runtime closure, "
            "not data access."
        ),
        "release_2020_mirrored": True,
        "release_2025_mirrored": True,
        "tim_files_found": len(tim_2020) + len(tim_2025),
        "par_files_found": len(par_2020) + len(par_2025),
        "nutimo_source_files_found": len(source_2025),
        "release_2020_doi": "10.5281/zenodo.3778978",
        "release_2025_doi": "10.5281/zenodo.13899771",
        "release_2020_tim_files": tim_2020,
        "release_2020_par_files": par_2020,
        "release_2025_tim_files": tim_2025,
        "release_2025_par_files": par_2025,
        "release_2025_code_highlights": code_highlights[:12],
        "release_2025_dependency_statement": (
            "Nutimo README names Boost, Tempo2, Minuit, and Acor; install_script.sh also expects "
            "manual third_party paths and compiles C++ plus Cython extensions."
        ),
        "local_command_paths": local_commands,
        "missing_local_dependencies": missing_local_dependencies,
        "next_bounded_action": (
            "perform a Nutimo build/runtime feasibility probe in the current workspace without "
            "turning Request 5 into an open-ended dependency provisioning project"
        ),
        "assets": [asdict(asset) for asset in assets],
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    write_svg(SUMMARY_SVG, summary)


if __name__ == "__main__":
    main()
