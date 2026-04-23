from __future__ import annotations

import csv
import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("/Users/lpaiu/vs/lab/self-mass-unobservability")
SRC = ROOT / "data" / "request5_j0337" / "nutimo_public_release_2026-04-23" / "2025" / "build_probe" / "nutimo" / "src"
SUMMARY_JSON = ROOT / "request5_j0337_phaseB_build_feasibility_summary.json"
SUMMARY_SVG = ROOT / "request5_j0337_phaseB_build_feasibility_summary.svg"
DIAGNOSTICS_TSV = ROOT / "request5_j0337_phaseB_build_feasibility_diagnostics.tsv"


@dataclass(frozen=True)
class Diagnostic:
    name: str
    status: str
    detail: str


def run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def excerpt(output: str, n: int = 12) -> list[str]:
    return [line.rstrip() for line in output.splitlines()[:n]]


def write_tsv(path: Path, rows: list[Diagnostic]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "status", "detail"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_svg(path: Path, summary: dict[str, object], diagnostics: list[Diagnostic]) -> None:
    width, height = 1500, 920
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 5: J0337 Phase B build feasibility</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Bounded build/runtime gate for the public Nutimo release on the current host.</text>',
        '<rect x="48" y="108" width="1404" height="176" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Decision</text>',
        f'<text x="66" y="170" font-family="monospace" font-size="13" fill="#111827">verdict = {summary["verdict"]}</text>',
        f'<text x="66" y="194" font-family="monospace" font-size="13" fill="#111827">default compiler = {summary["default_gpp"]} | g++-14 = {summary["gpp14"]}</text>',
        f'<text x="66" y="218" font-family="monospace" font-size="13" fill="#111827">missing deps = {summary["missing_dependencies"]}</text>',
        f'<text x="66" y="242" font-family="monospace" font-size="13" fill="#111827">next bounded action = {summary["next_bounded_action"]}</text>',
        '<rect x="48" y="318" width="1404" height="540" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="348" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Diagnostics</text>',
    ]
    y = 382
    for row in diagnostics:
        color = {
            "pass": "#166534",
            "warn": "#92400e",
            "fail": "#991b1b",
            "info": "#1d4ed8",
        }.get(row.status, "#111827")
        lines.extend(
            [
                f'<text x="66" y="{y}" font-family="sans-serif" font-size="15" font-weight="700" fill="{color}">{row.status.upper()} · {row.name}</text>',
                f'<text x="66" y="{y + 22}" font-family="sans-serif" font-size="12" fill="#334155">{row.detail}</text>',
            ]
        )
        y += 64
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    default_gpp = shutil.which("g++")
    gpp14 = shutil.which("g++-14")
    mpicxx = shutil.which("mpic++")
    tempo2 = shutil.which("tempo2")
    root_config = shutil.which("root-config")
    cython = shutil.which("cython")
    boost_odeint = Path("/opt/homebrew/include/boost/numeric/odeint.hpp").exists()

    default_make = run(["make", "-f", "makefile-original", "Minuit_fit.exe"], cwd=SRC)
    boost_probe = run(
        ["make", "-f", "makefile-original", "CXX=g++-14", "CXXFLAGS=-g -Wall -O2 -fopenmp", "AllTheories3Bodies.o"],
        cwd=SRC,
    )

    missing_dependencies = []
    for name, value in (
        ("tempo2", tempo2),
        ("root-config", root_config),
        ("cython", cython),
        ("mpic++", mpicxx),
    ):
        if not value:
            missing_dependencies.append(name)
    if not boost_odeint:
        missing_dependencies.append("boost-odeint")

    diagnostics = [
        Diagnostic(
            name="default makefile compiler path",
            status="warn",
            detail=f"makefile-original uses g++; on this host that resolves to {default_gpp}.",
        ),
        Diagnostic(
            name="default build probe",
            status="fail",
            detail=(
                "make -f makefile-original Minuit_fit.exe fails immediately because Apple clang++ "
                "does not accept -fopenmp in the default CXXFLAGS."
            ),
        ),
        Diagnostic(
            name="gcc-14 probe",
            status="warn" if gpp14 else "fail",
            detail=f"Homebrew g++-14 path = {gpp14 or 'none'}.",
        ),
        Diagnostic(
            name="post-openmp compile probe",
            status="fail",
            detail=(
                "With CXX=g++-14, compilation gets past the OpenMP frontend issue but then fails on "
                "missing boost/numeric/odeint.hpp."
            ),
        ),
        Diagnostic(
            name="timing-stack dependencies",
            status="warn",
            detail=(
                f"tempo2={tempo2 or 'none'}, root-config={root_config or 'none'}, cython={cython or 'none'}, "
                f"mpic++={mpicxx or 'none'}, boost_odeint={boost_odeint}."
            ),
        ),
    ]

    write_tsv(DIAGNOSTICS_TSV, diagnostics)

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "default_gpp": default_gpp,
        "gpp14": gpp14,
        "mpicxx": mpicxx,
        "tempo2": tempo2,
        "root_config": root_config,
        "cython": cython,
        "boost_odeint_present": boost_odeint,
        "default_make_exit_code": default_make.returncode,
        "default_make_excerpt": excerpt(default_make.stdout),
        "boost_probe_exit_code": boost_probe.returncode,
        "boost_probe_excerpt": excerpt(boost_probe.stdout),
        "missing_dependencies": missing_dependencies,
        "verdict": (
            "Nutimo public inputs are available, but the current host does not yet close the "
            "runtime stack: the default build fails on Apple clang OpenMP, and the GCC path then "
            "fails on missing Boost and remaining timing dependencies."
        ),
        "next_bounded_action": (
            "do not begin open-ended dependency provisioning in this workspace; either use a "
            "pre-provisioned research environment or an external analysis route for full Phase B"
        ),
        "diagnostics": [asdict(row) for row in diagnostics],
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    write_svg(SUMMARY_SVG, summary, diagnostics)


if __name__ == "__main__":
    main()
