from __future__ import annotations

import csv
import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SUMMARY_JSON = ROOT / "request4_llr_pep_environment_summary.json"
SUMMARY_SVG = ROOT / "request4_llr_pep_environment_summary.svg"
DIAGNOSTICS_TSV = ROOT / "request4_llr_pep_environment_diagnostics.tsv"


@dataclass(frozen=True)
class Diagnostic:
    name: str
    status: str
    detail: str


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def write_tsv(path: Path, rows: list[Diagnostic]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "status", "detail"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_svg(path: Path, summary: dict[str, object], diagnostics: list[Diagnostic]) -> None:
    width, height = 1500, 860
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: PEP environment probe</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Bounded check for a local x86_64 or legacy-compatible PEP build path.</text>',
        '<rect x="48" y="108" width="1404" height="162" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Decision</text>',
        f'<text x="66" y="170" font-family="monospace" font-size="13" fill="#111827">verdict = {summary["verdict"]}</text>',
        f'<text x="66" y="194" font-family="monospace" font-size="13" fill="#111827">rosetta_available = {summary["rosetta_available"]}</text>',
        f'<text x="66" y="218" font-family="monospace" font-size="13" fill="#111827">x86_homebrew_present = {summary["x86_homebrew_present"]} | x86_gfortran_count = {summary["x86_gfortran_count"]}</text>',
        f'<text x="66" y="242" font-family="monospace" font-size="13" fill="#111827">container_runtimes = {summary["detected_container_runtimes"]}</text>',
        '<rect x="48" y="304" width="1404" height="500" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="334" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Diagnostics</text>',
    ]
    y = 368
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
    rosetta = run(["/usr/bin/arch", "-x86_64", "/usr/bin/true"])
    x86_homebrew = Path("/usr/local/Homebrew").exists() or Path("/usr/local/bin/brew").exists()
    x86_gfortran_paths = sorted(str(path) for path in Path("/usr/local/bin").glob("gfortran*"))
    x86_gcc_paths = sorted(str(path) for path in Path("/usr/local/bin").glob("gcc*"))

    runtimes = {}
    for name in ("docker", "podman", "colima", "lima", "nerdctl", "orb"):
        path = shutil.which(name)
        runtimes[name] = path

    diagnostics = [
        Diagnostic(
            name="rosetta runtime",
            status="pass" if rosetta.returncode == 0 else "fail",
            detail=f"arch -x86_64 /usr/bin/true exit code = {rosetta.returncode}.",
        ),
        Diagnostic(
            name="x86 Homebrew prefix",
            status="warn" if not x86_homebrew else "pass",
            detail=f"/usr/local Homebrew present = {x86_homebrew}.",
        ),
        Diagnostic(
            name="x86 gfortran toolchain",
            status="warn" if not x86_gfortran_paths else "pass",
            detail=f"Detected /usr/local/bin gfortran paths = {x86_gfortran_paths or 'none'}.",
        ),
        Diagnostic(
            name="x86 gcc toolchain",
            status="warn" if not x86_gcc_paths else "pass",
            detail=f"Detected /usr/local/bin gcc paths = {x86_gcc_paths or 'none'}.",
        ),
        Diagnostic(
            name="container/VM runtime",
            status="warn" if not any(runtimes.values()) else "pass",
            detail=(
                "Detected runtime paths = "
                + ", ".join(f"{name}:{path or 'none'}" for name, path in runtimes.items())
            ),
        ),
    ]

    write_tsv(DIAGNOSTICS_TSV, diagnostics)

    detected_container_runtimes = [name for name, path in runtimes.items() if path]
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "rosetta_available": rosetta.returncode == 0,
        "x86_homebrew_present": x86_homebrew,
        "x86_gfortran_count": len(x86_gfortran_paths),
        "x86_gfortran_paths": x86_gfortran_paths,
        "x86_gcc_count": len(x86_gcc_paths),
        "x86_gcc_paths": x86_gcc_paths,
        "detected_container_runtimes": detected_container_runtimes,
        "runtime_paths": runtimes,
        "verdict": (
            "No bounded local x86_64 or legacy-compatible PEP build path is presently "
            "available on this host without provisioning a new toolchain or runtime."
        ),
        "project_consequence": (
            "Stop local PEP escalation here; any PEP hand-off now requires a separate "
            "pre-provisioned environment or an external analysis route."
        ),
        "diagnostics": [asdict(row) for row in diagnostics],
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    write_svg(SUMMARY_SVG, summary, diagnostics)


if __name__ == "__main__":
    main()
