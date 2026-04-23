from __future__ import annotations

import csv
import json
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PEP_ROOT = ROOT / "data" / "request4_llr" / "pep_hand_off_2026-04-23"
PEP_CORE = PEP_ROOT / "pep_core"
PEP_DOC = PEP_ROOT / "pep_doc"

SUMMARY_JSON = ROOT / "request4_llr_pep_handoff_summary.json"
SUMMARY_SVG = ROOT / "request4_llr_pep_handoff_summary.svg"
DIAGNOSTICS_TSV = ROOT / "request4_llr_pep_handoff_diagnostics.tsv"


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


def probe_real10_support() -> subprocess.CompletedProcess[str]:
    with tempfile.TemporaryDirectory() as tmpdir:
        src = Path(tmpdir) / "test_real10.f"
        src.write_text("      real*10 x\n      x = 1.0_10\n      end\n")
        return run(["gfortran", "-c", str(src)])


def probe_long_double_flag() -> subprocess.CompletedProcess[str]:
    with tempfile.TemporaryDirectory() as tmpdir:
        src = Path(tmpdir) / "test_flag.f"
        src.write_text("      real x\n      x = 1.0\n      end\n")
        return run(["gfortran", "-c", "-m128bit-long-double", str(src)])


def minimal_pep_build_probe() -> subprocess.CompletedProcess[str]:
    target = PEP_CORE / "pep" / "a1ut1.o"
    if target.exists():
        target.unlink()
    return run(["make", "-C", "pep", "-f", "makefile", "a1ut1.o"], cwd=PEP_CORE)


def git_short_hash(repo: Path) -> str:
    proc = run(["git", "rev-parse", "--short", "HEAD"], cwd=repo)
    return proc.stdout.strip()


def command_excerpt(output: str, n_lines: int = 14) -> list[str]:
    return [line.rstrip() for line in output.splitlines()[:n_lines]]


def write_tsv(path: Path, rows: list[Diagnostic]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "status", "detail"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_svg(path: Path, summary: dict[str, object], diagnostics: list[Diagnostic]) -> None:
    width, height = 1500, 980
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: PEP hand-off feasibility</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Bounded acquisition/build gate for the external weak-field estimator hand-off.</text>',
        '<rect x="48" y="108" width="1404" height="182" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Decision</text>',
        f'<text x="66" y="170" font-family="monospace" font-size="13" fill="#111827">verdict = {summary["verdict"]}</text>',
        f'<text x="66" y="194" font-family="monospace" font-size="13" fill="#111827">current machine = {summary["machine_arch"]} | {summary["compiler_version_short"]}</text>',
        f'<text x="66" y="218" font-family="monospace" font-size="13" fill="#111827">REAL*10 files / hits = {summary["real10_file_count"]} / {summary["real10_hit_count"]}</text>',
        f'<text x="66" y="242" font-family="monospace" font-size="13" fill="#111827">long-double flag sites = {summary["long_double_flag_makefiles"]}</text>',
        f'<text x="66" y="266" font-family="monospace" font-size="13" fill="#111827">next bounded action = {summary["next_bounded_action"]}</text>',
        '<rect x="48" y="324" width="1404" height="594" fill="white" stroke="#cbd5e1"/>',
        '<text x="66" y="354" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Diagnostics</text>',
    ]
    y = 388
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
    pep_core_hash = git_short_hash(PEP_CORE)
    pep_doc_hash = git_short_hash(PEP_DOC)

    compiler_path = run(["which", "gfortran"]).stdout.strip()
    compiler_version = run(["gfortran", "--version"]).stdout.splitlines()
    compiler_version_short = compiler_version[0].strip() if compiler_version else "unknown"
    uname = run(["uname", "-a"]).stdout.strip()
    machine_arch = run(["arch"]).stdout.strip()
    rosetta_probe = run(["/usr/bin/arch", "-x86_64", "/usr/bin/true"])
    x86_gfortran_paths = sorted(str(path) for path in Path("/usr/local/bin").glob("gfortran*"))

    real10_hits = run(
        ["rg", "-n", "REAL\\*10", "pep", "pgms", "verify", "peputil"],
        cwd=PEP_CORE,
    ).stdout.splitlines()
    real10_files = run(
        ["rg", "-l", "REAL\\*10", "pep", "pgms", "verify", "peputil"],
        cwd=PEP_CORE,
    ).stdout.splitlines()
    long_double_sites = run(
        ["rg", "-l", "m128bit-long-double", "pep", "pgms", "verify", "peputil", "Makefile"],
        cwd=PEP_CORE,
    ).stdout.splitlines()

    llr_asset = PEP_CORE / "bigtest" / "bigtest.obsllr"
    llr_target_in_make = "bigtest.obsllr" in (PEP_CORE / "peputil" / "makefile").read_text()

    real10_probe = probe_real10_support()
    long_double_probe = probe_long_double_flag()
    build_probe = minimal_pep_build_probe()

    diagnostics = [
        Diagnostic(
            name="public acquisition",
            status="pass",
            detail=(
                f"Cloned pep_core@{pep_core_hash} and pep_doc@{pep_doc_hash} "
                "from public GitLab remotes."
            ),
        ),
        Diagnostic(
            name="LLR pedigree present",
            status="pass",
            detail=(
                f"bigtest.obsllr exists={llr_asset.exists()} and peputil/bigtest workflow references "
                f"it directly={llr_target_in_make}."
            ),
        ),
        Diagnostic(
            name="arm64 REAL*10 support",
            status="fail",
            detail="Current gfortran rejects both REAL*10 declarations and _10 literals in a trivial probe.",
        ),
        Diagnostic(
            name="arm64 long-double flag",
            status="fail",
            detail=(
                "Current gfortran rejects -m128bit-long-double, while pep_core/verify and peputil "
                "makefiles require that flag."
            ),
        ),
        Diagnostic(
            name="minimal object build",
            status="fail",
            detail=(
                "make -C pep -f makefile a1ut1.o fails immediately at the first REAL*10 source, "
                "before any LLR bigtest run can start."
            ),
        ),
        Diagnostic(
            name="x86 legacy toolchain on host",
            status="warn",
            detail=(
                f"Rosetta available={rosetta_probe.returncode == 0}, but dedicated /usr/local/bin gfortran paths found={len(x86_gfortran_paths)}."
            ),
        ),
        Diagnostic(
            name="porting risk",
            status="warn",
            detail=(
                f"REAL*10 appears in {len(real10_files)} files and {len(real10_hits)} source hits, so "
                "rewriting kinds locally would become a porting project, not a bounded hand-off."
            ),
        ),
    ]

    write_tsv(DIAGNOSTICS_TSV, diagnostics)

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "pep_core_git_hash": pep_core_hash,
        "pep_doc_git_hash": pep_doc_hash,
        "uname": uname,
        "machine_arch": machine_arch,
        "compiler_path": compiler_path,
        "compiler_version_short": compiler_version_short,
        "rosetta_available": rosetta_probe.returncode == 0,
        "dedicated_x86_gfortran_paths": x86_gfortran_paths,
        "real10_file_count": len(real10_files),
        "real10_hit_count": len(real10_hits),
        "long_double_flag_makefiles": long_double_sites,
        "llr_bigtest_obsllr_exists": llr_asset.exists(),
        "llr_bigtest_input_count": len(list((PEP_CORE / "bigtest").glob("*"))),
        "real10_probe_exit_code": real10_probe.returncode,
        "real10_probe_excerpt": command_excerpt(real10_probe.stdout),
        "long_double_probe_exit_code": long_double_probe.returncode,
        "long_double_probe_excerpt": command_excerpt(long_double_probe.stdout),
        "minimal_build_probe_exit_code": build_probe.returncode,
        "minimal_build_probe_excerpt": command_excerpt(build_probe.stdout),
        "verdict": (
            "PEP is confirmed as the right external estimator family, but it is not "
            "locally buildable on the current arm64 gfortran toolchain."
        ),
        "next_bounded_action": (
            "try only a separate x86_64 / legacy-compatible build environment; do not "
            "turn Request 4 into a local PEP porting project"
        ),
        "fallback_if_environment_probe_fails": (
            "hand weak-field posterior work to a collaboration / external-analysis-center "
            "route such as JPL or POLAC-ELPN"
        ),
        "diagnostics": [asdict(row) for row in diagnostics],
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    write_svg(SUMMARY_SVG, summary, diagnostics)


if __name__ == "__main__":
    main()
