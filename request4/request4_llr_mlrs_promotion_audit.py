from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent
SHORTLIST_TSV = ROOT / "request4_llr_crd_representative_cases.tsv"
MLRS_LAB = REPO_ROOT / "data/request4_llr/mlrs_handshake_lab"
CONVERTER = (
    ROOT
    / "data/request4_llr/ilrs_pivot_scout/crd_sample_code_v2.01b/crd_conv_v1v2/crd_conv"
)
SPLITTER = (
    ROOT
    / "data/request4_llr/ilrs_pivot_scout/crd_sample_code_v2.01b/crd_split_c/crd_split"
)
ROOT_OUT = REPO_ROOT / "data/request4_llr/mlrs_promotion_audit_2026-04-23"
RAW_DIR = ROOT_OUT / "raw_fr2"
WORK_DIR = ROOT_OUT / "work"
SUMMARY_TSV = ROOT / "request4_llr_mlrs_promotion_audit.tsv"
SUMMARY_JSON = ROOT / "request4_llr_mlrs_promotion_audit_summary.json"
SUMMARY_SVG = ROOT / "request4_llr_mlrs_promotion_audit_summary.svg"
DOC_MD = ROOT / "REQUEST4_LLR_MLRS_PROMOTION_AUDIT.md"


@dataclass(frozen=True)
class RepresentativeCase:
    filename: str
    target: str
    yyyymm: str
    dominant_station: str
    row_count: int

    @property
    def year(self) -> str:
        return self.yyyymm[:4]

    @property
    def np2_filename(self) -> str:
        return self.filename

    @property
    def fr2_filename(self) -> str:
        return self.filename.replace(".np2", ".fr2")

    @property
    def fr2_url(self) -> str:
        return f"https://edc.dgfi.tum.de/pub/slr/data/fr_crd_v2/{self.target}/{self.year}/{self.fr2_filename}"


@dataclass(frozen=True)
class PromotionResult:
    filename: str
    target: str
    yyyymm: str
    dominant_station: str
    row_count: int
    fr2_url: str
    fr2_sha256: str
    monthly_h4_count: int
    split_pass_count: int
    chosen_pass_filename: str
    chosen_pass_h2: str
    chosen_pass_h3: str
    chosen_pass_c0: str
    v1_h2: str
    v1_h3: str
    normalized_c0: str
    recalc_passed: bool
    poisson_passed: bool
    frfil_has_93: int
    frfil_has_30: int
    frfil_has_12: int
    error_signature: str
    note: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_shortlist(limit: int = 8) -> list[RepresentativeCase]:
    cases: list[RepresentativeCase] = []
    with SHORTLIST_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cases.append(
                RepresentativeCase(
                    filename=row["filename"],
                    target=row["target"],
                    yyyymm=row["yyyymm"],
                    dominant_station=row["dominant_station"],
                    row_count=int(row["row_count"]),
                )
            )
            if len(cases) >= limit:
                break
    return cases


def ensure_download(case: RepresentativeCase) -> Path:
    target_dir = RAW_DIR / case.target / case.year
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / case.fr2_filename
    if path.exists():
        return path
    with urllib.request.urlopen(case.fr2_url, timeout=120) as response:
        payload = response.read()
    path.write_bytes(payload)
    return path


def looks_like_crd(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open("rb") as handle:
        return handle.read(32).lower().startswith(b"h1 crd")


def count_tags(path: Path, tags: tuple[str, ...]) -> dict[str, int]:
    counts = {tag.upper(): 0 for tag in tags}
    with path.open() as handle:
        for line in handle:
            tag = line[:2].upper()
            if tag in counts:
                counts[tag] += 1
    return counts


def run_command(args: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=cwd,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )


def split_monthly(case: RepresentativeCase, fr2_path: Path) -> tuple[Path, list[Path]]:
    case_dir = WORK_DIR / case.filename.replace(".np2", "")
    split_dir = case_dir / "split"
    if split_dir.exists():
        shutil.rmtree(split_dir)
    split_dir.mkdir(parents=True, exist_ok=True)
    result = run_command([str(SPLITTER), str(fr2_path)], split_dir)
    (split_dir / "split.log").write_text(result.stdout)
    passes = sorted(path for path in split_dir.glob("*.frd"))
    return split_dir, passes


def inspect_header(pass_file: Path) -> dict[str, str]:
    info = {"h2": "", "h3": "", "c0": "", "station_name": ""}
    with pass_file.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(("H2 ", "h2 ")):
                info["h2"] = line
                parts = line.split()
                if len(parts) > 1:
                    info["station_name"] = parts[1]
            elif line.startswith(("H3 ", "h3 ")):
                info["h3"] = line
            elif line.startswith(("C0 ", "c0 ")):
                info["c0"] = line
            if all(info[key] for key in ("h2", "h3", "c0")):
                break
    return info


def choose_pass(pass_files: list[Path], dominant_station: str) -> Path:
    if not pass_files:
        raise RuntimeError("No pass files produced by crd_split")
    exact_station: list[Path] = []
    for pass_file in pass_files:
        if inspect_header(pass_file)["station_name"] == dominant_station:
            exact_station.append(pass_file)
    pool = exact_station if exact_station else pass_files
    return sorted(pool)[0]


def convert_to_v1(split_pass: Path) -> Path:
    out_path = split_pass.with_suffix(".v1.frd")
    run = run_command([str(CONVERTER), str(split_pass), str(out_path)], split_pass.parent)
    (split_pass.parent / f"{split_pass.stem}.conv.log").write_text(run.stdout)
    if run.returncode != 0:
        raise RuntimeError(f"crd_conv failed for {split_pass.name}")
    return out_path


def normalize_mlrs_subset(v1_path: Path) -> tuple[Path, dict[str, str]]:
    out_path = v1_path.with_suffix(".mlrs.v1.frd")
    chosen_h2 = ""
    chosen_h3 = ""
    chosen_c0 = ""
    v1_h2 = ""
    v1_h3 = ""
    normalized_c0 = ""
    out_lines: list[str] = []
    for raw_line in v1_path.read_text().splitlines():
        line = raw_line
        if line.startswith(("H2 ", "h2 ")):
            chosen_h2 = line.strip() if not chosen_h2 else chosen_h2
            v1_h2 = line.strip()
        elif line.startswith(("H3 ", "h3 ")):
            chosen_h3 = line.strip() if not chosen_h3 else chosen_h3
            v1_h3 = line.strip()
        elif line.startswith(("C0 ", "c0 ")):
            chosen_c0 = line.strip() if not chosen_c0 else chosen_c0
            tokens = line.split()
            head = tokens[:3]
            ids = [token for token in tokens[3:] if token.lower() != "na"]
            ids = ids[:5]
            line = " ".join([head[0].lower(), head[1], head[2], *ids])
            normalized_c0 = line
        elif line.startswith(("C6 ", "c6 ", "C7 ", "c7 ", "41 ", "42 ")):
            continue
        out_lines.append(line.lower() if line[:2] in {"H1", "H2", "H3", "H4"} else line)
    out_path.write_text("\n".join(out_lines) + "\n")
    meta = {
        "chosen_h2": chosen_h2,
        "chosen_h3": chosen_h3,
        "chosen_c0": chosen_c0,
        "v1_h2": v1_h2,
        "v1_h3": v1_h3,
        "normalized_c0": normalized_c0,
    }
    return out_path, meta


def stage_and_run(case_name: str, normalized_input: Path) -> tuple[bool, bool, dict[str, int], str]:
    analysis_dir = MLRS_LAB / "data/analysis"
    frcal = analysis_dir / f"{case_name}.frcal"
    shutil.copy2(normalized_input, frcal)
    (analysis_dir / f"{case_name}.llc").write_text("")
    (analysis_dir / f"{case_name}.lsc").write_text("")
    for suffix in [
        "errorc",
        "frrec",
        "rsc",
        "npt",
        "frd",
        "frfin",
        "frfil",
        "frfil2",
        "frfil3",
    ]:
        path = analysis_dir / f"{case_name}.{suffix}"
        if path.exists():
            path.unlink()
    run = run_command(["csh", "-f", "bin/ldb_crd", case_name], MLRS_LAB)
    (analysis_dir / f"{case_name}.run.log").write_text(run.stdout)
    error_text = ""
    error_path = analysis_dir / f"{case_name}.errorc"
    if error_path.exists():
        error_text = error_path.read_text().strip()
    recalc_passed = "Recalc Failed." not in run.stdout and "Recalc Failed." not in error_text
    poisson_passed = "Poisson Failed." not in run.stdout and "Poisson Failed." not in error_text
    frfil = analysis_dir / f"{case_name}.frfil"
    counts = {"12": 0, "30": 0, "93": 0}
    if frfil.exists():
        counts = count_tags(frifil := frfil, ("12", "30", "93"))
    if error_text:
        signature = error_text.splitlines()[-1]
    elif "Poisson Failed." in run.stdout:
        signature = "Poisson Failed."
    elif run.returncode != 0:
        signature = f"ldb_crd exit {run.returncode}"
    else:
        signature = "ok"
    return recalc_passed, poisson_passed, counts, signature


def render_svg(results: list[PromotionResult]) -> None:
    width = 1200
    height = 180 + 36 * len(results)
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#f8fbff" />',
        '<text x="24" y="30" font-size="22" font-weight="700" fill="#102a43">Request 4 MLRS Promotion Audit</text>',
        '<text x="24" y="54" font-size="13" fill="#486581">monthly fr2 -> split -> v1 -> C0 normalization -> ldb_crd</text>',
    ]
    for idx, result in enumerate(results):
        y = 92 + idx * 36
        status = "recalc-pass" if result.recalc_passed else "recalc-fail"
        if result.recalc_passed and not result.poisson_passed:
            status = "recalc-only"
        lines.append(
            f'<text x="24" y="{y}" font-size="14" fill="#102a43">{result.filename} '
            f'[{result.dominant_station}/{result.target}]</text>'
        )
        lines.append(
            f'<text x="420" y="{y}" font-size="13" fill="#486581">split={result.split_pass_count} '
            f'93={result.frfil_has_93} 30={result.frfil_has_30} 12={result.frfil_has_12} '
            f'status={status}</text>'
        )
        lines.append(
            f'<text x="780" y="{y}" font-size="12" fill="#7b8794">{result.error_signature}</text>'
        )
    lines.append("</svg>")
    SUMMARY_SVG.write_text("\n".join(lines) + "\n")


def write_tsv(results: list[PromotionResult]) -> None:
    rows = [asdict(result) for result in results]
    with SUMMARY_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ROOT_OUT.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    results: list[PromotionResult] = []
    promotion_target = 3
    promoted = 0
    for case in load_shortlist():
        fr2_path = ensure_download(case)
        if not looks_like_crd(fr2_path):
            results.append(
                PromotionResult(
                    filename=case.filename,
                    target=case.target,
                    yyyymm=case.yyyymm,
                    dominant_station=case.dominant_station,
                    row_count=case.row_count,
                    fr2_url=case.fr2_url,
                    fr2_sha256=sha256_file(fr2_path),
                    monthly_h4_count=0,
                    split_pass_count=0,
                    chosen_pass_filename="",
                    chosen_pass_h2="",
                    chosen_pass_h3="",
                    chosen_pass_c0="",
                    v1_h2="",
                    v1_h3="",
                    normalized_c0="",
                    recalc_passed=False,
                    poisson_passed=False,
                    frfil_has_93=0,
                    frfil_has_30=0,
                    frfil_has_12=0,
                    error_signature="invalid_fullrate_payload",
                    note="public fr2 URL did not return a CRD payload",
                )
            )
            continue

        monthly_counts = count_tags(fr2_path, ("H4",))
        split_dir, pass_files = split_monthly(case, fr2_path)
        if not pass_files:
            results.append(
                PromotionResult(
                    filename=case.filename,
                    target=case.target,
                    yyyymm=case.yyyymm,
                    dominant_station=case.dominant_station,
                    row_count=case.row_count,
                    fr2_url=case.fr2_url,
                    fr2_sha256=sha256_file(fr2_path),
                    monthly_h4_count=monthly_counts["H4"],
                    split_pass_count=0,
                    chosen_pass_filename="",
                    chosen_pass_h2="",
                    chosen_pass_h3="",
                    chosen_pass_c0="",
                    v1_h2="",
                    v1_h3="",
                    normalized_c0="",
                    recalc_passed=False,
                    poisson_passed=False,
                    frfil_has_93=0,
                    frfil_has_30=0,
                    frfil_has_12=0,
                    error_signature="no_split_passes",
                    note="crd_split produced no single-pass files from the monthly fr2 payload",
                )
            )
            continue

        chosen_pass = choose_pass(pass_files, case.dominant_station)
        v1_path = convert_to_v1(chosen_pass)
        normalized_path, meta = normalize_mlrs_subset(v1_path)
        case_name = f"{case.target}_{case.yyyymm}_{case.dominant_station.lower()}_promo"
        recalc_passed, poisson_passed, fr_counts, signature = stage_and_run(case_name, normalized_path)
        note = (
            "recalc passed but no 93 stream reached Poisson"
            if recalc_passed and fr_counts["93"] == 0
            else "promotion advanced past parser stage"
        )
        results.append(
            PromotionResult(
                filename=case.filename,
                target=case.target,
                yyyymm=case.yyyymm,
                dominant_station=case.dominant_station,
                row_count=case.row_count,
                fr2_url=case.fr2_url,
                fr2_sha256=sha256_file(fr2_path),
                monthly_h4_count=monthly_counts["H4"],
                split_pass_count=len(pass_files),
                chosen_pass_filename=chosen_pass.name,
                chosen_pass_h2=meta["chosen_h2"],
                chosen_pass_h3=meta["chosen_h3"],
                chosen_pass_c0=meta["chosen_c0"],
                v1_h2=meta["v1_h2"],
                v1_h3=meta["v1_h3"],
                normalized_c0=meta["normalized_c0"],
                recalc_passed=recalc_passed,
                poisson_passed=poisson_passed,
                frfil_has_93=fr_counts["93"],
                frfil_has_30=fr_counts["30"],
                frfil_has_12=fr_counts["12"],
                error_signature=signature,
                note=note,
            )
        )
        promoted += 1
        if promoted >= promotion_target:
            break

    write_tsv(results)
    SUMMARY_JSON.write_text(
        json.dumps(
            {
                "generated_at_utc": subprocess.run(
                    ["python3", "-c", "from datetime import datetime,timezone; print(datetime.now(timezone.utc).isoformat())"],
                    check=True,
                    text=True,
                    stdout=subprocess.PIPE,
                ).stdout.strip(),
                "promotion_results": [asdict(result) for result in results],
            },
            indent=2,
        )
        + "\n"
    )
    render_svg(results)


if __name__ == "__main__":
    main()
