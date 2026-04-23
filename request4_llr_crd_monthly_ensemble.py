from __future__ import annotations

import csv
import hashlib
import json
import math
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from datetime import date, datetime, time, timedelta, timezone
from pathlib import Path
from statistics import median


SPEED_OF_LIGHT_M_S = 299_792_458.0


@dataclass(frozen=True)
class MonthlyCandidate:
    target: str
    yyyymm: str

    @property
    def year(self) -> str:
        return self.yyyymm[:4]

    @property
    def filename(self) -> str:
        return f"{self.target}_{self.yyyymm}.np2"

    @property
    def url(self) -> str:
        return (
            "https://edc.dgfi.tum.de/pub/slr/data/npt_crd_v2/"
            f"{self.target}/{self.year}/{self.filename}"
        )


@dataclass(frozen=True)
class EnsembleConfig:
    root_dir: str = "data/request4_llr/crd_monthly_ensemble_2026-04-23"
    rejected_dirname: str = "rejected_payloads"
    files_tsv: str = "request4_llr_crd_monthly_ensemble_files.tsv"
    rows_tsv: str = "request4_llr_crd_monthly_ensemble_normal_points.tsv"
    summary_json: str = "request4_llr_crd_monthly_ensemble_summary.json"
    summary_svg: str = "request4_llr_crd_monthly_ensemble_summary.svg"
    doc_md: str = "REQUEST4_LLR_CRD_MONTHLY_ENSEMBLE.md"
    targets: tuple[str, ...] = ("apollo15", "apollo11", "apollo14", "luna17", "luna21", "nglr1")
    months: tuple[str, ...] = (
        "202501",
        "202502",
        "202503",
        "202504",
        "202505",
        "202506",
        "202507",
        "202508",
        "202509",
        "202510",
        "202511",
        "202512",
        "202601",
        "202602",
        "202603",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_h4(tokens: list[str]) -> dict[str, object]:
    start_date = date(int(tokens[2]), int(tokens[3]), int(tokens[4]))
    stop_date = date(int(tokens[8]), int(tokens[9]), int(tokens[10]))
    start_sod = int(tokens[5]) * 3600 + int(tokens[6]) * 60 + int(tokens[7])
    stop_sod = int(tokens[11]) * 3600 + int(tokens[12]) * 60 + int(tokens[13])
    return {
        "start_date": start_date,
        "stop_date": stop_date,
        "start_sod": start_sod,
        "stop_sod": stop_sod,
    }


def sec_of_day_to_timestamp(base_date: date, base_start_sod: int, sec_of_day: float) -> str:
    day_offset = 0
    if sec_of_day + 12 * 3600 < base_start_sod:
        day_offset = 1
    day = base_date + timedelta(days=day_offset)
    whole_seconds = int(math.floor(sec_of_day))
    microseconds = int(round((sec_of_day - whole_seconds) * 1.0e6))
    if microseconds == 1_000_000:
        whole_seconds += 1
        microseconds = 0
    moment = datetime.combine(day, time(0, 0), tzinfo=timezone.utc) + timedelta(
        seconds=whole_seconds,
        microseconds=microseconds,
    )
    return moment.isoformat().replace("+00:00", "Z")


def parse_crd_np_file(path: Path, source_label: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    station_name = None
    station_code = None
    target_name = None
    target_id = None
    format_version = None
    h4 = None
    last_met = None

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            tokens = line.split()
            record_id = tokens[0].lower()
            if record_id == "h1":
                format_version = int(tokens[2])
            elif record_id == "h2":
                station_name = tokens[1]
                station_code = tokens[2]
            elif record_id == "h3":
                target_name = tokens[1]
                target_id = tokens[2]
            elif record_id == "h4":
                h4 = parse_h4(tokens)
            elif record_id == "20":
                last_met = {
                    "pressure_mbar": float(tokens[2]),
                    "temperature_k": float(tokens[3]),
                    "humidity_pct": float(tokens[4]),
                }
            elif record_id == "11":
                if h4 is None:
                    raise ValueError(f"11 record before h4 in {path}")
                sec_of_day = float(tokens[1])
                time_of_flight = float(tokens[2])
                signal_to_noise = None
                if len(tokens) > 13 and tokens[13].lower() != "na":
                    signal_to_noise = float(tokens[13])
                rows.append(
                    {
                        "source_label": source_label,
                        "source_file": path.name,
                        "crd_format_version": format_version,
                        "station_name": station_name,
                        "station_code": station_code,
                        "target_name": target_name,
                        "target_id": target_id,
                        "timestamp_utc": sec_of_day_to_timestamp(
                            h4["start_date"],
                            h4["start_sod"],
                            sec_of_day,
                        ),
                        "sec_of_day": sec_of_day,
                        "time_of_flight_s": time_of_flight,
                        "one_way_range_m": 0.5 * SPEED_OF_LIGHT_M_S * time_of_flight,
                        "epoch_event": int(tokens[4]),
                        "np_window_length_s": float(tokens[5]),
                        "num_ranges": int(tokens[6]),
                        "bin_rms_ps": parse_optional_float(tokens[7]),
                        "bin_skew": parse_optional_float(tokens[8]),
                        "bin_kurtosis": parse_optional_float(tokens[9]),
                        "bin_pmm_ps": parse_optional_float(tokens[10]),
                        "return_rate_hz": parse_optional_float(tokens[11]),
                        "detector_channel": int(tokens[12]),
                        "signal_to_noise": signal_to_noise,
                        "pressure_mbar": last_met["pressure_mbar"] if last_met else None,
                        "temperature_k": last_met["temperature_k"] if last_met else None,
                        "humidity_pct": last_met["humidity_pct"] if last_met else None,
                    }
                )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_optional_float(token: str) -> float | None:
    if token.lower() == "na":
        return None
    return float(token)


def fetch_head(candidate: MonthlyCandidate) -> dict[str, object]:
    request = urllib.request.Request(candidate.url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            return {
                "ok": True,
                "status": int(response.status),
                "content_length": int(response.headers.get("Content-Length", "0") or 0),
                "content_type": response.headers.get("Content-Type"),
            }
    except urllib.error.HTTPError as exc:
        return {"ok": False, "status": int(exc.code), "reason": str(exc.reason)}
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "status": None, "reason": repr(exc)}


def download_file(candidate: MonthlyCandidate, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(candidate.url, timeout=60) as response:
        payload = response.read()
    path.write_bytes(payload)


def looks_like_crd_payload(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open("rb") as handle:
        prefix = handle.read(64)
    return prefix.lower().startswith(b"h1 crd")


def build_candidates(config: EnsembleConfig) -> list[MonthlyCandidate]:
    return [MonthlyCandidate(target=target, yyyymm=yyyymm) for target in config.targets for yyyymm in config.months]


def render_summary_svg(summary: dict[str, object], path: Path) -> None:
    by_target = summary["by_target"]
    width = 860
    height = 360
    left = 180
    top = 44
    row_height = 34
    max_value = max(value["normal_point_rows"] for value in by_target.values()) or 1
    bar_max = 620
    rows = []
    for idx, (target, stats) in enumerate(by_target.items()):
        y = top + idx * row_height
        bar = int(round(bar_max * stats["normal_point_rows"] / max_value))
        rows.append(
            f'<text x="16" y="{y + 18}" font-size="15" fill="#102a43">{target}</text>'
            f'<rect x="{left}" y="{y}" width="{bar}" height="18" rx="4" fill="#1f6feb" />'
            f'<text x="{left + bar + 10}" y="{y + 14}" font-size="13" fill="#334e68">'
            f'{stats["normal_point_rows"]} rows / {stats["files"]} files</text>'
        )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
  <rect width="100%" height="100%" fill="#f8fbff" />
  <text x="16" y="24" font-size="20" font-weight="700" fill="#102a43">Request 4 Lunar CRD Monthly Ensemble</text>
  <text x="16" y="334" font-size="13" fill="#486581">files={summary['downloaded_files']}  normal_points={summary['total_normal_points']}  coverage={summary['coverage_start_utc']} .. {summary['coverage_stop_utc']}</text>
  {''.join(rows)}
</svg>
"""
    path.write_text(svg)


def main() -> None:
    config = EnsembleConfig()
    root = Path(config.root_dir)
    raw_dir = root / "raw"
    rejected_dir = root / config.rejected_dirname
    root.mkdir(parents=True, exist_ok=True)
    rejected_dir.mkdir(parents=True, exist_ok=True)

    candidates = build_candidates(config)
    file_rows: list[dict[str, object]] = []
    point_rows: list[dict[str, object]] = []

    for candidate in candidates:
        path = raw_dir / candidate.filename
        if looks_like_crd_payload(path):
            head = {
                "ok": True,
                "status": 200,
                "content_length": path.stat().st_size,
                "content_type": "text/plain",
                "reason": "local_cache",
            }
        else:
            head = fetch_head(candidate)
        file_row = {
            "target": candidate.target,
            "yyyymm": candidate.yyyymm,
            "year": candidate.year,
            "filename": candidate.filename,
            "url": candidate.url,
            "head_ok": head.get("ok"),
            "head_status": head.get("status"),
            "head_content_length": head.get("content_length"),
            "head_reason": head.get("reason"),
            "downloaded": False,
            "sha256": None,
            "size_bytes": None,
            "row_count": 0,
            "coverage_start_utc": None,
            "coverage_stop_utc": None,
            "stations": "",
        }
        if head.get("ok"):
            if not looks_like_crd_payload(path):
                download_file(candidate, path)
            if looks_like_crd_payload(path):
                rows = parse_crd_np_file(path, f"{candidate.target}_{candidate.yyyymm}")
                timestamps = [row["timestamp_utc"] for row in rows]
                stations = sorted({row["station_name"] for row in rows if row["station_name"]})
                file_row.update(
                    {
                        "downloaded": True,
                        "sha256": sha256_file(path),
                        "size_bytes": path.stat().st_size,
                        "row_count": len(rows),
                        "coverage_start_utc": min(timestamps) if timestamps else None,
                        "coverage_stop_utc": max(timestamps) if timestamps else None,
                        "stations": ",".join(stations),
                    }
                )
                point_rows.extend(rows)
            else:
                rejected_path = rejected_dir / candidate.filename
                path.replace(rejected_path)
                file_row["head_reason"] = "invalid_payload"
        file_rows.append(file_row)

    downloaded = [row for row in file_rows if row["downloaded"]]
    if not downloaded:
        raise SystemExit("no lunar monthly CRD files were downloaded")

    write_tsv(Path(config.files_tsv), file_rows)
    write_tsv(Path(config.rows_tsv), point_rows)

    by_target: dict[str, dict[str, int]] = {}
    by_station: dict[str, int] = {}
    for row in downloaded:
        target_stats = by_target.setdefault(row["target"], {"files": 0, "normal_point_rows": 0})
        target_stats["files"] += 1
        target_stats["normal_point_rows"] += int(row["row_count"])
    for row in point_rows:
        by_station[row["station_name"]] = by_station.get(row["station_name"], 0) + 1

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "dataset_root": config.root_dir,
        "downloaded_files": len(downloaded),
        "scanned_candidates": len(file_rows),
        "total_normal_points": len(point_rows),
        "coverage_start_utc": min(row["timestamp_utc"] for row in point_rows),
        "coverage_stop_utc": max(row["timestamp_utc"] for row in point_rows),
        "targets": sorted(by_target),
        "stations": sorted(by_station),
        "by_target": by_target,
        "by_station": dict(sorted(by_station.items(), key=lambda item: (-item[1], item[0]))),
        "median_one_way_range_m": median(row["one_way_range_m"] for row in point_rows),
        "downloaded_files_by_month": [row["yyyymm"] for row in downloaded],
        "max_file_row_count": max(int(row["row_count"]) for row in downloaded),
        "min_file_row_count": min(int(row["row_count"]) for row in downloaded),
        "downloaded_file_records": downloaded,
    }
    Path(config.summary_json).write_text(json.dumps(summary, indent=2))
    render_summary_svg(summary, Path(config.summary_svg))

    print("=== Request 4 lunar monthly CRD ensemble ===")
    print(f"downloaded_files = {summary['downloaded_files']} / {summary['scanned_candidates']}")
    print(f"total_normal_points = {summary['total_normal_points']}")
    print(f"coverage = {summary['coverage_start_utc']} .. {summary['coverage_stop_utc']}")
    print(f"targets = {summary['targets']}")
    print(f"stations = {summary['stations']}")
    print(f"median_one_way_range_m = {summary['median_one_way_range_m']:.3f}")
    print(f"max_file_row_count = {summary['max_file_row_count']}")
    print(f"min_file_row_count = {summary['min_file_row_count']}")
    print()
    print(f"wrote {config.files_tsv}")
    print(f"wrote {config.rows_tsv}")
    print(f"wrote {config.summary_json}")
    print(f"wrote {config.summary_svg}")


if __name__ == "__main__":
    main()
