from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parent
FILES_TSV = ROOT / "request4_llr_crd_monthly_ensemble_files.tsv"
ROWS_TSV = ROOT / "request4_llr_crd_monthly_ensemble_normal_points.tsv"
OUT_COVERAGE_TSV = ROOT / "request4_llr_crd_monthly_coverage.tsv"
OUT_SHORTLIST_TSV = ROOT / "request4_llr_crd_representative_cases.tsv"
OUT_SUMMARY_JSON = ROOT / "request4_llr_crd_coverage_summary.json"
OUT_SUMMARY_SVG = ROOT / "request4_llr_crd_coverage_summary.svg"


@dataclass(frozen=True)
class FileMetric:
    filename: str
    target: str
    yyyymm: str
    row_count: int
    dominant_station: str
    dominant_station_rows: int
    dominant_station_fraction: float
    station_count: int
    station_mix: str
    coverage_start_utc: str
    coverage_stop_utc: str
    coverage_span_days: float
    median_num_ranges: float
    median_bin_rms_ps: float | None
    median_bin_pmm_ps: float | None
    median_return_rate_hz: float | None
    median_one_way_range_m: float
    diversity_vector_norm: float


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_timestamp(token: str) -> datetime:
    if token.endswith("Z"):
        token = token[:-1] + "+00:00"
    return datetime.fromisoformat(token)


def parse_optional_float(token: str) -> float | None:
    if token in ("", "None", "na", "NA"):
        return None
    return float(token)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def safe_median(values: list[float]) -> float | None:
    if not values:
        return None
    return float(median(values))


def compute_file_metrics(files: list[dict[str, str]], rows: list[dict[str, str]]) -> list[FileMetric]:
    by_file: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_file[row["source_file"]].append(row)

    metrics: list[FileMetric] = []
    for file_row in files:
        if file_row["downloaded"] != "True":
            continue
        filename = file_row["filename"]
        file_points = by_file[filename]
        if not file_points:
            continue
        station_counts = Counter(point["station_name"] for point in file_points)
        dominant_station, dominant_rows = station_counts.most_common(1)[0]
        times = [parse_timestamp(point["timestamp_utc"]) for point in file_points]
        span_days = (max(times) - min(times)).total_seconds() / 86400.0 if len(times) > 1 else 0.0
        median_num_ranges = float(median(int(point["num_ranges"]) for point in file_points))
        median_bin_rms = safe_median(
            [value for value in (parse_optional_float(point["bin_rms_ps"]) for point in file_points) if value is not None]
        )
        median_bin_pmm = safe_median(
            [value for value in (parse_optional_float(point["bin_pmm_ps"]) for point in file_points) if value is not None]
        )
        median_return_rate = safe_median(
            [value for value in (parse_optional_float(point["return_rate_hz"]) for point in file_points) if value is not None]
        )
        median_range = float(median(float(point["one_way_range_m"]) for point in file_points))
        dominant_fraction = dominant_rows / len(file_points)
        diversity_vector_norm = math.sqrt(
            math.log1p(len(file_points)) ** 2
            + (len(station_counts) / 4.0) ** 2
            + (1.0 - dominant_fraction) ** 2
        )
        metrics.append(
            FileMetric(
                filename=filename,
                target=file_row["target"],
                yyyymm=file_row["yyyymm"],
                row_count=len(file_points),
                dominant_station=dominant_station,
                dominant_station_rows=dominant_rows,
                dominant_station_fraction=dominant_fraction,
                station_count=len(station_counts),
                station_mix=",".join(f"{station}:{count}" for station, count in station_counts.most_common()),
                coverage_start_utc=min(point["timestamp_utc"] for point in file_points),
                coverage_stop_utc=max(point["timestamp_utc"] for point in file_points),
                coverage_span_days=span_days,
                median_num_ranges=median_num_ranges,
                median_bin_rms_ps=median_bin_rms,
                median_bin_pmm_ps=median_bin_pmm,
                median_return_rate_hz=median_return_rate,
                median_one_way_range_m=median_range,
                diversity_vector_norm=diversity_vector_norm,
            )
        )
    return sorted(metrics, key=lambda metric: (metric.yyyymm, metric.target, metric.filename))


def metric_feature_vector(metric: FileMetric, station_order: list[str], target_order: list[str]) -> list[float]:
    month = int(metric.yyyymm[4:6])
    month_phase = 2.0 * math.pi * (month - 1) / 12.0
    features: list[float] = []
    features.extend(1.0 if metric.dominant_station == station else 0.0 for station in station_order)
    features.extend(1.0 if metric.target == target else 0.0 for target in target_order)
    features.extend(
        [
            math.sin(month_phase),
            math.cos(month_phase),
            math.log1p(metric.row_count),
            metric.station_count / 4.0,
            1.0 - metric.dominant_station_fraction,
            (metric.median_bin_rms_ps or 0.0) / 500.0,
            metric.median_num_ranges / 1000.0,
        ]
    )
    return features


def euclidean(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def combo_winners(metrics: list[FileMetric]) -> dict[tuple[str, str], FileMetric]:
    winners: dict[tuple[str, str], FileMetric] = {}
    for metric in metrics:
        key = (metric.dominant_station, metric.target)
        incumbent = winners.get(key)
        if incumbent is None or (
            metric.row_count,
            metric.station_count,
            -int(metric.yyyymm),
        ) > (
            incumbent.row_count,
            incumbent.station_count,
            -int(incumbent.yyyymm),
        ):
            winners[key] = metric
    return winners


def shortlist_representatives(metrics: list[FileMetric], shortlist_size: int = 8, min_rows: int = 5) -> list[FileMetric]:
    station_order = sorted({metric.dominant_station for metric in metrics})
    target_order = sorted({metric.target for metric in metrics})
    winners = combo_winners(metrics)
    candidate_pool = sorted(
        [metric for metric in winners.values() if metric.row_count >= min_rows],
        key=lambda metric: (
            -metric.row_count,
            -metric.station_count,
            metric.dominant_station_fraction,
            metric.target,
            metric.yyyymm,
        ),
    )
    if not candidate_pool:
        return []

    selected = [candidate_pool[0]]
    remaining = candidate_pool[1:]
    while remaining and len(selected) < min(shortlist_size, len(candidate_pool)):
        selected_vectors = [metric_feature_vector(metric, station_order, target_order) for metric in selected]
        best_metric = None
        best_score = -1.0
        for metric in remaining:
            vector = metric_feature_vector(metric, station_order, target_order)
            min_distance = min(euclidean(vector, selected_vector) for selected_vector in selected_vectors)
            score = min_distance + 0.015 * metric.row_count
            if score > best_score:
                best_score = score
                best_metric = metric
        assert best_metric is not None
        selected.append(best_metric)
        remaining = [metric for metric in remaining if metric.filename != best_metric.filename]
    return selected


def render_svg(summary: dict[str, object], path: Path) -> None:
    shortlist = summary["shortlist"]
    width = 980
    height = 420
    lines: list[str] = []
    for idx, item in enumerate(shortlist):
        y = 88 + idx * 34
        lines.append(
            f'<text x="24" y="{y}" font-size="14" fill="#102a43">{idx + 1}. '
            f'{item["filename"]}  [{item["dominant_station"]}/{item["target"]}]</text>'
        )
        lines.append(
            f'<text x="420" y="{y}" font-size="13" fill="#486581">'
            f'rows={item["row_count"]}  station_mix={item["station_mix"]}</text>'
        )
    bars: list[str] = []
    max_rows = max(item["normal_point_rows"] for item in summary["by_target"].values()) or 1
    for idx, (target, stats) in enumerate(summary["by_target"].items()):
        y = 64 + idx * 28
        bar = int(round(260 * stats["normal_point_rows"] / max_rows))
        bars.append(f'<text x="24" y="{y}" font-size="14" fill="#102a43">{target}</text>')
        bars.append(f'<rect x="120" y="{y-12}" width="{bar}" height="14" rx="4" fill="#1f6feb" />')
        bars.append(
            f'<text x="{128 + bar}" y="{y}" font-size="12" fill="#486581">'
            f'{stats["normal_point_rows"]} rows</text>'
        )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
  <rect width="100%" height="100%" fill="#f8fbff" />
  <text x="24" y="30" font-size="22" font-weight="700" fill="#102a43">Request 4 CRD Coverage Map</text>
  <text x="24" y="52" font-size="13" fill="#486581">valid files={summary["valid_files"]}  normal points={summary["total_normal_points"]}  stations={",".join(summary["stations"])}</text>
  {''.join(bars)}
  <text x="24" y="332" font-size="18" font-weight="700" fill="#102a43">Representative Case Shortlist</text>
  {''.join(lines)}
</svg>
"""
    path.write_text(svg)


def main() -> None:
    files = load_tsv(FILES_TSV)
    rows = load_tsv(ROWS_TSV)
    metrics = compute_file_metrics(files, rows)
    coverage_rows = [
        {
            "filename": metric.filename,
            "target": metric.target,
            "yyyymm": metric.yyyymm,
            "row_count": metric.row_count,
            "dominant_station": metric.dominant_station,
            "dominant_station_rows": metric.dominant_station_rows,
            "dominant_station_fraction": f"{metric.dominant_station_fraction:.6f}",
            "station_count": metric.station_count,
            "station_mix": metric.station_mix,
            "coverage_start_utc": metric.coverage_start_utc,
            "coverage_stop_utc": metric.coverage_stop_utc,
            "coverage_span_days": f"{metric.coverage_span_days:.6f}",
            "median_num_ranges": f"{metric.median_num_ranges:.1f}",
            "median_bin_rms_ps": "" if metric.median_bin_rms_ps is None else f"{metric.median_bin_rms_ps:.3f}",
            "median_bin_pmm_ps": "" if metric.median_bin_pmm_ps is None else f"{metric.median_bin_pmm_ps:.3f}",
            "median_return_rate_hz": "" if metric.median_return_rate_hz is None else f"{metric.median_return_rate_hz:.6f}",
            "median_one_way_range_m": f"{metric.median_one_way_range_m:.6f}",
            "diversity_vector_norm": f"{metric.diversity_vector_norm:.6f}",
        }
        for metric in metrics
    ]
    write_tsv(OUT_COVERAGE_TSV, coverage_rows)

    shortlist = shortlist_representatives(metrics, shortlist_size=8)
    shortlist_rows = [
        {
            "rank": idx + 1,
            "filename": metric.filename,
            "target": metric.target,
            "yyyymm": metric.yyyymm,
            "dominant_station": metric.dominant_station,
            "row_count": metric.row_count,
            "station_count": metric.station_count,
            "dominant_station_fraction": f"{metric.dominant_station_fraction:.6f}",
            "station_mix": metric.station_mix,
            "coverage_start_utc": metric.coverage_start_utc,
            "coverage_stop_utc": metric.coverage_stop_utc,
            "median_num_ranges": f"{metric.median_num_ranges:.1f}",
            "median_bin_rms_ps": "" if metric.median_bin_rms_ps is None else f"{metric.median_bin_rms_ps:.3f}",
            "median_bin_pmm_ps": "" if metric.median_bin_pmm_ps is None else f"{metric.median_bin_pmm_ps:.3f}",
            "median_return_rate_hz": "" if metric.median_return_rate_hz is None else f"{metric.median_return_rate_hz:.6f}",
            "selection_reason": "greedy_diversity_over_combo_winners",
        }
        for idx, metric in enumerate(shortlist)
    ]
    write_tsv(OUT_SHORTLIST_TSV, shortlist_rows)

    by_target: dict[str, dict[str, int]] = {}
    by_station: dict[str, int] = {}
    by_dominant_station: dict[str, dict[str, int]] = defaultdict(lambda: {"files": 0, "normal_point_rows": 0})
    for metric in metrics:
        target_stats = by_target.setdefault(metric.target, {"files": 0, "normal_point_rows": 0})
        target_stats["files"] += 1
        target_stats["normal_point_rows"] += metric.row_count
        station_stats = by_dominant_station[metric.dominant_station]
        station_stats["files"] += 1
        station_stats["normal_point_rows"] += metric.row_count
    for row in rows:
        by_station[row["station_name"]] = by_station.get(row["station_name"], 0) + 1

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "valid_files": len(metrics),
        "total_normal_points": len(rows),
        "coverage_start_utc": min(metric.coverage_start_utc for metric in metrics),
        "coverage_stop_utc": max(metric.coverage_stop_utc for metric in metrics),
        "targets": sorted(by_target),
        "stations": sorted(by_station),
        "by_target": by_target,
        "by_station": dict(sorted(by_station.items(), key=lambda item: (-item[1], item[0]))),
        "by_dominant_station": dict(sorted(by_dominant_station.items())),
        "shortlist": shortlist_rows,
        "coverage_gaps": [
            {
                "dominant_station": metric.dominant_station,
                "target": metric.target,
                "filename": metric.filename,
                "yyyymm": metric.yyyymm,
                "row_count": metric.row_count,
                "station_mix": metric.station_mix,
            }
            for metric in sorted(combo_winners(metrics).values(), key=lambda metric: (metric.row_count, metric.dominant_station, metric.target))
            if metric.row_count < 5
        ],
        "top_row_count_files": [
            {
                "filename": metric.filename,
                "target": metric.target,
                "yyyymm": metric.yyyymm,
                "dominant_station": metric.dominant_station,
                "row_count": metric.row_count,
                "station_mix": metric.station_mix,
            }
            for metric in sorted(metrics, key=lambda metric: (-metric.row_count, metric.filename))[:10]
        ],
    }
    OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    render_svg(summary, OUT_SUMMARY_SVG)

    print("=== Request 4 CRD coverage map ===")
    print(f"valid_files = {summary['valid_files']}")
    print(f"total_normal_points = {summary['total_normal_points']}")
    print(f"coverage = {summary['coverage_start_utc']} .. {summary['coverage_stop_utc']}")
    print(f"targets = {summary['targets']}")
    print(f"stations = {summary['stations']}")
    print("shortlist =")
    for row in shortlist_rows:
        print(
            f"  {row['rank']}. {row['filename']} "
            f"[{row['dominant_station']}/{row['target']}] rows={row['row_count']} "
            f"mix={row['station_mix']}"
        )
    print()
    print(f"wrote {OUT_COVERAGE_TSV.name}")
    print(f"wrote {OUT_SHORTLIST_TSV.name}")
    print(f"wrote {OUT_SUMMARY_JSON.name}")
    print(f"wrote {OUT_SUMMARY_SVG.name}")


if __name__ == "__main__":
    main()
