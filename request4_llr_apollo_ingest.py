from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import urllib.request
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path


SPEED_OF_LIGHT_M_S = 299_792_458.0


@dataclass(frozen=True)
class ApolloBatch:
    label: str
    url: str
    release_date: str
    coverage_start: str
    coverage_end: str
    expected_count: int
    notes: str


@dataclass(frozen=True)
class ApolloRelease:
    release_id: str = "apollo_legacy_text_release_2026-04-22"
    source_page_url: str = "https://tmurphy.physics.ucsd.edu/apollo/norm_pts.html"
    source_page_note: str = (
        "Pinned from the APOLLO normal-point page as accessed on 2026-04-22. "
        "This is the older legacy text format described on that page, not the "
        "newer CRD stream mirrored through CDDIS."
    )
    batches: tuple[ApolloBatch, ...] = (
        ApolloBatch(
            label="group_a",
            url="https://tmurphy.physics.ucsd.edu/apollo/121106a_apollo_np.txt",
            release_date="2012-11-06",
            coverage_start="2006-04-07",
            coverage_end="2010-10-30",
            expected_count=941,
            notes="Initial long baseline release from the APOLLO archive page.",
        ),
        ApolloBatch(
            label="group_b",
            url="https://tmurphy.physics.ucsd.edu/apollo/121106b_apollo_np.txt",
            release_date="2012-11-06",
            coverage_start="2010-12-01",
            coverage_end="2012-04-06",
            expected_count=506,
            notes="Second batch from the APOLLO archive page.",
        ),
        ApolloBatch(
            label="group_c_replacement",
            url="https://tmurphy.physics.ucsd.edu/apollo/131219c_apollo_np.txt",
            release_date="2013-12-19",
            coverage_start="2012-04-06",
            coverage_end="2013-09-01",
            expected_count=361,
            notes="Replacement group c file listed on the APOLLO archive page.",
        ),
        ApolloBatch(
            label="group_d_complete",
            url="https://tmurphy.physics.ucsd.edu/apollo/160725d_apollo_np.txt",
            release_date="2016-07-25",
            coverage_start="2013-09-30",
            coverage_end="2016-07-25",
            expected_count=760,
            notes="Complete group d file superseding earlier d extensions.",
        ),
        ApolloBatch(
            label="group_e_complete",
            url="https://tmurphy.physics.ucsd.edu/apollo/160912e_apollo_np2.txt",
            release_date="2016-09-12",
            coverage_start="2016-09-12",
            coverage_end="2019-08-18",
            expected_count=915,
            notes="Group e complete file as listed on the APOLLO archive page.",
        ),
        ApolloBatch(
            label="group_f",
            url="https://tmurphy.physics.ucsd.edu/apollo/190823f_apollo_np2.txt",
            release_date="2020-12-27",
            coverage_start="2019-08-23",
            coverage_end="2020-12-27",
            expected_count=355,
            notes="Latest legacy-text batch listed at the top of the APOLLO archive page.",
        ),
    )


REFLECTOR_NAMES = {
    0: "Apollo11",
    1: "Luna17",
    2: "Apollo14",
    3: "Apollo15",
    4: "Luna21",
}


def fetch_bytes(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def parse_legacy_line(line: str, batch_label: str, source_name: str) -> dict[str, object]:
    if len(line) < 79:
        raise ValueError(f"line too short ({len(line)}): {line!r}")
    if not line.startswith("51"):
        raise ValueError(f"unexpected record type in {line!r}")

    year = int(line[2:6].strip())
    month = int(line[6:8].strip())
    day = int(line[8:10].strip())
    hour = int(line[10:12].strip())
    minute = int(line[12:14].strip())
    second_1e7 = int(line[14:23].strip())
    second_total = second_1e7 / 1.0e7
    whole_seconds = int(math.floor(second_total))
    microsecond = int(round((second_total - whole_seconds) * 1.0e6))
    if microsecond == 1_000_000:
        whole_seconds += 1
        microsecond = 0
    timestamp = datetime(year, month, day, hour, minute, tzinfo=timezone.utc) + timedelta(
        seconds=whole_seconds,
        microseconds=microsecond,
    )

    round_trip_s = int(line[23:37].strip()) / 1.0e13
    reflector_id = int(line[37:38].strip())
    station_id = int(line[38:43].strip())
    photon_count = int(line[43:46].strip())
    uncertainty_ps = int(line[46:52].strip()) / 10.0
    snr = int(line[52:55].strip()) / 10.0
    quality = line[55:56]
    pressure_mbar = int(line[56:62].strip()) / 100.0
    temperature_c = int(line[62:66].strip()) / 10.0
    humidity_pct = int(line[66:68].strip())
    wavelength_angstrom = int(line[68:73].strip())
    wavelength_code = line[73:74]
    np_duration_s = int(line[74:78].strip())
    duration_code = line[78:79]

    one_way_range_m = 0.5 * SPEED_OF_LIGHT_M_S * round_trip_s

    return {
        "batch_label": batch_label,
        "source_file": source_name,
        "record_type": "51",
        "timestamp_utc": timestamp.isoformat().replace("+00:00", "Z"),
        "year": year,
        "month": month,
        "day": day,
        "hour": hour,
        "minute": minute,
        "second_of_minute": second_total,
        "round_trip_s": round_trip_s,
        "one_way_range_m": one_way_range_m,
        "reflector_id": reflector_id,
        "reflector_name": REFLECTOR_NAMES.get(reflector_id, f"unknown_{reflector_id}"),
        "station_id": station_id,
        "photon_count": photon_count,
        "uncertainty_ps": uncertainty_ps,
        "snr": snr,
        "quality": quality,
        "pressure_mbar": pressure_mbar,
        "temperature_c": temperature_c,
        "humidity_pct": humidity_pct,
        "wavelength_angstrom": wavelength_angstrom,
        "wavelength_code": wavelength_code,
        "np_duration_s": np_duration_s,
        "duration_code": duration_code,
        "raw_line": line,
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "batch_label",
        "source_file",
        "record_type",
        "timestamp_utc",
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second_of_minute",
        "round_trip_s",
        "one_way_range_m",
        "reflector_id",
        "reflector_name",
        "station_id",
        "photon_count",
        "uncertainty_ps",
        "snr",
        "quality",
        "pressure_mbar",
        "temperature_c",
        "humidity_pct",
        "wavelength_angstrom",
        "wavelength_code",
        "np_duration_s",
        "duration_code",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})


def summarize_rows(rows: list[dict[str, object]], release: ApolloRelease, retrieval_date: str, downloads: list[dict[str, object]]) -> dict[str, object]:
    by_batch: dict[str, list[dict[str, object]]] = defaultdict(list)
    by_reflector = Counter()
    by_quality = Counter()
    by_year = Counter()
    station_ids = Counter()
    for row in rows:
        by_batch[str(row["batch_label"])].append(row)
        by_reflector[str(row["reflector_name"])] += 1
        by_quality[str(row["quality"])] += 1
        by_year[int(row["year"])] += 1
        station_ids[int(row["station_id"])] += 1

    def stats(values: list[float]) -> dict[str, float]:
        arr = np.asarray(values, dtype=float)
        return {
            "min": float(np.min(arr)),
            "max": float(np.max(arr)),
            "mean": float(np.mean(arr)),
            "median": float(np.median(arr)),
        }

    import numpy as np

    sorted_rows = sorted(rows, key=lambda item: item["timestamp_utc"])
    batch_summaries = []
    for batch in release.batches:
        batch_rows = by_batch[batch.label]
        batch_summaries.append(
            {
                "label": batch.label,
                "release_date": batch.release_date,
                "coverage_start": batch.coverage_start,
                "coverage_end": batch.coverage_end,
                "expected_count": batch.expected_count,
                "parsed_count": len(batch_rows),
                "count_matches_page": len(batch_rows) == batch.expected_count,
                "first_timestamp_utc": batch_rows[0]["timestamp_utc"] if batch_rows else None,
                "last_timestamp_utc": batch_rows[-1]["timestamp_utc"] if batch_rows else None,
                "reflector_counts": dict(Counter(str(row["reflector_name"]) for row in batch_rows)),
            }
        )

    return {
        "release": asdict(release),
        "retrieval_date": retrieval_date,
        "downloads": downloads,
        "total_count": len(rows),
        "first_timestamp_utc": sorted_rows[0]["timestamp_utc"],
        "last_timestamp_utc": sorted_rows[-1]["timestamp_utc"],
        "station_ids": dict(sorted(station_ids.items())),
        "reflector_counts": dict(by_reflector),
        "quality_counts": dict(by_quality),
        "year_counts": dict(sorted(by_year.items())),
        "photon_count_stats": stats([float(row["photon_count"]) for row in rows]),
        "uncertainty_ps_stats": stats([float(row["uncertainty_ps"]) for row in rows]),
        "snr_stats": stats([float(row["snr"]) for row in rows]),
        "one_way_range_m_stats": stats([float(row["one_way_range_m"]) for row in rows]),
        "np_duration_s_stats": stats([float(row["np_duration_s"]) for row in rows]),
        "batch_summaries": batch_summaries,
    }


def _svg_bar(parts: list[str], x: float, y: float, width: float, height: float, values: list[tuple[str, int]], title: str, fill: str) -> None:
    parts.append(f'<rect x="{x}" y="{y}" width="{width}" height="{height}" fill="white" stroke="#cbd5e1" stroke-width="1"/>')
    parts.append(f'<text x="{x + 14}" y="{y + 24}" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">{html.escape(title)}</text>')
    if not values:
        return
    max_value = max(value for _, value in values)
    bar_area_x = x + 40
    bar_area_y = y + 46
    bar_area_w = width - 60
    bar_area_h = height - 70
    slot_w = bar_area_w / max(len(values), 1)
    for index, (label, value) in enumerate(values):
        bar_h = 0.0 if max_value == 0 else bar_area_h * value / max_value
        bar_x = bar_area_x + index * slot_w + 6
        bar_y = bar_area_y + bar_area_h - bar_h
        bar_w = max(slot_w - 12, 8)
        parts.append(f'<rect x="{bar_x:.1f}" y="{bar_y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" fill="{fill}"/>')
        parts.append(f'<text x="{bar_x + bar_w/2:.1f}" y="{bar_area_y + bar_area_h + 16:.1f}" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">{html.escape(str(label))}</text>')
        parts.append(f'<text x="{bar_x + bar_w/2:.1f}" y="{bar_y - 4:.1f}" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#0f172a">{value}</text>')


def write_summary_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1500, 980
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: APOLLO real-data ingest scaffold</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Pinned release: APOLLO legacy text normal-point batches a-f from the public archive page, accessed 2026-04-22.</text>',
        '<rect x="48" y="110" width="1404" height="130" fill="white" stroke="#cbd5e1" stroke-width="1"/>',
        f'<text x="68" y="140" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Release summary</text>',
        f'<text x="68" y="170" font-family="monospace" font-size="13" fill="#111827">total normal points = {summary["total_count"]}</text>',
        f'<text x="68" y="194" font-family="monospace" font-size="13" fill="#111827">coverage = {summary["first_timestamp_utc"]} .. {summary["last_timestamp_utc"]}</text>',
        f'<text x="68" y="218" font-family="monospace" font-size="13" fill="#111827">station ids = {summary["station_ids"]}</text>',
        f'<text x="520" y="170" font-family="monospace" font-size="13" fill="#111827">uncertainty_ps median = {summary["uncertainty_ps_stats"]["median"]:.2f}</text>',
        f'<text x="520" y="194" font-family="monospace" font-size="13" fill="#111827">photon_count median = {summary["photon_count_stats"]["median"]:.0f}</text>',
        f'<text x="520" y="218" font-family="monospace" font-size="13" fill="#111827">np_duration_s median = {summary["np_duration_s_stats"]["median"]:.0f}</text>',
        f'<text x="980" y="170" font-family="monospace" font-size="13" fill="#111827">one_way_range_m mean = {summary["one_way_range_m_stats"]["mean"]:.3f}</text>',
        f'<text x="980" y="194" font-family="monospace" font-size="13" fill="#111827">one_way_range_m min/max = {summary["one_way_range_m_stats"]["min"]:.3f} / {summary["one_way_range_m_stats"]["max"]:.3f}</text>',
        f'<text x="980" y="218" font-family="monospace" font-size="13" fill="#111827">retrieval_date = {summary["retrieval_date"]}</text>',
    ]
    _svg_bar(parts, 48, 270, 430, 300, list(summary["reflector_counts"].items()), "Reflector counts", "#2563eb")
    _svg_bar(parts, 534, 270, 430, 300, list(summary["quality_counts"].items()), "Quality grades", "#059669")
    year_items = [(str(k), int(v)) for k, v in summary["year_counts"].items()]
    _svg_bar(parts, 1020, 270, 432, 300, year_items, "Counts by year", "#dc2626")

    parts.append('<rect x="48" y="612" width="1404" height="300" fill="white" stroke="#cbd5e1" stroke-width="1"/>')
    parts.append('<text x="68" y="640" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Batch checks</text>')
    for idx, batch in enumerate(summary["batch_summaries"]):
        y = 674 + idx * 34
        match = "yes" if batch["count_matches_page"] else "NO"
        parts.append(
            f'<text x="68" y="{y}" font-family="monospace" font-size="13" fill="#111827">'
            f'{html.escape(batch["label"])}: parsed={batch["parsed_count"]} expected={batch["expected_count"]} match={match} '
            f'coverage={html.escape(str(batch["first_timestamp_utc"]))}..{html.escape(str(batch["last_timestamp_utc"]))}'
            '</text>'
        )
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def build_release(output_dir: str | Path = ".") -> dict[str, object]:
    root = Path(output_dir)
    release = ApolloRelease()
    retrieval_date = date.today().isoformat()
    release_root = root / "data" / "request4_llr" / release.release_id
    raw_root = release_root / "raw"
    raw_root.mkdir(parents=True, exist_ok=True)

    downloads: list[dict[str, object]] = []
    source_page_payload = fetch_bytes(release.source_page_url)
    source_page_path = raw_root / "norm_pts.html"
    source_page_path.write_bytes(source_page_payload)
    downloads.append(
        {
            "label": "source_page",
            "url": release.source_page_url,
            "path": str(source_page_path),
            "sha256": sha256_bytes(source_page_payload),
            "bytes": len(source_page_payload),
        }
    )

    rows: list[dict[str, object]] = []
    for batch in release.batches:
        payload = fetch_bytes(batch.url)
        filename = Path(batch.url).name
        target = raw_root / filename
        target.write_bytes(payload)
        downloads.append(
            {
                "label": batch.label,
                "url": batch.url,
                "path": str(target),
                "sha256": sha256_bytes(payload),
                "bytes": len(payload),
                "expected_count": batch.expected_count,
            }
        )
        text = payload.decode("utf-8")
        for line in text.splitlines():
            if not line.strip():
                continue
            rows.append(parse_legacy_line(line.rstrip("\n"), batch.label, filename))

    rows.sort(key=lambda item: (item["timestamp_utc"], item["reflector_id"], item["one_way_range_m"]))

    manifest = {
        "release": asdict(release),
        "retrieval_date": retrieval_date,
        "downloads": downloads,
    }
    (release_root / "manifest.json").write_text(json.dumps(manifest, indent=2))
    write_tsv(root / "request4_llr_apollo_normal_points.tsv", rows)
    summary = summarize_rows(rows, release, retrieval_date, downloads)
    (root / "request4_llr_apollo_release_summary.json").write_text(json.dumps(summary, indent=2))
    write_summary_svg(root / "request4_llr_apollo_release_summary.svg", summary)
    return summary


def print_summary(summary: dict[str, object]) -> None:
    print("=== Request 4 APOLLO ingest scaffold ===")
    print(f"total normal points = {summary['total_count']}")
    print(f"coverage = {summary['first_timestamp_utc']} .. {summary['last_timestamp_utc']}")
    print(f"station ids = {summary['station_ids']}")
    print(f"reflector counts = {summary['reflector_counts']}")
    print(f"quality counts = {summary['quality_counts']}")
    print(
        "uncertainty_ps median = "
        f"{summary['uncertainty_ps_stats']['median']:.2f}, "
        f"photon_count median = {summary['photon_count_stats']['median']:.0f}, "
        f"np_duration_s median = {summary['np_duration_s_stats']['median']:.0f}"
    )
    print()
    for batch in summary["batch_summaries"]:
        print(
            f"{batch['label']}: parsed={batch['parsed_count']} expected={batch['expected_count']} "
            f"match={batch['count_matches_page']} "
            f"coverage={batch['first_timestamp_utc']}..{batch['last_timestamp_utc']}"
        )


def main(output_dir: str | Path = ".") -> None:
    summary = build_release(output_dir)
    print_summary(summary)
    print()
    print(f"wrote {Path(output_dir) / 'request4_llr_apollo_normal_points.tsv'}")
    print(f"wrote {Path(output_dir) / 'request4_llr_apollo_release_summary.json'}")
    print(f"wrote {Path(output_dir) / 'request4_llr_apollo_release_summary.svg'}")


if __name__ == "__main__":
    main()
