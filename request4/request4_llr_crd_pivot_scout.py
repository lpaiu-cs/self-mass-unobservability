from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import tarfile
import urllib.request
from dataclasses import asdict, dataclass
from datetime import date, datetime, time, timedelta, timezone
from pathlib import Path
from statistics import mean, median


SPEED_OF_LIGHT_M_S = 299_792_458.0
ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent


@dataclass(frozen=True)
class PivotAsset:
    label: str
    url: str
    local_name: str
    role: str


@dataclass(frozen=True)
class PivotConfig:
    root_dir: str = str(REPO_ROOT / "data" / "request4_llr" / "ilrs_pivot_scout")
    summary_json: str = "request4_llr_crd_pivot_summary.json"
    summary_svg: str = "request4_llr_crd_pivot_summary.svg"
    sample_tsv: str = "request4_llr_crd_lunar_sample_normal_points.tsv"
    manifest_json: str = str(REPO_ROOT / "data" / "request4_llr" / "ilrs_pivot_scout" / "manifest.json")
    assets: tuple[PivotAsset, ...] = (
        PivotAsset(
            label="crd_sample_code_v2.01b",
            url="https://ilrs.gsfc.nasa.gov/docs/2021/crd_v2.01b.tgz",
            local_name="crd_v2.01b.tgz",
            role="ILRS public CRD v2 sample code and converters",
        ),
        PivotAsset(
            label="crd_test_data_v2.01",
            url="https://ilrs.gsfc.nasa.gov/docs/2019/testcrd_for_2.01.tgz",
            local_name="testcrd_for_2.01.tgz",
            role="ILRS public CRD v2 test data",
        ),
        PivotAsset(
            label="orbitNP_release_1.2.1",
            url="https://ilrs.gsfc.nasa.gov/technology/software/orbitNP_release_1.2.1.tar.gz",
            local_name="orbitNP_release_1.2.1.tar.gz",
            role="Public SLR normal-point and residual software",
        ),
        PivotAsset(
            label="MLRS_Lunar_Code_v1.0",
            url="https://ilrs.gsfc.nasa.gov/technology/software/MLRS_Lunar_Code_v1.0.tgz",
            local_name="MLRS_Lunar_Code_v1.0.tgz",
            role="Public lunar prediction, residual, and normal-point code bundle",
        ),
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def head_status(url: str) -> dict[str, object]:
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return {
                "ok": True,
                "status": int(response.status),
                "content_type": response.headers.get("Content-Type"),
                "content_length": response.headers.get("Content-Length"),
            }
    except urllib.error.HTTPError as exc:
        return {"ok": False, "status": int(exc.code), "reason": str(exc.reason)}
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "status": None, "reason": repr(exc)}


def ensure_asset(path: Path, url: str) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=60) as response:
        payload = response.read()
    path.write_bytes(payload)


def ensure_extracted(root: Path, archive_name: str) -> None:
    archive_path = root / archive_name
    if archive_name == "crd_v2.01b.tgz":
        sentinel = root / "crd_sample_code_v2.01b"
    elif archive_name == "testcrd_for_2.01.tgz":
        sentinel = root / "crd_v2.01_test_data"
    elif archive_name == "orbitNP_release_1.2.1.tar.gz":
        sentinel = root / "orbitNP_release_1.2.1"
    elif archive_name == "MLRS_Lunar_Code_v1.0.tgz":
        sentinel = root / "MLRS_Lunar_Code"
    else:
        raise ValueError(f"unhandled archive: {archive_name}")
    if sentinel.exists():
        return
    with tarfile.open(archive_path, "r:gz") as tar:
        tar.extractall(root)


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
                detector_channel = int(tokens[12])
                signal_to_noise = None
                if len(tokens) > 13 and tokens[13].lower() != "na":
                    signal_to_noise = float(tokens[13])
                row = {
                    "source_label": source_label,
                    "source_file": path.name,
                    "crd_format_version": format_version,
                    "station_name": station_name,
                    "station_code": station_code,
                    "target_name": target_name,
                    "target_id": target_id,
                    "timestamp_utc": sec_of_day_to_timestamp(h4["start_date"], h4["start_sod"], sec_of_day),
                    "sec_of_day": sec_of_day,
                    "time_of_flight_s": time_of_flight,
                    "one_way_range_m": 0.5 * SPEED_OF_LIGHT_M_S * time_of_flight,
                    "epoch_event": int(tokens[4]),
                    "np_window_length_s": float(tokens[5]),
                    "num_ranges": int(tokens[6]),
                    "bin_rms_ps": float(tokens[7]),
                    "bin_skew": float(tokens[8]),
                    "bin_kurtosis": float(tokens[9]),
                    "bin_pmm_ps": float(tokens[10]),
                    "return_rate_hz": float(tokens[11]),
                    "detector_channel": detector_channel,
                    "signal_to_noise": signal_to_noise,
                    "pressure_mbar": last_met["pressure_mbar"] if last_met else None,
                    "temperature_k": last_met["temperature_k"] if last_met else None,
                    "humidity_pct": last_met["humidity_pct"] if last_met else None,
                }
                rows.append(row)
    return rows


def write_sample_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1440, 900
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: CRD pivot scout</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Public ILRS software/data availability and LLR CRD sample ingest status.</text>',
        '<rect x="48" y="108" width="1344" height="170" fill="white" stroke="#cbd5e1" stroke-width="1"/>',
        '<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Pivot status</text>',
        f'<text x="66" y="168" font-family="monospace" font-size="13" fill="#111827">CDDIS CRD root status = {summary["cddis_crd_root"]["status"]}</text>',
        f'<text x="66" y="192" font-family="monospace" font-size="13" fill="#111827">CDDIS note = {summary["cddis_crd_root"].get("reason", "OK")}</text>',
        f'<text x="66" y="216" font-family="monospace" font-size="13" fill="#111827">public ILRS assets = {summary["asset_counts"]["public_ok"]} OK / {summary["asset_counts"]["public_total"]} checked</text>',
        f'<text x="66" y="240" font-family="monospace" font-size="13" fill="#111827">public lunar CRD sample points = {summary["lunar_sample"]["row_count"]}</text>',
        f'<text x="66" y="264" font-family="monospace" font-size="13" fill="#111827">lunar targets = {", ".join(summary["lunar_sample"]["targets"])}</text>',
        '<rect x="48" y="320" width="1344" height="230" fill="white" stroke="#cbd5e1" stroke-width="1"/>',
        '<text x="66" y="350" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Asset scan</text>',
    ]
    y = 380
    for asset in summary["assets"]:
        lines.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="12" fill="#111827">{asset["label"]}: status={asset["head"]["status"]} local={asset["local_name"]}</text>'
        )
        y += 22
    lines.extend(
        [
            '<rect x="48" y="590" width="1344" height="240" fill="white" stroke="#cbd5e1" stroke-width="1"/>',
            '<text x="66" y="620" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Interpretation</text>',
        ]
    )
    for idx, note in enumerate(summary["notes"]):
        lines.append(f'<text x="66" y="{652 + idx * 24}" font-family="sans-serif" font-size="13" fill="#374151">{note}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines))


def build_summary(output_dir: str | Path = ".") -> dict[str, object]:
    root = Path(output_dir)
    config = PivotConfig()
    pivot_root = root / config.root_dir
    pivot_root.mkdir(parents=True, exist_ok=True)

    assets_summary = []
    public_ok = 0
    for asset in config.assets:
        local_path = pivot_root / asset.local_name
        ensure_asset(local_path, asset.url)
        ensure_extracted(pivot_root, asset.local_name)
        head = head_status(asset.url)
        if head.get("ok"):
            public_ok += 1
        assets_summary.append(
            {
                "label": asset.label,
                "url": asset.url,
                "local_name": asset.local_name,
                "role": asset.role,
                "head": head,
                "sha256": sha256_file(local_path),
                "size_bytes": local_path.stat().st_size,
            }
        )

    cddis_crd_root = head_status("https://cddis.nasa.gov/archive/slr/data/npt_crd_v2/")

    lunar_files = [
        pivot_root / "MLRS_Lunar_Code/data/analysis/ref/s25y11d016t0231#103.npt",
        pivot_root / "MLRS_Lunar_Code/data/analysis/ref/s25y11d042t0234#103.npt",
    ]
    all_rows: list[dict[str, object]] = []
    for sample_path in lunar_files:
        all_rows.extend(parse_crd_np_file(sample_path, "mlrs_public_lunar_sample"))
    all_rows.sort(key=lambda row: row["timestamp_utc"])
    write_sample_tsv(root / config.sample_tsv, all_rows)

    one_way_ranges = [float(row["one_way_range_m"]) for row in all_rows]
    num_ranges = [int(row["num_ranges"]) for row in all_rows]
    signal_to_noise_values = [row["signal_to_noise"] for row in all_rows if row["signal_to_noise"] is not None]
    targets = sorted({str(row["target_name"]) for row in all_rows})
    stations = sorted({str(row["station_name"]) for row in all_rows})

    manifest = {
        "root_dir": config.root_dir,
        "assets": assets_summary,
        "sample_files": [str(path.relative_to(root)) for path in lunar_files],
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    }
    manifest_path = root / config.manifest_json
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2))

    summary = {
        "config": asdict(config),
        "asset_counts": {"public_ok": public_ok, "public_total": len(config.assets)},
        "cddis_crd_root": cddis_crd_root,
        "assets": assets_summary,
        "lunar_sample": {
            "row_count": len(all_rows),
            "coverage_start": all_rows[0]["timestamp_utc"],
            "coverage_end": all_rows[-1]["timestamp_utc"],
            "targets": targets,
            "stations": stations,
            "one_way_range_m_stats": {
                "min": min(one_way_ranges),
                "max": max(one_way_ranges),
                "mean": mean(one_way_ranges),
                "median": median(one_way_ranges),
            },
            "num_ranges_stats": {
                "min": min(num_ranges),
                "max": max(num_ranges),
                "mean": mean(num_ranges),
                "median": median(num_ranges),
            },
            "signal_to_noise_stats": (
                {
                    "min": min(signal_to_noise_values),
                    "max": max(signal_to_noise_values),
                    "mean": mean(signal_to_noise_values),
                    "median": median(signal_to_noise_values),
                }
                if signal_to_noise_values
                else None
            ),
        },
        "software_scan": {
            "orbitNP": {
                "role": "SLR normal-point generation and residual flattening, not a lunar verdict estimator.",
                "readme_path": f"{config.root_dir}/orbitNP_release_1.2.1/README",
            },
            "mlrs_lunar_code": {
                "role": "Public lunar prediction, residual, and normal-point stack with CRD sample data.",
                "sample_np_path": f"{config.root_dir}/MLRS_Lunar_Code/data/analysis/ref/s25y11d042t0234#103.npt",
                "build_notes": [
                    "bundle contains llr_npt and llr_npt_crd code paths",
                    "makefiles and scripts hard-code g77, /prod/bin, and csh-era assumptions",
                    "gfortran is available locally, so modernization is plausible but not zero-cost",
                ],
            },
            "crd_sample_code": {
                "role": "Public CRD parser/checker/converter starter kit",
                "readme_path": f"{config.root_dir}/crd_sample_code_v2.01b/README_CRD",
            },
        },
        "notes": [
            "CDDIS CRD remains the canonical Request 4 data path, but direct archive access currently returns HTTP 401 and requires Earthdata Login.",
            "Public ILRS tarballs are directly accessible and sufficient to stand up a CRD v2 parser and sample ingest path in the workspace.",
            "The MLRS Lunar Code bundle is the closest public LLR-specific software hand-off point found so far, but it needs local build modernization before it can serve as a drop-in estimator.",
            "orbitNP is useful as a reference for residual flattening and normal-point logic, but it is an SLR tool, not the final weak-field LLR verdict path.",
            "So the healthy next move is a canonical-path bridge: keep the CRD parser/ingest side live and either modernize the public MLRS stack locally or connect to a more mature external estimator stack.",
        ],
    }

    (root / config.summary_json).write_text(json.dumps(summary, indent=2))
    write_svg(root / config.summary_svg, summary)
    return summary


def print_summary(summary: dict[str, object]) -> None:
    print("=== Request 4 CRD pivot scout ===")
    print(f"public ILRS assets OK = {summary['asset_counts']['public_ok']} / {summary['asset_counts']['public_total']}")
    print(f"CDDIS CRD root status = {summary['cddis_crd_root']['status']}")
    print(f"lunar sample rows = {summary['lunar_sample']['row_count']}")
    print(f"coverage = {summary['lunar_sample']['coverage_start']} .. {summary['lunar_sample']['coverage_end']}")
    print(f"targets = {summary['lunar_sample']['targets']}")
    print(f"stations = {summary['lunar_sample']['stations']}")
    print(f"one-way range median = {summary['lunar_sample']['one_way_range_m_stats']['median']:.3f} m")


if __name__ == "__main__":
    summary = build_summary()
    print_summary(summary)
    print()
    print("wrote request4_llr_crd_lunar_sample_normal_points.tsv")
    print("wrote request4_llr_crd_pivot_summary.json")
    print("wrote request4_llr_crd_pivot_summary.svg")
