from __future__ import annotations

import csv
import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SYNODIC_MONTH_DAYS = 29.530588
SIDEREAL_MONTH_DAYS = 27.321661
ANOMALISTIC_MONTH_DAYS = 27.554550
TROPICAL_YEAR_DAYS = 365.242190


@dataclass(frozen=True)
class BaselineScaffoldConfig:
    input_tsv: str = "request4_llr_apollo_normal_points.tsv"
    basis_tsv: str = "request4_llr_apollo_baseline_basis.tsv"
    summary_json: str = "request4_llr_apollo_baseline_scaffold_summary.json"
    summary_svg: str = "request4_llr_apollo_baseline_scaffold_summary.svg"
    reference_batch: str = "group_a"
    reference_reflector: str = "Apollo15"


def parse_timestamp(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00")).astimezone(timezone.utc)


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def centered(values: np.ndarray) -> tuple[np.ndarray, float, float]:
    mean = float(np.mean(values))
    std = float(np.std(values))
    if std == 0.0:
        return values * 0.0, mean, std
    return (values - mean) / std, mean, std


def build_basis(rows: list[dict[str, str]], config: BaselineScaffoldConfig) -> tuple[list[dict[str, object]], dict[str, object]]:
    timestamps = np.array([parse_timestamp(row["timestamp_utc"]).timestamp() for row in rows], dtype=float)
    t0 = float(np.min(timestamps))
    t_days = (timestamps - t0) / 86400.0
    t_center = float(np.mean(t_days))
    t_span = float(np.max(t_days) - np.min(t_days))
    t_centered = t_days - t_center
    t_year = t_centered / TROPICAL_YEAR_DAYS

    pressure = np.array([float(row["pressure_mbar"]) for row in rows], dtype=float)
    temperature = np.array([float(row["temperature_c"]) for row in rows], dtype=float)
    humidity = np.array([float(row["humidity_pct"]) for row in rows], dtype=float)
    duration = np.array([float(row["np_duration_s"]) for row in rows], dtype=float)
    photon_count = np.array([float(row["photon_count"]) for row in rows], dtype=float)
    zero_photon_count = int(np.sum(photon_count <= 0.0))
    # The legacy APOLLO release contains one zero-photon record. Floor at 1
    # before taking the log so the basis remains finite.
    log_photon = np.log10(np.clip(photon_count, 1.0, None))
    uncertainty = np.array([float(row["uncertainty_ps"]) for row in rows], dtype=float)

    pressure_z, pressure_mean, pressure_std = centered(pressure)
    temperature_z, temperature_mean, temperature_std = centered(temperature)
    humidity_z, humidity_mean, humidity_std = centered(humidity)
    duration_z, duration_mean, duration_std = centered(duration)
    log_photon_z, log_photon_mean, log_photon_std = centered(log_photon)
    uncertainty_z, uncertainty_mean, uncertainty_std = centered(uncertainty)

    batch_labels = [row["batch_label"] for row in rows]
    reflector_names = [row["reflector_name"] for row in rows]
    unique_batches = sorted(set(batch_labels))
    unique_reflectors = sorted(set(reflector_names))
    batch_columns = [label for label in unique_batches if label != config.reference_batch]
    reflector_columns = [name for name in unique_reflectors if name != config.reference_reflector]

    synodic_phase = 2.0 * math.pi * t_days / SYNODIC_MONTH_DAYS
    sidereal_phase = 2.0 * math.pi * t_days / SIDEREAL_MONTH_DAYS
    anomalistic_phase = 2.0 * math.pi * t_days / ANOMALISTIC_MONTH_DAYS
    annual_phase = 2.0 * math.pi * t_days / TROPICAL_YEAR_DAYS

    output_rows: list[dict[str, object]] = []
    for idx, row in enumerate(rows):
        basis_row: dict[str, object] = {
            "timestamp_utc": row["timestamp_utc"],
            "batch_label": row["batch_label"],
            "reflector_name": row["reflector_name"],
            "one_way_range_m": float(row["one_way_range_m"]),
            "uncertainty_ps": float(row["uncertainty_ps"]),
            "t_days": float(t_days[idx]),
            "t_year_centered": float(t_year[idx]),
            "t_year_centered_sq": float(t_year[idx] ** 2),
            "pressure_z": float(pressure_z[idx]),
            "temperature_z": float(temperature_z[idx]),
            "humidity_z": float(humidity_z[idx]),
            "np_duration_z": float(duration_z[idx]),
            "log_photon_count_z": float(log_photon_z[idx]),
            "uncertainty_z": float(uncertainty_z[idx]),
            "annual_sin": float(math.sin(annual_phase[idx])),
            "annual_cos": float(math.cos(annual_phase[idx])),
            "synodic_sin": float(math.sin(synodic_phase[idx])),
            "synodic_cos": float(math.cos(synodic_phase[idx])),
            "sidereal_sin": float(math.sin(sidereal_phase[idx])),
            "sidereal_cos": float(math.cos(sidereal_phase[idx])),
            "anomalistic_sin": float(math.sin(anomalistic_phase[idx])),
            "anomalistic_cos": float(math.cos(anomalistic_phase[idx])),
        }
        for label in batch_columns:
            basis_row[f"batch_{label}"] = 1 if row["batch_label"] == label else 0
        for name in reflector_columns:
            basis_row[f"reflector_{name}"] = 1 if row["reflector_name"] == name else 0
        output_rows.append(basis_row)

    summary = {
        "config": asdict(config),
        "row_count": len(rows),
        "coverage_start": rows[0]["timestamp_utc"],
        "coverage_end": rows[-1]["timestamp_utc"],
        "reference_batch": config.reference_batch,
        "reference_reflector": config.reference_reflector,
        "batch_counts": dict(Counter(batch_labels)),
        "reflector_counts": dict(Counter(reflector_names)),
        "basis_columns": list(output_rows[0].keys()),
        "batch_indicator_columns": [f"batch_{label}" for label in batch_columns],
        "reflector_indicator_columns": [f"reflector_{name}" for name in reflector_columns],
        "time_origin_utc": datetime.fromtimestamp(t0, tz=timezone.utc).isoformat().replace("+00:00", "Z"),
        "time_span_days": t_span,
        "zscore_stats": {
            "pressure": {"mean": pressure_mean, "std": pressure_std},
            "temperature": {"mean": temperature_mean, "std": temperature_std},
            "humidity": {"mean": humidity_mean, "std": humidity_std},
            "np_duration": {"mean": duration_mean, "std": duration_std},
            "log_photon_count": {"mean": log_photon_mean, "std": log_photon_std},
            "uncertainty_ps": {"mean": uncertainty_mean, "std": uncertainty_std},
        },
        "zero_photon_count": zero_photon_count,
        "harmonic_period_days": {
            "annual": TROPICAL_YEAR_DAYS,
            "synodic_month": SYNODIC_MONTH_DAYS,
            "sidereal_month": SIDEREAL_MONTH_DAYS,
            "anomalistic_month": ANOMALISTIC_MONTH_DAYS,
        },
        "notes": [
            "This is a nuisance/design scaffold for the APOLLO-only branch, not a GR residual fit.",
            "Batch indicators are included because APOLLO release groups may coincide with reduction-version boundaries.",
            "Reflector indicators provide a first scaffold for reflector-specific biases.",
            "Annual and lunar-month harmonics are surrogate basis functions only; they are not a physical Earth-Moon ephemeris.",
            "One zero-photon record is floored at photon_count=1 before taking log10 so the design matrix stays finite.",
        ],
        "stop_rule": (
            "If a narrow APOLLO-only baseline GR residual layer cannot be stabilized on top of this scaffold, "
            "the next move should be to attach to an existing LLR estimator or pivot to the CRD/ILRS path."
        ),
    }
    return output_rows, summary


def write_basis_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _svg_bar(parts: list[str], x: float, y: float, width: float, height: float, values: list[tuple[str, int]], title: str, fill: str) -> None:
    parts.append(f'<rect x="{x}" y="{y}" width="{width}" height="{height}" fill="white" stroke="#cbd5e1" stroke-width="1"/>')
    parts.append(f'<text x="{x + 12}" y="{y + 24}" font-family="sans-serif" font-size="15" font-weight="700" fill="#111827">{title}</text>')
    if not values:
        return
    max_value = max(value for _, value in values)
    inner_x = x + 36
    inner_y = y + 42
    inner_w = width - 52
    inner_h = height - 66
    slot_w = inner_w / len(values)
    for idx, (label, value) in enumerate(values):
        bar_h = 0.0 if max_value == 0 else inner_h * value / max_value
        bar_x = inner_x + idx * slot_w + 5
        bar_y = inner_y + inner_h - bar_h
        bar_w = max(slot_w - 10, 8)
        parts.append(f'<rect x="{bar_x:.1f}" y="{bar_y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" fill="{fill}"/>')
        parts.append(f'<text x="{bar_x + bar_w/2:.1f}" y="{inner_y + inner_h + 14:.1f}" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">{label}</text>')
        parts.append(f'<text x="{bar_x + bar_w/2:.1f}" y="{bar_y - 4:.1f}" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#0f172a">{value}</text>')


def write_summary_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1440, 900
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: APOLLO baseline nuisance scaffold</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">APOLLO-only design matrix scaffold with batch, reflector, meteo, and lunar/annual surrogate basis terms.</text>',
        '<rect x="48" y="108" width="1344" height="132" fill="white" stroke="#cbd5e1" stroke-width="1"/>',
        '<text x="68" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Scaffold summary</text>',
        f'<text x="68" y="168" font-family="monospace" font-size="13" fill="#111827">row_count = {summary["row_count"]}</text>',
        f'<text x="68" y="192" font-family="monospace" font-size="13" fill="#111827">coverage = {summary["coverage_start"]} .. {summary["coverage_end"]}</text>',
        f'<text x="68" y="216" font-family="monospace" font-size="13" fill="#111827">reference batch = {summary["reference_batch"]}; reference reflector = {summary["reference_reflector"]}</text>',
        f'<text x="760" y="168" font-family="monospace" font-size="13" fill="#111827">batch columns = {len(summary["batch_indicator_columns"])}</text>',
        f'<text x="760" y="192" font-family="monospace" font-size="13" fill="#111827">reflector columns = {len(summary["reflector_indicator_columns"])}</text>',
        f'<text x="760" y="216" font-family="monospace" font-size="13" fill="#111827">time span days = {summary["time_span_days"]:.2f}</text>',
    ]
    _svg_bar(parts, 48, 280, 420, 260, list(summary["batch_counts"].items()), "Batch counts", "#2563eb")
    _svg_bar(parts, 510, 280, 420, 260, list(summary["reflector_counts"].items()), "Reflector counts", "#059669")
    harm = [(k, int(v)) for k, v in {"annual": 2, "synodic": 2, "sidereal": 2, "anomalistic": 2}.items()]
    _svg_bar(parts, 972, 280, 420, 260, harm, "Harmonic blocks", "#dc2626")
    parts.append('<rect x="48" y="586" width="1344" height="250" fill="white" stroke="#cbd5e1" stroke-width="1"/>')
    parts.append('<text x="68" y="614" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Notes</text>')
    for idx, note in enumerate(summary["notes"]):
        parts.append(f'<text x="68" y="{646 + idx*24}" font-family="sans-serif" font-size="13" fill="#374151">{note}</text>')
    parts.append(f'<text x="68" y="758" font-family="sans-serif" font-size="13" fill="#111827">stop rule: {summary["stop_rule"]}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def build_scaffold(output_dir: str | Path = ".") -> dict[str, object]:
    root = Path(output_dir)
    config = BaselineScaffoldConfig()
    input_path = root / config.input_tsv
    rows = load_rows(input_path)
    rows.sort(key=lambda row: row["timestamp_utc"])
    basis_rows, summary = build_basis(rows, config)
    write_basis_tsv(root / config.basis_tsv, basis_rows)
    (root / config.summary_json).write_text(json.dumps(summary, indent=2))
    write_summary_svg(root / config.summary_svg, summary)
    return summary


def print_summary(summary: dict[str, object]) -> None:
    print("=== Request 4 APOLLO baseline nuisance scaffold ===")
    print(f"row_count = {summary['row_count']}")
    print(f"coverage = {summary['coverage_start']} .. {summary['coverage_end']}")
    print(f"reference batch = {summary['reference_batch']}")
    print(f"reference reflector = {summary['reference_reflector']}")
    print(f"batch columns = {len(summary['batch_indicator_columns'])}")
    print(f"reflector columns = {len(summary['reflector_indicator_columns'])}")
    print(f"time span days = {summary['time_span_days']:.2f}")
    print(f"batch counts = {summary['batch_counts']}")
    print(f"reflector counts = {summary['reflector_counts']}")


def main(output_dir: str | Path = ".") -> None:
    summary = build_scaffold(output_dir)
    print_summary(summary)
    root = Path(output_dir)
    print()
    print(f"wrote {root / BaselineScaffoldConfig().basis_tsv}")
    print(f"wrote {root / BaselineScaffoldConfig().summary_json}")
    print(f"wrote {root / BaselineScaffoldConfig().summary_svg}")


if __name__ == "__main__":
    main()
