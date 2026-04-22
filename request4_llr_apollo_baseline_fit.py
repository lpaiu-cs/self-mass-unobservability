from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from skyfield.api import Loader


SPEED_OF_LIGHT_M_S = 299_792_458.0
TROPICAL_YEAR_DAYS = 365.242190
SYNODIC_MONTH_DAYS = 29.530588
SIDEREAL_MONTH_DAYS = 27.321661
ANOMALISTIC_MONTH_DAYS = 27.554550
SOLAR_DAY_DAYS = 1.0
SIDEREAL_DAY_DAYS = 0.9972695663


@dataclass(frozen=True)
class BaselineFitConfig:
    basis_tsv: str = "request4_llr_apollo_baseline_basis.tsv"
    residuals_tsv: str = "request4_llr_apollo_baseline_residuals.tsv"
    coefficients_tsv: str = "request4_llr_apollo_baseline_coefficients.tsv"
    summary_json: str = "request4_llr_apollo_baseline_fit_summary.json"
    summary_svg: str = "request4_llr_apollo_baseline_fit_summary.svg"
    skyfield_cache_dir: str = "data/request4_llr/skyfield_cache"
    ephemeris_name: str = "de421.bsp"
    annual_harmonics: int = 3
    monthly_harmonics: int = 6
    day_harmonics: int = 4


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def uncertainty_ps_to_one_way_sigma_m(values_ps: np.ndarray) -> np.ndarray:
    return 0.5 * SPEED_OF_LIGHT_M_S * values_ps * 1.0e-12


def parse_timestamp_columns(rows: list[dict[str, str]]) -> tuple[np.ndarray, ...]:
    years = np.array([int(row["timestamp_utc"][0:4]) for row in rows], dtype=int)
    months = np.array([int(row["timestamp_utc"][5:7]) for row in rows], dtype=int)
    days = np.array([int(row["timestamp_utc"][8:10]) for row in rows], dtype=int)
    hours = np.array([int(row["timestamp_utc"][11:13]) for row in rows], dtype=int)
    minutes = np.array([int(row["timestamp_utc"][14:16]) for row in rows], dtype=int)
    seconds = np.array([float(row["timestamp_utc"][17:-1]) for row in rows], dtype=float)
    return years, months, days, hours, minutes, seconds


def build_nominal_geocentric_range(rows: list[dict[str, str]], root: Path, config: BaselineFitConfig) -> np.ndarray:
    loader = Loader(str(root / config.skyfield_cache_dir))
    eph = loader(config.ephemeris_name)
    ts = loader.timescale()
    years, months, days, hours, minutes, seconds = parse_timestamp_columns(rows)
    t = ts.utc(years, months, days, hours, minutes, seconds)
    earth = eph["earth"]
    moon = eph["moon"]
    return (moon.at(t) - earth.at(t)).distance().m


def collect_indicator_columns(rows: list[dict[str, str]], prefix: str, exclude: str) -> list[str]:
    names = []
    for key in rows[0].keys():
        if key.startswith(prefix) and key != exclude:
            names.append(key)
    return names


def build_design_matrix(rows: list[dict[str, str]], config: BaselineFitConfig) -> tuple[np.ndarray, list[str]]:
    t_days = np.array([float(row["t_days"]) for row in rows], dtype=float)
    t_year = np.array([float(row["t_year_centered"]) for row in rows], dtype=float)
    names: list[str] = ["const", "t_year_centered", "t_year_centered_sq"]
    columns = [np.ones(len(rows)), t_year, t_year**2]

    periodic_blocks = [
        ("annual", TROPICAL_YEAR_DAYS, config.annual_harmonics),
        ("synodic", SYNODIC_MONTH_DAYS, config.monthly_harmonics),
        ("sidereal_month", SIDEREAL_MONTH_DAYS, config.monthly_harmonics),
        ("anomalistic", ANOMALISTIC_MONTH_DAYS, config.monthly_harmonics),
        ("solar_day", SOLAR_DAY_DAYS, config.day_harmonics),
        ("sidereal_day", SIDEREAL_DAY_DAYS, config.day_harmonics),
    ]
    for label, period_days, max_harmonic in periodic_blocks:
        for harmonic in range(1, max_harmonic + 1):
            phase = 2.0 * math.pi * harmonic * t_days / period_days
            columns.append(np.sin(phase))
            columns.append(np.cos(phase))
            names.append(f"{label}_{harmonic}_sin")
            names.append(f"{label}_{harmonic}_cos")

    annual_frequency = 2.0 * math.pi / TROPICAL_YEAR_DAYS
    monthly_bases = [
        ("synodic", 2.0 * math.pi / SYNODIC_MONTH_DAYS),
        ("sidereal_month", 2.0 * math.pi / SIDEREAL_MONTH_DAYS),
        ("anomalistic", 2.0 * math.pi / ANOMALISTIC_MONTH_DAYS),
    ]
    for label, base_frequency in monthly_bases:
        for harmonic in range(1, 4):
            for suffix, frequency in (
                ("minus_annual", harmonic * base_frequency - annual_frequency),
                ("plus_annual", harmonic * base_frequency + annual_frequency),
            ):
                phase = frequency * t_days
                columns.append(np.sin(phase))
                columns.append(np.cos(phase))
                names.append(f"{label}_{harmonic}_{suffix}_sin")
                names.append(f"{label}_{harmonic}_{suffix}_cos")

    for key in (
        "pressure_z",
        "temperature_z",
        "humidity_z",
        "np_duration_z",
        "log_photon_count_z",
        "uncertainty_z",
    ):
        columns.append(np.array([float(row[key]) for row in rows], dtype=float))
        names.append(key)

    batch_columns = collect_indicator_columns(rows, "batch_", "batch_label")
    reflector_columns = collect_indicator_columns(rows, "reflector_", "reflector_name")
    for key in batch_columns + reflector_columns:
        columns.append(np.array([float(row[key]) for row in rows], dtype=float))
        names.append(key)

    return np.column_stack(columns), names


def weighted_lstsq(X: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> tuple[np.ndarray, dict[str, float | int]]:
    safe_sigma = np.clip(sigma, 1.0e-6, None)
    inv_sigma = 1.0 / safe_sigma
    Xw = X * inv_sigma[:, None]
    yw = y * inv_sigma
    beta, residuals, rank, singular_values = np.linalg.lstsq(Xw, yw, rcond=1.0e-12)
    diagnostics = {
        "rank": int(rank),
        "singular_value_max": float(np.max(singular_values)),
        "singular_value_min": float(np.min(singular_values)),
        "condition_number": float(np.max(singular_values) / np.min(singular_values)),
        "weighted_residual_sum_sq": float(np.sum(residuals)) if len(residuals) else 0.0,
    }
    return beta, diagnostics


def residual_stats(residual_m: np.ndarray, sigma_m: np.ndarray) -> dict[str, float]:
    weights = 1.0 / np.clip(sigma_m, 1.0e-6, None) ** 2
    return {
        "rms_m": float(np.sqrt(np.mean(residual_m**2))),
        "wrms_m": float(np.sqrt(np.sum(weights * residual_m**2) / np.sum(weights))),
        "median_abs_m": float(np.median(np.abs(residual_m))),
        "p95_abs_m": float(np.percentile(np.abs(residual_m), 95.0)),
        "max_abs_m": float(np.max(np.abs(residual_m))),
    }


def grouped_residual_stats(
    residual_m: np.ndarray,
    sigma_m: np.ndarray,
    labels: list[str],
) -> dict[str, dict[str, float | int]]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for idx, label in enumerate(labels):
        grouped[label].append(idx)

    output: dict[str, dict[str, float | int]] = {}
    for label, indices in grouped.items():
        idx = np.array(indices, dtype=int)
        subset_residual = residual_m[idx]
        subset_sigma = sigma_m[idx]
        weights = 1.0 / np.clip(subset_sigma, 1.0e-6, None) ** 2
        output[label] = {
            "count": int(len(indices)),
            "weighted_mean_m": float(np.sum(weights * subset_residual) / np.sum(weights)),
            "wrms_m": float(np.sqrt(np.sum(weights * subset_residual**2) / np.sum(weights))),
            "median_abs_m": float(np.median(np.abs(subset_residual))),
        }
    return output


def write_residuals_tsv(
    path: Path,
    rows: list[dict[str, str]],
    nominal_range_m: np.ndarray,
    nominal_residual_m: np.ndarray,
    surrogate_correction_m: np.ndarray,
    fitted_residual_m: np.ndarray,
    sigma_m: np.ndarray,
) -> None:
    fieldnames = [
        "timestamp_utc",
        "batch_label",
        "reflector_name",
        "one_way_range_m",
        "nominal_geocentric_range_m",
        "nominal_residual_m",
        "surrogate_correction_m",
        "fitted_residual_m",
        "uncertainty_ps",
        "one_way_sigma_m",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for idx, row in enumerate(rows):
            writer.writerow(
                {
                    "timestamp_utc": row["timestamp_utc"],
                    "batch_label": row["batch_label"],
                    "reflector_name": row["reflector_name"],
                    "one_way_range_m": row["one_way_range_m"],
                    "nominal_geocentric_range_m": nominal_range_m[idx],
                    "nominal_residual_m": nominal_residual_m[idx],
                    "surrogate_correction_m": surrogate_correction_m[idx],
                    "fitted_residual_m": fitted_residual_m[idx],
                    "uncertainty_ps": row["uncertainty_ps"],
                    "one_way_sigma_m": sigma_m[idx],
                }
            )


def write_coefficients_tsv(path: Path, names: list[str], beta: np.ndarray) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["term", "coefficient_m"], delimiter="\t")
        writer.writeheader()
        for name, value in zip(names, beta, strict=True):
            writer.writerow({"term": name, "coefficient_m": float(value)})


def write_summary_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1440, 920
    nominal = summary["nominal_residual_stats"]
    fitted = summary["fitted_residual_stats"]
    batch_stats = summary["batch_stats"]
    top_batches = sorted(batch_stats.items(), key=lambda item: item[1]["wrms_m"], reverse=True)[:6]

    def block(parts: list[str], x: int, y: int, w: int, h: int, title: str, lines: list[str]) -> None:
        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" stroke="#cbd5e1" stroke-width="1"/>')
        parts.append(f'<text x="{x + 18}" y="{y + 28}" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">{title}</text>')
        for idx, line in enumerate(lines):
            parts.append(
                f'<text x="{x + 18}" y="{y + 56 + idx * 22}" font-family="monospace" font-size="13" fill="#111827">{line}</text>'
            )

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 4: APOLLO baseline surrogate fit</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">DE421 Earth-Moon center nominal plus APOLLO-only batch/reflector/meteo/harmonic residual layer.</text>',
    ]
    block(
        parts,
        48,
        110,
        420,
        170,
        "Nominal residual",
        [
            f"rms_m = {nominal['rms_m']:.3f}",
            f"wrms_m = {nominal['wrms_m']:.3f}",
            f"median_abs_m = {nominal['median_abs_m']:.3f}",
            f"p95_abs_m = {nominal['p95_abs_m']:.3f}",
            f"max_abs_m = {nominal['max_abs_m']:.3f}",
        ],
    )
    block(
        parts,
        510,
        110,
        420,
        170,
        "Fitted residual",
        [
            f"rms_m = {fitted['rms_m']:.3f}",
            f"wrms_m = {fitted['wrms_m']:.3f}",
            f"median_abs_m = {fitted['median_abs_m']:.3f}",
            f"p95_abs_m = {fitted['p95_abs_m']:.3f}",
            f"max_abs_m = {fitted['max_abs_m']:.3f}",
        ],
    )
    block(
        parts,
        972,
        110,
        420,
        170,
        "Assessment",
        [
            f"status = {summary['assessment']['status']}",
            f"improvement_rms = {summary['assessment']['improvement_factor_rms']:.3f}x",
            f"improvement_wrms = {summary['assessment']['improvement_factor_wrms']:.3f}x",
            f"design_cols = {summary['design_matrix']['column_count']}",
            f"condition_number = {summary['fit_diagnostics']['condition_number']:.3e}",
        ],
    )
    parts.append('<rect x="48" y="330" width="1344" height="220" fill="white" stroke="#cbd5e1" stroke-width="1"/>')
    parts.append('<text x="66" y="360" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Worst batch residual blocks</text>')
    for idx, (label, stats) in enumerate(top_batches):
        x = 66 + idx * 220
        parts.append(f'<text x="{x}" y="390" font-family="sans-serif" font-size="13" font-weight="700" fill="#111827">{label}</text>')
        parts.append(f'<text x="{x}" y="414" font-family="monospace" font-size="12" fill="#111827">count = {stats["count"]}</text>')
        parts.append(f'<text x="{x}" y="436" font-family="monospace" font-size="12" fill="#111827">wrms_m = {stats["wrms_m"]:.1f}</text>')
        parts.append(f'<text x="{x}" y="458" font-family="monospace" font-size="12" fill="#111827">median_abs_m = {stats["median_abs_m"]:.1f}</text>')
        parts.append(f'<text x="{x}" y="480" font-family="monospace" font-size="12" fill="#111827">wmean_m = {stats["weighted_mean_m"]:.1f}</text>')
    parts.append('<rect x="48" y="590" width="1344" height="270" fill="white" stroke="#cbd5e1" stroke-width="1"/>')
    parts.append('<text x="66" y="620" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Notes</text>')
    for idx, note in enumerate(summary["notes"]):
        parts.append(f'<text x="66" y="{652 + idx * 24}" font-family="sans-serif" font-size="13" fill="#374151">{note}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def build_baseline_fit(output_dir: str | Path = ".") -> dict[str, object]:
    root = Path(output_dir)
    config = BaselineFitConfig()
    rows = load_rows(root / config.basis_tsv)
    nominal_range_m = build_nominal_geocentric_range(rows, root, config)
    observed_range_m = np.array([float(row["one_way_range_m"]) for row in rows], dtype=float)
    uncertainty_ps = np.array([float(row["uncertainty_ps"]) for row in rows], dtype=float)
    sigma_m = uncertainty_ps_to_one_way_sigma_m(uncertainty_ps)
    nominal_residual_m = observed_range_m - nominal_range_m

    X, names = build_design_matrix(rows, config)
    beta, diagnostics = weighted_lstsq(X, nominal_residual_m, sigma_m)
    surrogate_correction_m = X @ beta
    fitted_residual_m = nominal_residual_m - surrogate_correction_m

    nominal_stats = residual_stats(nominal_residual_m, sigma_m)
    fitted_stats = residual_stats(fitted_residual_m, sigma_m)
    batch_stats = grouped_residual_stats(fitted_residual_m, sigma_m, [row["batch_label"] for row in rows])
    reflector_stats = grouped_residual_stats(fitted_residual_m, sigma_m, [row["reflector_name"] for row in rows])

    assessment_status = "not_closed"
    if fitted_stats["median_abs_m"] < 1.0e4 and fitted_stats["wrms_m"] < 1.0e4:
        assessment_status = "narrow_surrogate_closed"
    elif fitted_stats["median_abs_m"] < 1.0e5 and fitted_stats["wrms_m"] < 1.0e5:
        assessment_status = "partially_closed"

    summary = {
        "config": asdict(config),
        "row_count": len(rows),
        "coverage_start": rows[0]["timestamp_utc"],
        "coverage_end": rows[-1]["timestamp_utc"],
        "nominal_model": {
            "description": "Skyfield DE421 Earth-Moon center distance only.",
            "ephemeris_name": config.ephemeris_name,
            "cache_dir": config.skyfield_cache_dir,
        },
        "design_matrix": {
            "column_count": int(X.shape[1]),
            "terms": names,
        },
        "fit_diagnostics": diagnostics,
        "nominal_residual_stats": nominal_stats,
        "fitted_residual_stats": fitted_stats,
        "batch_stats": batch_stats,
        "reflector_stats": reflector_stats,
        "assessment": {
            "status": assessment_status,
            "improvement_factor_rms": nominal_stats["rms_m"] / fitted_stats["rms_m"],
            "improvement_factor_wrms": nominal_stats["wrms_m"] / fitted_stats["wrms_m"],
        },
        "notes": [
            "This is a narrow APOLLO-only baseline surrogate, not a production LLR estimator.",
            "The nominal model is geocentric Earth-Moon center distance from DE421 and omits full station, reflector, troposphere, orientation, and relativistic light-time modeling.",
            "The residual layer includes batch offsets, reflector biases, meteo surrogates, and annual/monthly/day harmonic blocks only.",
            "Residuals remaining at O(10^5-10^6) m mean the bespoke APOLLO-only branch is still far from a credible weak-field verdict path.",
            "At that point the healthy next move is to attach to an existing LLR estimator or pivot to the CRD/ILRS canonical path rather than keep growing this surrogate indefinitely.",
        ],
    }

    write_residuals_tsv(
        root / config.residuals_tsv,
        rows,
        nominal_range_m,
        nominal_residual_m,
        surrogate_correction_m,
        fitted_residual_m,
        sigma_m,
    )
    write_coefficients_tsv(root / config.coefficients_tsv, names, beta)
    (root / config.summary_json).write_text(json.dumps(summary, indent=2))
    write_summary_svg(root / config.summary_svg, summary)
    return summary


def print_summary(summary: dict[str, object]) -> None:
    nominal = summary["nominal_residual_stats"]
    fitted = summary["fitted_residual_stats"]
    assessment = summary["assessment"]
    print("=== Request 4 APOLLO baseline surrogate fit ===")
    print(f"row_count = {summary['row_count']}")
    print(f"coverage = {summary['coverage_start']} .. {summary['coverage_end']}")
    print(f"nominal rms_m = {nominal['rms_m']:.3f}")
    print(f"nominal wrms_m = {nominal['wrms_m']:.3f}")
    print(f"fitted rms_m = {fitted['rms_m']:.3f}")
    print(f"fitted wrms_m = {fitted['wrms_m']:.3f}")
    print(f"fitted median_abs_m = {fitted['median_abs_m']:.3f}")
    print(f"improvement factor (rms) = {assessment['improvement_factor_rms']:.3f}x")
    print(f"assessment = {assessment['status']}")


if __name__ == "__main__":
    summary = build_baseline_fit()
    print_summary(summary)
    print()
    print("wrote request4_llr_apollo_baseline_residuals.tsv")
    print("wrote request4_llr_apollo_baseline_coefficients.tsv")
    print("wrote request4_llr_apollo_baseline_fit_summary.json")
    print("wrote request4_llr_apollo_baseline_fit_summary.svg")
