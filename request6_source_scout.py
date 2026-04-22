from __future__ import annotations

import itertools
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

import request6_clock_sector as r6
import request6_lever_arm_audit as audit


@dataclass(frozen=True)
class SourceCandidate:
    name: str
    system_type: str
    readiness: str
    nongamma_strategy: str
    mp_msun: float
    mp_sigma_msun: float
    mc_msun: float
    mc_sigma_msun: float
    pb_days: float
    eccentricity: float
    gamma_obs_s: float
    gamma_sigma_s: float
    source_urls: tuple[str, ...]
    notes: str


CANDIDATES: tuple[SourceCandidate, ...] = (
    SourceCandidate(
        name="PSR B1913+16",
        system_type="DNS",
        readiness="moderate",
        nongamma_strategy="omdot + Pbdot (radiative; Galactic/kinematic corrections required)",
        mp_msun=1.4398,
        mp_sigma_msun=0.0002,
        mc_msun=1.3886,
        mc_sigma_msun=0.0002,
        pb_days=0.322997448911,
        eccentricity=0.6171334,
        gamma_obs_s=4.2992e-3,
        gamma_sigma_s=8.0e-7,
        source_urls=("https://arxiv.org/abs/1011.0718",),
        notes="Best precision anchor in the scanned sample. It does not bracket s* by itself, but its gamma precision is already at the audit target scale.",
    ),
    SourceCandidate(
        name="PSR J1757-1854",
        system_type="DNS",
        readiness="clean",
        nongamma_strategy="omdot + h3 + varsigma (Shapiro-based); Pbdot available as a cross-check",
        mp_msun=1.3384,
        mp_sigma_msun=0.0009,
        mc_msun=1.3946,
        mc_sigma_msun=0.0009,
        pb_days=0.18353783587,
        eccentricity=0.6058142,
        gamma_obs_s=3.587e-3,
        gamma_sigma_s=1.2e-5,
        source_urls=("https://arxiv.org/abs/1711.07697",),
        notes="Methodologically attractive because the leave-one-out mass route can be kept almost entirely inside non-gamma PK information.",
    ),
    SourceCandidate(
        name="PSR J1756-2251",
        system_type="DNS",
        readiness="clean",
        nongamma_strategy="omdot + r + s (Shapiro-based); Pbdot measured but systematics discussed in the paper",
        mp_msun=1.341,
        mp_sigma_msun=0.007,
        mc_msun=1.230,
        mc_sigma_msun=0.007,
        pb_days=0.31963390143,
        eccentricity=0.1805694,
        gamma_obs_s=1.148e-3,
        gamma_sigma_s=9.0e-6,
        source_urls=("https://arxiv.org/abs/1406.5507",),
        notes="Useful clean Shapiro cross-check. The source table labels gamma in ms but the numerical value is on the 10^-3 s scale; the scout uses the physically consistent seconds interpretation.",
    ),
    SourceCandidate(
        name="PSR J1913+1102",
        system_type="DNS",
        readiness="moderate",
        nongamma_strategy="omdot + Pbdot (radiative); gamma then used only as the clock observable",
        mp_msun=1.62,
        mp_sigma_msun=0.03,
        mc_msun=1.27,
        mc_sigma_msun=0.03,
        pb_days=0.2062523345,
        eccentricity=0.089531,
        gamma_obs_s=4.71e-4,
        gamma_sigma_s=1.5e-5,
        source_urls=("https://arxiv.org/abs/2007.04175",),
        notes="The best high-s lever-arm source in the scanned set, but the gamma precision is far weaker than B1913+16.",
    ),
    SourceCandidate(
        name="PSR J1141-6545",
        system_type="NS-WD",
        readiness="conditional",
        nongamma_strategy="omdot + Pbdot + Shapiro-like geometry, but xdot / spin-orbit coupling must be modelled",
        mp_msun=1.27,
        mp_sigma_msun=0.01,
        mc_msun=1.02,
        mc_sigma_msun=0.01,
        pb_days=0.1976509593,
        eccentricity=0.171884,
        gamma_obs_s=7.73e-4,
        gamma_sigma_s=1.1e-5,
        source_urls=("https://arxiv.org/abs/0804.0956", "https://arxiv.org/abs/2001.11405"),
        notes="Best low-s candidate found in the scan. It gives the desired opposite-side lever arm, but only if the WD-driven xdot nuisance is carried explicitly.",
    ),
    SourceCandidate(
        name="PSR J1906+0746",
        system_type="compact companion (DNS or WD)",
        readiness="tricky",
        nongamma_strategy="omdot + Pbdot + possible xdot, but gamma and xdot are correlated",
        mp_msun=1.316,
        mp_sigma_msun=0.005,
        mc_msun=1.297,
        mc_sigma_msun=0.005,
        pb_days=0.165993045,
        eccentricity=0.0853022,
        gamma_obs_s=4.59e-4,
        gamma_sigma_s=2.0e-6,
        source_urls=("https://arxiv.org/abs/2602.05947",),
        notes="Interesting because it sits slightly below s*, but the paper explicitly reports a ~3.5 sigma mass shift once xdot is fitted, so this is not a front-line clean test case.",
    ),
)


def load_request6_summary(path: str | Path = "request6_clock_sector_summary.json") -> dict[str, object]:
    return json.loads(Path(path).read_text())


def load_lever_arm_audit(path: str | Path = "request6_lever_arm_audit_summary.json") -> dict[str, object]:
    return json.loads(Path(path).read_text())


def mean_s_from_mass(mp_msun: float) -> tuple[float, float, float]:
    lower, upper = r6.mass_conditioned_s_prior(np.array([mp_msun], dtype=float))
    s_low = float(lower[0])
    s_high = float(upper[0])
    return s_low, s_high, 0.5 * (s_low + s_high)


def readiness_penalty(readiness: str) -> float:
    return {
        "clean": 1.0,
        "moderate": 0.8,
        "conditional": 0.45,
        "tricky": 0.2,
    }[readiness]


def candidate_row(candidate: SourceCandidate, sstar: float) -> dict[str, object]:
    s_low, s_high, s_mid = mean_s_from_mass(candidate.mp_msun)
    xc = candidate.mc_msun / (candidate.mp_msun + candidate.mc_msun)
    sigma_delta = candidate.gamma_sigma_s / candidate.gamma_obs_s
    ds = s_mid - sstar
    return {
        "name": candidate.name,
        "readiness": candidate.readiness,
        "system_type": candidate.system_type,
        "nongamma_strategy": candidate.nongamma_strategy,
        "mp_msun": candidate.mp_msun,
        "mp_sigma_msun": candidate.mp_sigma_msun,
        "mc_msun": candidate.mc_msun,
        "mc_sigma_msun": candidate.mc_sigma_msun,
        "pb_days": candidate.pb_days,
        "eccentricity": candidate.eccentricity,
        "gamma_obs_s": candidate.gamma_obs_s,
        "gamma_sigma_s": candidate.gamma_sigma_s,
        "gamma_fractional_sigma": sigma_delta,
        "xc": xc,
        "s_low": s_low,
        "s_high": s_high,
        "s_mid": s_mid,
        "delta_s_from_sstar": ds,
        "direction": "positive" if ds > 0.0 else "negative" if ds < 0.0 else "zero",
        "source_urls": list(candidate.source_urls),
        "notes": candidate.notes,
    }


def scout_metrics(
    row: dict[str, object],
    base_rows: list[dict[str, float | str]],
    sstar: float,
    current_kappa95: float,
) -> dict[str, object]:
    ds = float(row["delta_s_from_sstar"])
    xc = float(row["xc"])
    sigma_delta = float(row["gamma_fractional_sigma"])
    slope_proxy = (ds**2) / ((1.0 + xc) ** 2 * sigma_delta**2) if sigma_delta > 0.0 else math.inf
    rows = list(base_rows) + [
        {
            "name": str(row["name"]),
            "sbar": float(row["s_mid"]),
            "xc": xc,
            "gamma_gr_s": float(row["gamma_obs_s"]),
            "sigma_gamma_s": float(row["gamma_sigma_s"]),
            "sigma_delta": sigma_delta,
            "delta_gamma_obs": 0.0,
        }
    ]
    fisher = audit.fisher_eta_kappa(rows, sstar)
    new_kappa95 = float(fisher["abs95_kappa_star"])
    improvement = current_kappa95 / new_kappa95 if new_kappa95 > 0.0 else math.inf
    composite = slope_proxy * readiness_penalty(str(row["readiness"]))
    return {
        "slope_proxy": slope_proxy,
        "new_abs95_kappa_star": new_kappa95,
        "improvement_factor": improvement,
        "composite_priority_score": composite,
    }


def _combo_entry(
    combo: tuple[dict[str, object], ...],
    base_rows: list[dict[str, float | str]],
    sstar: float,
    current_kappa95: float,
) -> dict[str, object]:
    rows = list(base_rows)
    for item in combo:
        rows.append(
            {
                "name": str(item["name"]),
                "sbar": float(item["s_mid"]),
                "xc": float(item["xc"]),
                "gamma_gr_s": float(item["gamma_obs_s"]),
                "sigma_gamma_s": float(item["gamma_sigma_s"]),
                "sigma_delta": float(item["gamma_fractional_sigma"]),
                "delta_gamma_obs": 0.0,
            }
        )
    fit = audit.fisher_eta_kappa(rows, sstar)
    kappa95 = float(fit["abs95_kappa_star"])
    return {
        "names": [str(item["name"]) for item in combo],
        "new_abs95_kappa_star": kappa95,
        "improvement_factor": current_kappa95 / kappa95 if kappa95 > 0.0 else math.inf,
        "directions": [str(item["direction"]) for item in combo],
        "readiness": [str(item["readiness"]) for item in combo],
    }


def best_pair_summaries(
    rows: list[dict[str, object]],
    base_rows: list[dict[str, float | str]],
    sstar: float,
    current_kappa95: float,
) -> dict[str, object]:
    entries = [_combo_entry(combo, base_rows, sstar, current_kappa95) for combo in itertools.combinations(rows, 2)]
    entries.sort(key=lambda item: item["new_abs95_kappa_star"])
    cross_side = [
        entry
        for entry in entries
        if "positive" in entry["directions"] and "negative" in entry["directions"]
    ]
    non_tricky_cross_side = [
        entry
        for entry in cross_side
        if "tricky" not in entry["readiness"]
    ]
    return {
        "best_pair_any": entries[0],
        "best_pair_cross_side": cross_side[0],
        "best_pair_cross_side_non_tricky": non_tricky_cross_side[0],
    }


def write_candidate_table(path: str | Path, rows: list[dict[str, object]]) -> Path:
    header = [
        "name",
        "readiness",
        "system_type",
        "mp_msun",
        "mc_msun",
        "s_mid",
        "delta_s_from_sstar",
        "gamma_obs_ms",
        "gamma_sigma_ms",
        "gamma_fractional_sigma",
        "xc",
        "new_abs95_kappa_star",
        "improvement_factor",
        "nongamma_strategy",
    ]
    lines = ["\t".join(header)]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    str(row["name"]),
                    str(row["readiness"]),
                    str(row["system_type"]),
                    f"{float(row['mp_msun']):.6f}",
                    f"{float(row['mc_msun']):.6f}",
                    f"{float(row['s_mid']):.6f}",
                    f"{float(row['delta_s_from_sstar']):+.6f}",
                    f"{1.0e3 * float(row['gamma_obs_s']):.6f}",
                    f"{1.0e3 * float(row['gamma_sigma_s']):.6f}",
                    f"{float(row['gamma_fractional_sigma']):.6e}",
                    f"{float(row['xc']):.6f}",
                    f"{float(row['new_abs95_kappa_star']):.6e}",
                    f"{float(row['improvement_factor']):.6f}",
                    str(row["nongamma_strategy"]),
                ]
            )
        )
    output = Path(path)
    output.write_text("\n".join(lines) + "\n")
    return output


def _map_x(value: float, xmin: float, xmax: float, x0: float, width: float) -> float:
    return x0 + width * (value - xmin) / (xmax - xmin)


def _map_y_log(value: float, ymin: float, ymax: float, y0: float, height: float) -> float:
    return y0 + height * (math.log10(ymax) - math.log10(value)) / (math.log10(ymax) - math.log10(ymin))


def write_summary_svg(path: str | Path, summary: dict[str, object]) -> Path:
    rows = summary["candidate_rows_ranked"]
    lever = summary["lever_arm_target_reference"]
    x_values = [float(row["delta_s_from_sstar"]) for row in rows]
    y_values = [float(row["gamma_fractional_sigma"]) for row in rows]
    xmin = min(-0.02, min(x_values) - 0.005)
    xmax = max(0.07, max(x_values) + 0.005)
    ymin = min(1.0e-4, min(y_values) * 0.8)
    ymax = max(5.0e-2, max(y_values) * 1.2)

    x0 = 90.0
    y0 = 150.0
    width = 760.0
    height = 340.0
    colors = {
        "clean": "#0f766e",
        "moderate": "#b45309",
        "conditional": "#7c3aed",
        "tricky": "#b91c1c",
    }

    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="980" height="620" viewBox="0 0 980 620">',
        '<rect width="980" height="620" fill="#f8fafc"/>',
        '<text x="48" y="56" font-family="sans-serif" font-size="24" font-weight="700" fill="#0f172a">Request 6 Source Scout</text>',
        '<text x="48" y="84" font-family="sans-serif" font-size="14" fill="#334155">Actual gamma-measured systems were projected onto the current clock-only basis around s* to see which sources add slope information rather than just more of the same amplitude direction.</text>',
        f'<text x="48" y="110" font-family="monospace" font-size="13" fill="#475569">current |kappa_*|_95 = {summary["current_basis_gaussian"]["abs95_kappa_star"]:.3e}; audit target: |Delta s| ~ {lever["target_abs_ds"]:.3f}, sigma_delta ~ {lever["target_sigma_delta"]:.1e}</text>',
        f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{width:.1f}" height="{height:.1f}" fill="white" stroke="#0f172a" stroke-width="1.2"/>',
    ]

    for x_tick in np.linspace(xmin, xmax, 7):
        x = _map_x(float(x_tick), xmin, xmax, x0, width)
        parts.append(f'<line x1="{x:.1f}" y1="{y0:.1f}" x2="{x:.1f}" y2="{y0 + height:.1f}" stroke="#e2e8f0" stroke-width="1"/>')
        parts.append(f'<text x="{x:.1f}" y="{y0 + height + 24:.1f}" text-anchor="middle" font-family="monospace" font-size="12" fill="#475569">{x_tick:+.3f}</text>')
    for y_tick in [1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2]:
        if ymin <= y_tick <= ymax:
            y = _map_y_log(y_tick, ymin, ymax, y0, height)
            parts.append(f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x0 + width:.1f}" y2="{y:.1f}" stroke="#e2e8f0" stroke-width="1"/>')
            parts.append(f'<text x="{x0 - 10:.1f}" y="{y + 4:.1f}" text-anchor="end" font-family="monospace" font-size="12" fill="#475569">{y_tick:.0e}</text>')

    x_pos = _map_x(lever["target_abs_ds"], xmin, xmax, x0, width)
    x_neg = _map_x(-lever["target_abs_ds"], xmin, xmax, x0, width)
    y_target = _map_y_log(lever["target_sigma_delta"], ymin, ymax, y0, height)
    parts.append(f'<line x1="{x_pos:.1f}" y1="{y0:.1f}" x2="{x_pos:.1f}" y2="{y0 + height:.1f}" stroke="#2563eb" stroke-width="1.5" stroke-dasharray="5 4"/>')
    parts.append(f'<line x1="{x_neg:.1f}" y1="{y0:.1f}" x2="{x_neg:.1f}" y2="{y0 + height:.1f}" stroke="#2563eb" stroke-width="1.5" stroke-dasharray="5 4"/>')
    parts.append(f'<line x1="{x0:.1f}" y1="{y_target:.1f}" x2="{x0 + width:.1f}" y2="{y_target:.1f}" stroke="#2563eb" stroke-width="1.5" stroke-dasharray="5 4"/>')
    parts.append(f'<text x="{x0 + width - 4:.1f}" y="{y_target - 8:.1f}" text-anchor="end" font-family="sans-serif" font-size="12" fill="#1d4ed8">single-source precision target</text>')

    for row in rows:
        x = _map_x(float(row["delta_s_from_sstar"]), xmin, xmax, x0, width)
        y = _map_y_log(float(row["gamma_fractional_sigma"]), ymin, ymax, y0, height)
        color = colors[str(row["readiness"])]
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6.2" fill="{color}" stroke="white" stroke-width="1.4"/>')
        parts.append(f'<text x="{x + 8:.1f}" y="{y - 8:.1f}" font-family="sans-serif" font-size="12" fill="#0f172a">{row["name"]}</text>')

    parts.append(f'<text x="{x0 + width/2:.1f}" y="{y0 + height + 54:.1f}" text-anchor="middle" font-family="sans-serif" font-size="14" fill="#0f172a">delta s = s_candidate - s*</text>')
    parts.append(f'<text x="28" y="{y0 + height/2:.1f}" transform="rotate(-90 28 {y0 + height/2:.1f})" text-anchor="middle" font-family="sans-serif" font-size="14" fill="#0f172a">fractional gamma uncertainty (proxy for sigma_delta)</text>')

    legend_y = 540.0
    legend_x = 96.0
    for idx, label in enumerate(["clean", "moderate", "conditional", "tricky"]):
        parts.append(f'<circle cx="{legend_x + 140*idx:.1f}" cy="{legend_y:.1f}" r="6" fill="{colors[label]}"/>')
        parts.append(f'<text x="{legend_x + 140*idx + 12:.1f}" y="{legend_y + 4:.1f}" font-family="sans-serif" font-size="13" fill="#0f172a">{label}</text>')

    parts.append("</svg>")
    output = Path(path)
    output.write_text("\n".join(parts))
    return output


def build_summary(output_dir: str | Path = ".") -> dict[str, object]:
    output_root = Path(output_dir)
    request6_summary = load_request6_summary(output_root / "request6_clock_sector_summary.json")
    audit_summary = load_lever_arm_audit(output_root / "request6_lever_arm_audit_summary.json")

    base_rows = audit.current_rows(request6_summary)
    current_basis = audit.fisher_eta_kappa(
        base_rows,
        float(request6_summary["effective_combinations"]["linearized_basis"]["clock_only"]["sstar"]),
    )
    sstar = float(current_basis["sstar"])
    current_kappa95 = float(current_basis["abs95_kappa_star"])

    rows = []
    for candidate in CANDIDATES:
        row = candidate_row(candidate, sstar)
        row.update(scout_metrics(row, base_rows, sstar, current_kappa95))
        rows.append(row)

    rows.sort(
        key=lambda item: (
            -float(item["composite_priority_score"]),
            float(item["new_abs95_kappa_star"]),
        )
    )

    pair_summary = best_pair_summaries(rows, base_rows, sstar, current_kappa95)
    best_single = min(rows, key=lambda item: float(item["new_abs95_kappa_star"]))
    best_positive = min(
        [row for row in rows if float(row["delta_s_from_sstar"]) > 0.0],
        key=lambda item: float(item["new_abs95_kappa_star"]),
    )
    best_negative = min(
        [row for row in rows if float(row["delta_s_from_sstar"]) < 0.0],
        key=lambda item: float(item["new_abs95_kappa_star"]),
    )

    target_pair = min(
        audit_summary["representative_pair_rows"],
        key=lambda item: abs(float(item["sigma_delta"]) - 2.0e-4),
    )

    summary = {
        "current_basis_gaussian": current_basis,
        "current_rows": base_rows,
        "lever_arm_target_reference": {
            "target_abs_ds": float(target_pair["best_abs_ds"]),
            "target_sigma_delta": float(target_pair["sigma_delta"]),
        },
        "candidate_rows_ranked": rows,
        "best_single_candidate": {
            "name": str(best_single["name"]),
            "new_abs95_kappa_star": float(best_single["new_abs95_kappa_star"]),
            "improvement_factor": float(best_single["improvement_factor"]),
        },
        "best_positive_side_candidate": {
            "name": str(best_positive["name"]),
            "delta_s_from_sstar": float(best_positive["delta_s_from_sstar"]),
            "new_abs95_kappa_star": float(best_positive["new_abs95_kappa_star"]),
        },
        "best_negative_side_candidate": {
            "name": str(best_negative["name"]),
            "delta_s_from_sstar": float(best_negative["delta_s_from_sstar"]),
            "new_abs95_kappa_star": float(best_negative["new_abs95_kappa_star"]),
        },
        "pair_summaries": pair_summary,
        "sources": [asdict(candidate) for candidate in CANDIDATES],
    }

    (output_root / "request6_source_scout_summary.json").write_text(json.dumps(summary, indent=2))
    write_candidate_table(output_root / "request6_source_scout_candidates.tsv", rows)
    write_summary_svg(output_root / "request6_source_scout_summary.svg", summary)
    return summary


def print_summary(summary: dict[str, object]) -> None:
    print("=== Request 6 source scout ===")
    current = summary["current_basis_gaussian"]
    target = summary["lever_arm_target_reference"]
    print(
        f"current |kappa_*|_95 = {current['abs95_kappa_star']:.3e} at s* = {current['sstar']:.6f}; "
        f"audit target ~ (|Delta s|={target['target_abs_ds']:.3f}, sigma_delta={target['target_sigma_delta']:.1e})"
    )
    print()
    print("=== Candidate ranking ===")
    for row in summary["candidate_rows_ranked"]:
        print(
            f"{row['name']}: readiness={row['readiness']}, "
            f"delta_s={row['delta_s_from_sstar']:+.4f}, "
            f"sigma_delta={row['gamma_fractional_sigma']:.3e}, "
            f"new |kappa_*|_95={row['new_abs95_kappa_star']:.3e}, "
            f"improvement={row['improvement_factor']:.2f}x"
        )
    print()
    print("=== Best combinations ===")
    best_single = summary["best_single_candidate"]
    print(
        f"best single: {best_single['name']} -> "
        f"|kappa_*|_95={best_single['new_abs95_kappa_star']:.3e} "
        f"({best_single['improvement_factor']:.2f}x)"
    )
    best_pair = summary["pair_summaries"]["best_pair_any"]
    print(
        f"best pair(any): {', '.join(best_pair['names'])} -> "
        f"|kappa_*|_95={best_pair['new_abs95_kappa_star']:.3e}"
    )
    cross = summary["pair_summaries"]["best_pair_cross_side_non_tricky"]
    print(
        f"best cross-side non-tricky pair: {', '.join(cross['names'])} -> "
        f"|kappa_*|_95={cross['new_abs95_kappa_star']:.3e}"
    )


def main(output_dir: str | Path = ".") -> None:
    summary = build_summary(output_dir)
    print_summary(summary)


if __name__ == "__main__":
    main()
