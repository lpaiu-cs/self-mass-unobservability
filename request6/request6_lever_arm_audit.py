from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

import request6_clock_sector as r6


@dataclass(frozen=True)
class LeverArmConfig:
    candidate_s_min: float = 0.08
    candidate_s_max: float = 0.20
    candidate_s_points: int = 161
    candidate_sigma_delta_min: float = 5.0e-5
    candidate_sigma_delta_max: float = 2.0e-3
    candidate_sigma_delta_points: int = 161
    candidate_xc: float = 0.50
    target_kappa_abs95: float = 1.0e-2
    secondary_kappa_abs95: float = 5.0e-3


def load_request6_summary(summary_path: str | Path = "request6_clock_sector_summary.json") -> dict[str, object]:
    return json.loads(Path(summary_path).read_text())


def current_rows(summary: dict[str, object]) -> list[dict[str, float | str]]:
    pulls = {item["name"]: item for item in summary["models"]["clock_only"]["pulls"]}
    rows = []
    for system in r6.SYSTEMS:
        meta = summary["system_meta"][system.name]
        pull = pulls[system.name]
        sigma_total = math.sqrt(system.gamma_sigma_s**2 + pull["gamma_model_std_s"] ** 2)
        rows.append(
            {
                "name": system.name,
                "sbar": float(meta["self_gravity_mean"]),
                "xc": float(meta["xc_mean"]),
                "gamma_gr_s": float(meta["gamma_pred_gr_s"]),
                "sigma_gamma_s": float(sigma_total),
                "sigma_delta": float(sigma_total / meta["gamma_pred_gr_s"]),
                "delta_gamma_obs": float(meta["gamma_obs_over_gr_minus_one"]),
            }
        )
    return rows


def fisher_eta_kappa(rows: list[dict[str, float | str]], sstar: float) -> dict[str, object]:
    fisher = np.zeros((2, 2), dtype=float)
    for row in rows:
        sbar = float(row["sbar"])
        xc = float(row["xc"])
        sigma_delta = float(row["sigma_delta"])
        ds = sbar - sstar
        coeff = 1.0 / (1.0 + xc)
        vec = np.array([coeff, coeff * ds], dtype=float)
        fisher += (1.0 / sigma_delta**2) * np.outer(vec, vec)

    cov = np.linalg.pinv(fisher, hermitian=True)
    std_eta = float(math.sqrt(max(cov[0, 0], 0.0)))
    std_kappa = float(math.sqrt(max(cov[1, 1], 0.0)))
    corr = 0.0
    if std_eta > 0.0 and std_kappa > 0.0:
        corr = float(cov[0, 1] / (std_eta * std_kappa))
    return {
        "sstar": float(sstar),
        "covariance": cov.tolist(),
        "std_eta_star": std_eta,
        "std_kappa_star": std_kappa,
        "abs95_eta_star": 1.96 * std_eta,
        "abs95_kappa_star": 1.96 * std_kappa,
        "correlation": corr,
    }


def slope_information_proxy(rows: list[dict[str, float | str]], sstar: float) -> float:
    return float(
        sum(
            ((float(row["sbar"]) - sstar) ** 2) / ((1.0 + float(row["xc"])) ** 2 * float(row["sigma_delta"]) ** 2)
            for row in rows
        )
    )


def scenario_grid(
    base_rows: list[dict[str, float | str]],
    sstar: float,
    config: LeverArmConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    s_grid = np.linspace(config.candidate_s_min, config.candidate_s_max, config.candidate_s_points)
    sigma_grid = np.geomspace(config.candidate_sigma_delta_min, config.candidate_sigma_delta_max, config.candidate_sigma_delta_points)
    current = fisher_eta_kappa(base_rows, sstar)
    current_kappa95 = float(current["abs95_kappa_star"])

    kappa95 = np.zeros((sigma_grid.size, s_grid.size), dtype=float)
    improvement = np.zeros_like(kappa95)
    for i, sigma_delta in enumerate(sigma_grid):
        for j, s_candidate in enumerate(s_grid):
            rows = list(base_rows) + [
                {
                    "name": "candidate",
                    "sbar": float(s_candidate),
                    "xc": config.candidate_xc,
                    "gamma_gr_s": 1.0,
                    "sigma_gamma_s": float(sigma_delta),
                    "sigma_delta": float(sigma_delta),
                    "delta_gamma_obs": 0.0,
                }
            ]
            fit = fisher_eta_kappa(rows, sstar)
            kappa95[i, j] = float(fit["abs95_kappa_star"])
            improvement[i, j] = current_kappa95 / kappa95[i, j] if kappa95[i, j] > 0.0 else 0.0
    return s_grid, sigma_grid, kappa95, improvement


def symmetric_pair_grid(
    base_rows: list[dict[str, float | str]],
    sstar: float,
    config: LeverArmConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ds_max = min(sstar - config.candidate_s_min, config.candidate_s_max - sstar)
    ds_grid = np.linspace(0.0, ds_max, config.candidate_s_points)
    sigma_grid = np.geomspace(config.candidate_sigma_delta_min, config.candidate_sigma_delta_max, config.candidate_sigma_delta_points)
    current = fisher_eta_kappa(base_rows, sstar)
    current_kappa95 = float(current["abs95_kappa_star"])

    kappa95 = np.zeros((sigma_grid.size, ds_grid.size), dtype=float)
    improvement = np.zeros_like(kappa95)
    for i, sigma_delta in enumerate(sigma_grid):
        for j, ds in enumerate(ds_grid):
            rows = list(base_rows)
            for sign in (-1.0, 1.0):
                rows.append(
                    {
                        "name": f"candidate_{sign:+.0f}",
                        "sbar": float(sstar + sign * ds),
                        "xc": config.candidate_xc,
                        "gamma_gr_s": 1.0,
                        "sigma_gamma_s": float(sigma_delta),
                        "sigma_delta": float(sigma_delta),
                        "delta_gamma_obs": 0.0,
                    }
                )
            fit = fisher_eta_kappa(rows, sstar)
            kappa95[i, j] = float(fit["abs95_kappa_star"])
            improvement[i, j] = current_kappa95 / kappa95[i, j] if kappa95[i, j] > 0.0 else 0.0
    return ds_grid, sigma_grid, kappa95, improvement


def threshold_table(
    s_grid: np.ndarray,
    sigma_grid: np.ndarray,
    kappa95: np.ndarray,
    sstar: float,
    targets: list[float],
    mode: str = "single",
) -> list[dict[str, float | None]]:
    table = []
    for sigma_idx, sigma_delta in enumerate(sigma_grid):
        row = {"sigma_delta": float(sigma_delta)}
        distances = np.abs(s_grid - sstar) if mode == "single" else s_grid
        for target in targets:
            mask = kappa95[sigma_idx] <= target
            if not np.any(mask):
                row[f"reach_{target:.0e}"] = None
            else:
                row[f"reach_{target:.0e}"] = float(np.min(distances[mask]))
        table.append(row)
    return table


def pick_representative_rows(
    s_grid: np.ndarray,
    sigma_grid: np.ndarray,
    kappa95: np.ndarray,
    improvement: np.ndarray,
    sstar: float,
    mode: str = "single",
) -> list[dict[str, float]]:
    picks = []
    for sigma_delta in [1.0e-3, 5.0e-4, 2.0e-4, 1.0e-4]:
        sigma_idx = int(np.argmin(np.abs(np.log(sigma_grid) - np.log(sigma_delta))))
        best_idx = int(np.argmin(kappa95[sigma_idx]))
        picks.append(
            {
                "sigma_delta": float(sigma_grid[sigma_idx]),
                "best_s": float(s_grid[best_idx] if mode == "single" else sstar + s_grid[best_idx]),
                "best_abs_ds": float(abs(s_grid[best_idx] - sstar) if mode == "single" else s_grid[best_idx]),
                "kappa95": float(kappa95[sigma_idx, best_idx]),
                "improvement": float(improvement[sigma_idx, best_idx]),
            }
        )
    return picks


def _svg_heatmap(
    parts: list[str],
    x0: float,
    y0: float,
    w: float,
    h: float,
    xgrid: np.ndarray,
    ygrid: np.ndarray,
    values: np.ndarray,
    title: str,
    logy: bool = False,
) -> None:
    parts.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" fill="none" stroke="#111827" stroke-width="1.2"/>')
    parts.append(f'<text x="{x0:.1f}" y="{y0 - 10:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">{title}</text>')
    norm = values / np.max(values)
    for iy in range(ygrid.size - 1):
        for ix in range(xgrid.size - 1):
            frac = float(norm[iy, ix])
            color = f"rgb({248 - int(150*frac)},{250 - int(110*frac)},{252 - int(210*frac)})"
            x = x0 + w * ix / (xgrid.size - 1)
            y = y0 + h * iy / (ygrid.size - 1)
            rw = w / (xgrid.size - 1)
            rh = h / (ygrid.size - 1)
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{rw:.2f}" height="{rh:.2f}" fill="{color}" stroke="none"/>')
    for frac, label in [(0.0, f"{xgrid[0]:.2f}"), (0.5, f"{xgrid[len(xgrid)//2]:.2f}"), (1.0, f"{xgrid[-1]:.2f}")]:
        x = x0 + frac * w
        parts.append(f'<line x1="{x:.1f}" y1="{y0+h:.1f}" x2="{x:.1f}" y2="{y0+h+6:.1f}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x-18:.1f}" y="{y0+h+20:.1f}" font-family="monospace" font-size="11" fill="#374151">{label}</text>')
    y_labels = [ygrid[0], ygrid[len(ygrid)//2], ygrid[-1]]
    if logy:
        y_labels = [ygrid[-1], ygrid[len(ygrid)//2], ygrid[0]]
    for frac, label_val in zip([0.0, 0.5, 1.0], y_labels):
        y = y0 + frac * h
        parts.append(f'<line x1="{x0-6:.1f}" y1="{y:.1f}" x2="{x0:.1f}" y2="{y:.1f}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x0-82:.1f}" y="{y+4:.1f}" font-family="monospace" font-size="11" fill="#374151">{label_val:.1e}</text>')
    parts.append(f'<text x="{x0+w/2-28:.1f}" y="{y0+h+38:.1f}" font-family="sans-serif" font-size="12" fill="#111827">candidate s</text>')
    parts.append(f'<text x="{x0-64:.1f}" y="{y0+h/2:.1f}" font-family="sans-serif" font-size="12" fill="#111827" transform="rotate(-90 {x0-64:.1f},{y0+h/2:.1f})">sigma_delta</text>')


def _svg_rows(
    parts: list[str],
    x0: float,
    y0: float,
    rows: list[dict[str, float | str]],
    summary: dict[str, object],
) -> None:
    parts.append(f'<text x="{x0:.1f}" y="{y0:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">Current systems</text>')
    yy = y0 + 28
    for row in rows:
        parts.append(
            f'<text x="{x0:.1f}" y="{yy:.1f}" font-family="monospace" font-size="12" fill="#111827">'
            f'{row["name"]}: s={float(row["sbar"]):.4f}, Xc={float(row["xc"]):.3f}, sigma_delta={float(row["sigma_delta"]):.2e}'
            '</text>'
        )
        yy += 18
    diag = summary["effective_combinations"]["sbar_diagnostics"]
    parts.append(f'<text x="{x0:.1f}" y="{yy+10:.1f}" font-family="monospace" font-size="12" fill="#111827">current slope proxy = {diag["slope_information_proxy"]:.3e}</text>')


def write_summary_svg(
    path: str | Path,
    rows: list[dict[str, float | str]],
    summary: dict[str, object],
    s_grid: np.ndarray,
    sigma_grid: np.ndarray,
    kappa95: np.ndarray,
    improvement: np.ndarray,
    representative_rows: list[dict[str, float]],
    ds_grid_pair: np.ndarray,
    kappa95_pair: np.ndarray,
    improvement_pair: np.ndarray,
    representative_pair_rows: list[dict[str, float]],
) -> Path:
    output_path = Path(path)
    width, height = 1400, 980
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 6 Follow-up: Lever-Arm Audit</text>',
        '<text x="48" y="72" font-family="sans-serif" font-size="14" fill="#4b5563">Linearized Gaussian audit for the local basis (eta_*, kappa_*) around the current s*.</text>',
    ]
    _svg_heatmap(parts, 90, 140, 500, 300, s_grid, sigma_grid[::-1], improvement[::-1], "Single-candidate improvement in |kappa_*|_95", logy=True)
    _svg_heatmap(parts, 770, 140, 500, 300, ds_grid_pair, sigma_grid[::-1], improvement_pair[::-1], "Symmetric-pair improvement in |kappa_*|_95", logy=True)
    _svg_rows(parts, 110, 590, rows, summary)
    yy = 710
    parts.append(f'<text x="110" y="{yy:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">Single-candidate precision tiers</text>')
    yy += 24
    for item in representative_rows:
        parts.append(
            f'<text x="110" y="{yy:.1f}" font-family="monospace" font-size="12" fill="#111827">'
            f'sigma_delta={item["sigma_delta"]:.1e}: best s={item["best_s"]:.3f}, '
            f'|ds|={item["best_abs_ds"]:.3f}, |kappa_*|_95={item["kappa95"]:.3e}, '
            f'improvement={item["improvement"]:.2f}x'
            '</text>'
        )
        yy += 18
    yy += 20
    parts.append(f'<text x="110" y="{yy:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">Symmetric-pair precision tiers</text>')
    yy += 24
    for item in representative_pair_rows:
        parts.append(
            f'<text x="110" y="{yy:.1f}" font-family="monospace" font-size="12" fill="#111827">'
            f'sigma_delta={item["sigma_delta"]:.1e}: best |ds|={item["best_abs_ds"]:.3f}, '
            f'|kappa_*|_95={item["kappa95"]:.3e}, improvement={item["improvement"]:.2f}x'
            '</text>'
        )
        yy += 18
    parts.append("</svg>")
    output_path.write_text("\n".join(parts))
    return output_path


def run_lever_arm_audit(output_dir: str | Path = ".") -> dict[str, object]:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    config = LeverArmConfig()
    request6_summary = load_request6_summary(output_root / "request6_clock_sector_summary.json")
    rows = current_rows(request6_summary)
    sstar = float(request6_summary["effective_combinations"]["linearized_basis"]["clock_only"]["sstar"])
    current_fit = fisher_eta_kappa(rows, sstar)
    current_proxy = slope_information_proxy(rows, sstar)

    s_grid, sigma_grid, kappa95, improvement = scenario_grid(rows, sstar, config)
    ds_grid_pair, sigma_grid_pair, kappa95_pair, improvement_pair = symmetric_pair_grid(rows, sstar, config)
    thresholds = threshold_table(
        s_grid,
        sigma_grid,
        kappa95,
        sstar,
        [config.target_kappa_abs95, config.secondary_kappa_abs95],
        mode="single",
    )
    thresholds_pair = threshold_table(
        ds_grid_pair,
        sigma_grid_pair,
        kappa95_pair,
        sstar,
        [config.target_kappa_abs95, config.secondary_kappa_abs95],
        mode="pair",
    )
    representative = pick_representative_rows(s_grid, sigma_grid, kappa95, improvement, sstar)
    representative_pair = pick_representative_rows(ds_grid_pair, sigma_grid_pair, kappa95_pair, improvement_pair, sstar, mode="pair")

    summary = {
        "config": asdict(config),
        "current_rows": rows,
        "current_basis_gaussian": current_fit,
        "current_slope_information_proxy": current_proxy,
        "sstar": sstar,
        "candidate_s_grid": s_grid.tolist(),
        "candidate_sigma_delta_grid": sigma_grid.tolist(),
        "kappa95_grid": kappa95.tolist(),
        "improvement_grid": improvement.tolist(),
        "symmetric_pair_abs_ds_grid": ds_grid_pair.tolist(),
        "symmetric_pair_sigma_delta_grid": sigma_grid_pair.tolist(),
        "symmetric_pair_kappa95_grid": kappa95_pair.tolist(),
        "symmetric_pair_improvement_grid": improvement_pair.tolist(),
        "threshold_table": thresholds,
        "threshold_table_symmetric_pair": thresholds_pair,
        "representative_rows": representative,
        "representative_pair_rows": representative_pair,
    }

    summary_path = output_root / "request6_lever_arm_audit_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    svg_path = write_summary_svg(
        output_root / "request6_lever_arm_audit_summary.svg",
        rows,
        request6_summary,
        s_grid,
        sigma_grid,
        kappa95,
        improvement,
        representative,
        ds_grid_pair,
        kappa95_pair,
        improvement_pair,
        representative_pair,
    )
    return {"summary": summary, "summary_path": summary_path, "svg_path": svg_path}


def print_summary(summary: dict[str, object]) -> None:
    current = summary["current_basis_gaussian"]
    print("=== Current linearized basis ===")
    print(
        f"s* = {summary['sstar']:.6f}, "
        f"|eta_*|_95 ~= {current['abs95_eta_star']:.3e}, "
        f"|kappa_*|_95 ~= {current['abs95_kappa_star']:.3e}"
    )
    print(f"slope proxy = {summary['current_slope_information_proxy']:.3e}")
    print()
    print("=== Representative single-candidate scenarios ===")
    for item in summary["representative_rows"]:
        print(
            f"sigma_delta = {item['sigma_delta']:.1e}: "
            f"best s = {item['best_s']:.3f}, "
            f"|ds| = {item['best_abs_ds']:.3f}, "
            f"|kappa_*|_95 = {item['kappa95']:.3e}, "
            f"improvement = {item['improvement']:.2f}x"
        )
    print()
    print("=== Representative symmetric-pair scenarios ===")
    for item in summary["representative_pair_rows"]:
        print(
            f"sigma_delta = {item['sigma_delta']:.1e}: "
            f"best |ds| = {item['best_abs_ds']:.3f}, "
            f"|kappa_*|_95 = {item['kappa95']:.3e}, "
            f"improvement = {item['improvement']:.2f}x"
        )


def main(output_dir: str | Path = ".") -> None:
    result = run_lever_arm_audit(output_dir)
    print_summary(result["summary"])
    print()
    print(f"wrote {result['summary_path']}")
    print(f"wrote {result['svg_path']}")


if __name__ == "__main__":
    main()
