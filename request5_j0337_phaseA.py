from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class PhaseABounds:
    name: str
    delta_95: float


@dataclass(frozen=True)
class EOSPrior:
    name: str
    s_ns_min: float
    s_ns_max: float


@dataclass(frozen=True)
class PhaseAConfig:
    sigma1_min: float = -5.0e-5
    sigma1_max: float = 5.0e-5
    sigma1_points: int = 201
    sigma2_min: float = -4.0e-4
    sigma2_max: float = 4.0e-4
    sigma2_points: int = 241
    s_ns_points: int = 201
    s_wd: float = 1.0e-4


def gaussian_sigma_from_95(delta_95: float) -> float:
    return delta_95 / 1.96


def sigma_grids(config: PhaseAConfig) -> tuple[np.ndarray, np.ndarray]:
    sigma1 = np.linspace(config.sigma1_min, config.sigma1_max, config.sigma1_points)
    sigma2 = np.linspace(config.sigma2_min, config.sigma2_max, config.sigma2_points)
    return sigma1, sigma2


def delta_smu(sigma1: np.ndarray, sigma2: np.ndarray, s_ns: np.ndarray, s_wd: float) -> np.ndarray:
    return sigma1[..., None] * (s_ns - s_wd) + sigma2[..., None] * (s_ns**2 - s_wd**2)


def posterior_grid(
    config: PhaseAConfig,
    bound: PhaseABounds,
    prior: EOSPrior,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sigma1, sigma2 = sigma_grids(config)
    s_ns = np.linspace(prior.s_ns_min, prior.s_ns_max, config.s_ns_points)
    sigma1_3d = sigma1[:, None]
    sigma2_3d = sigma2[None, :]
    sigma_delta = gaussian_sigma_from_95(bound.delta_95)

    delta = delta_smu(sigma1_3d, sigma2_3d, s_ns, config.s_wd)
    log_like = -0.5 * (delta / sigma_delta) ** 2
    posterior = np.mean(np.exp(log_like), axis=2)
    posterior /= np.sum(posterior)
    return sigma1, sigma2, posterior


def marginal_stats(grid: np.ndarray, marginal: np.ndarray) -> dict[str, float]:
    marginal = marginal / np.sum(marginal)
    mean = float(np.sum(grid * marginal))
    std = float(np.sqrt(np.sum((grid - mean) ** 2 * marginal)))
    mode = float(grid[int(np.argmax(marginal))])

    order = np.argsort(np.abs(grid))
    cumulative = np.cumsum(marginal[order])
    idx_95 = int(np.searchsorted(cumulative, 0.95, side="left"))
    bound_95 = float(np.abs(grid[order[min(idx_95, grid.size - 1)]]))
    return {
        "mean": mean,
        "std": std,
        "mode": mode,
        "abs_95": bound_95,
    }


def summarize_posterior(sigma1: np.ndarray, sigma2: np.ndarray, posterior: np.ndarray) -> dict[str, dict[str, float]]:
    p1 = np.sum(posterior, axis=1)
    p2 = np.sum(posterior, axis=0)
    zero_s2 = int(np.argmin(np.abs(sigma2)))
    zero_s1 = int(np.argmin(np.abs(sigma1)))
    return {
        "sigma1": marginal_stats(sigma1, p1),
        "sigma2": marginal_stats(sigma2, p2),
        "sigma1_given_sigma2_eq_0": marginal_stats(sigma1, posterior[:, zero_s2]),
        "sigma2_given_sigma1_eq_0": marginal_stats(sigma2, posterior[zero_s1, :]),
    }


def write_posterior_table(path: str | Path, sigma1: np.ndarray, sigma2: np.ndarray, posterior: np.ndarray) -> Path:
    output_path = Path(path)
    rows = []
    for i, s1 in enumerate(sigma1):
        for j, s2 in enumerate(sigma2):
            rows.append([s1, s2, posterior[i, j]])
    np.savetxt(
        output_path,
        np.array(rows),
        header="sigma1\tsigma2\tposterior",
        comments="",
        delimiter="\t",
    )
    return output_path


def _svg_heatmap(
    parts: list[str],
    x0: float,
    y0: float,
    w: float,
    h: float,
    sigma1: np.ndarray,
    sigma2: np.ndarray,
    posterior: np.ndarray,
    title: str,
) -> None:
    parts.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" fill="none" stroke="#111827" stroke-width="1.2"/>')
    parts.append(f'<text x="{x0:.1f}" y="{y0 - 10:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">{title}</text>')
    pnorm = posterior / np.max(posterior)
    for i in range(sigma1.size - 1):
        for j in range(sigma2.size - 1):
            frac = float(pnorm[i, j])
            color = f"rgb({255 - int(120*frac)},{247 - int(110*frac)},{255 - int(180*frac)})"
            x = x0 + w * i / (sigma1.size - 1)
            y = y0 + h * (1.0 - (j + 1) / (sigma2.size - 1))
            rw = w / (sigma1.size - 1)
            rh = h / (sigma2.size - 1)
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{rw:.2f}" height="{rh:.2f}" fill="{color}" stroke="none"/>')
    for frac, label in [(0.0, f"{sigma1[0]:.1e}"), (0.5, "0"), (1.0, f"{sigma1[-1]:.1e}")]:
        x = x0 + frac * w
        parts.append(f'<line x1="{x:.1f}" y1="{y0+h:.1f}" x2="{x:.1f}" y2="{y0+h+6:.1f}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x-18:.1f}" y="{y0+h+20:.1f}" font-family="monospace" font-size="11" fill="#374151">{label}</text>')
    for frac, label in [(0.0, f"{sigma2[0]:.1e}"), (0.5, "0"), (1.0, f"{sigma2[-1]:.1e}")]:
        y = y0 + h * (1.0 - frac)
        parts.append(f'<line x1="{x0-6:.1f}" y1="{y:.1f}" x2="{x0:.1f}" y2="{y:.1f}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x0-72:.1f}" y="{y+4:.1f}" font-family="monospace" font-size="11" fill="#374151">{label}</text>')
    parts.append(f'<text x="{x0 + w/2 - 20:.1f}" y="{y0+h+40:.1f}" font-family="sans-serif" font-size="12" fill="#111827">sigma1</text>')
    parts.append(f'<text x="{x0-60:.1f}" y="{y0+h/2:.1f}" font-family="sans-serif" font-size="12" fill="#111827" transform="rotate(-90 {x0-60:.1f},{y0+h/2:.1f})">sigma2</text>')


def write_summary_svg(path: str | Path, summary: dict[str, object]) -> Path:
    output_path = Path(path)
    optimistic = summary["bounds"]["optimistic"]["wide"]
    conservative = summary["bounds"]["conservative"]["wide"]

    width, height = 1280, 900
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 5 Phase A: J0337 Published-Bound Translation</text>',
        '<text x="48" y="72" font-family="sans-serif" font-size="14" fill="#4b5563">Posterior mapping from published Delta bounds to (sigma1, sigma2), marginalized over s_NS priors.</text>',
    ]
    _svg_heatmap(
        parts,
        110,
        130,
        400,
        320,
        np.array(summary["sigma1_grid"]),
        np.array(summary["sigma2_grid"]),
        np.array(optimistic["posterior"]),
        "Optimistic bound: |Delta| < 1.5e-6 (95%)",
    )
    _svg_heatmap(
        parts,
        700,
        130,
        400,
        320,
        np.array(summary["sigma1_grid"]),
        np.array(summary["sigma2_grid"]),
        np.array(conservative["posterior"]),
        "Conservative bound: |Delta| < 2.3e-6 (95%)",
    )
    parts.append('<text x="90" y="545" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">EOS sensitivity summary (95% bounds on |sigma|)</text>')
    y = 585
    parts.append('<text x="90" y="585" font-family="monospace" font-size="13" fill="#374151">prior                optimistic sigma1   optimistic sigma2   conservative sigma1   conservative sigma2</text>')
    for key in ["low", "wide", "high"]:
        opt = summary["bounds"]["optimistic"][key]["stats"]
        cons = summary["bounds"]["conservative"][key]["stats"]
        parts.append(
            f'<text x="90" y="{y+32:.1f}" font-family="monospace" font-size="13" fill="#111827">'
            f'{key:<8} [{summary["priors"][key]["s_ns_min"]:.2f},{summary["priors"][key]["s_ns_max"]:.2f}]   '
            f'{opt["sigma1"]["abs_95"]:.2e}         {opt["sigma2"]["abs_95"]:.2e}         '
            f'{cons["sigma1"]["abs_95"]:.2e}         {cons["sigma2"]["abs_95"]:.2e}</text>'
        )
        y += 34
    parts.append('<text x="90" y="760" font-family="sans-serif" font-size="13" fill="#4b5563">Interpretation: sigma1 is constrained at the 10^-5 level, while sigma2 remains broader but still lands in the 10^-4 regime once strong-field s_NS ~ 0.1-0.2 is used.</text>')
    parts.append("</svg>")
    output_path.write_text("\n".join(parts))
    return output_path


def run_phase_a(output_dir: str | Path = ".") -> dict[str, object]:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    config = PhaseAConfig()
    sigma1, sigma2 = sigma_grids(config)

    bounds = {
        "optimistic": PhaseABounds("optimistic", 1.5e-6),
        "conservative": PhaseABounds("conservative", 2.3e-6),
    }
    priors = {
        "low": EOSPrior("low", 0.10, 0.15),
        "wide": EOSPrior("wide", 0.10, 0.20),
        "high": EOSPrior("high", 0.15, 0.20),
    }

    summary: dict[str, object] = {
        "config": asdict(config),
        "sigma1_grid": sigma1.tolist(),
        "sigma2_grid": sigma2.tolist(),
        "priors": {name: asdict(prior) for name, prior in priors.items()},
        "bounds": {},
    }

    for bound_name, bound in bounds.items():
        bound_summary: dict[str, object] = {}
        for prior_name, prior in priors.items():
            g1, g2, posterior = posterior_grid(config, bound, prior)
            stats = summarize_posterior(g1, g2, posterior)
            payload = {
                "delta_95": bound.delta_95,
                "gaussian_sigma": gaussian_sigma_from_95(bound.delta_95),
                "stats": stats,
            }
            if prior_name == "wide":
                payload["posterior"] = posterior.tolist()
                write_posterior_table(output_root / f"request5_j0337_phaseA_posterior_{bound_name}.tsv", g1, g2, posterior)
            bound_summary[prior_name] = payload
        summary["bounds"][bound_name] = bound_summary

    summary_path = output_root / "request5_j0337_phaseA_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    svg_path = write_summary_svg(output_root / "request5_j0337_phaseA_summary.svg", summary)
    return {"summary": summary, "summary_path": summary_path, "svg_path": svg_path}


def print_summary(summary: dict[str, object]) -> None:
    print("=== J0337 Phase A translation ===")
    for bound_name in ["optimistic", "conservative"]:
        wide = summary["bounds"][bound_name]["wide"]["stats"]
        print(
            f"{bound_name}: "
            f"sigma1 = {wide['sigma1']['mean']:.3e} ± {wide['sigma1']['std']:.3e}, "
            f"|sigma1|_95 = {wide['sigma1']['abs_95']:.3e}; "
            f"sigma2 = {wide['sigma2']['mean']:.3e} ± {wide['sigma2']['std']:.3e}, "
            f"|sigma2|_95 = {wide['sigma2']['abs_95']:.3e}; "
            f"conditional(|sigma1|,|sigma2|)_95 = "
            f"({wide['sigma1_given_sigma2_eq_0']['abs_95']:.3e}, "
            f"{wide['sigma2_given_sigma1_eq_0']['abs_95']:.3e})"
        )
    print()
    print("=== EOS sensitivity (95% bounds) ===")
    for prior_name in ["low", "wide", "high"]:
        opt = summary["bounds"]["optimistic"][prior_name]["stats"]
        cons = summary["bounds"]["conservative"][prior_name]["stats"]
        prior = summary["priors"][prior_name]
        print(
            f"{prior_name} [{prior['s_ns_min']:.2f}, {prior['s_ns_max']:.2f}]: "
            f"opt(|s1|,|s2|)=({opt['sigma1']['abs_95']:.3e}, {opt['sigma2']['abs_95']:.3e}), "
            f"cons(|s1|,|s2|)=({cons['sigma1']['abs_95']:.3e}, {cons['sigma2']['abs_95']:.3e})"
        )


def main(output_dir: str | Path = ".") -> None:
    result = run_phase_a(output_dir)
    print_summary(result["summary"])
    print()
    print(f"wrote {result['summary_path']}")
    print(f"wrote {result['svg_path']}")


if __name__ == "__main__":
    main()
