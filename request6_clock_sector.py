from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


T_SUN = 4.925490947e-6  # s
SECONDS_PER_DAY = 86_400.0
SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY


@dataclass(frozen=True)
class ClockConfig:
    observable_model: str = "potential_only"
    zeta1_min: float = -4.0e-2
    zeta1_max: float = 4.0e-2
    zeta1_points: int = 241
    zeta2_min: float = -2.5e-1
    zeta2_max: float = 2.5e-1
    zeta2_points: int = 241
    mass_samples: int = 4096
    batch_size: int = 128
    rng_seed: int = 602214


@dataclass(frozen=True)
class PulsarClockSystem:
    name: str
    pb_days: float
    pb_sigma_days: float
    eccentricity: float
    eccentricity_sigma: float
    projected_semimajor_axis_s: float
    projected_semimajor_axis_sigma_s: float
    gamma_obs_s: float
    gamma_sigma_s: float
    inference_kind: str
    source_label: str
    source_url: str
    inference_note: str
    shapiro_range_us: float | None = None
    shapiro_range_sigma_us: float | None = None
    shapiro_shape: float | None = None
    shapiro_shape_sigma: float | None = None
    shapiro_log_shape: float | None = None
    shapiro_log_shape_sigma: float | None = None
    periastron_advance_deg_yr: float | None = None
    periastron_advance_sigma_deg_yr: float | None = None
    mass_ratio: float | None = None
    mass_ratio_sigma: float | None = None


SYSTEMS = [
    PulsarClockSystem(
        name="PSR B1534+12",
        pb_days=0.420737298879,
        pb_sigma_days=2.0e-12,
        eccentricity=0.27367752,
        eccentricity_sigma=7.0e-8,
        projected_semimajor_axis_s=3.7294636,
        projected_semimajor_axis_sigma_s=6.0e-7,
        gamma_obs_s=2.0708e-3,
        gamma_sigma_s=5.0e-7,
        inference_kind="omdot_plus_s",
        source_label="Fonseca et al. (2014)",
        source_url="https://arxiv.org/abs/1402.4836",
        inference_note="Leave-one-out masses inferred from non-gamma PK parameters dotomega, s, and r.",
        shapiro_range_us=6.6,
        shapiro_range_sigma_us=0.2,
        shapiro_shape=0.9772,
        shapiro_shape_sigma=0.0016,
        periastron_advance_deg_yr=1.7557950,
        periastron_advance_sigma_deg_yr=1.9e-6,
    ),
    PulsarClockSystem(
        name="PSR J0737-3039A/B",
        pb_days=0.1022515592973,
        pb_sigma_days=1.0e-12,
        eccentricity=0.087777023,
        eccentricity_sigma=6.1e-8,
        projected_semimajor_axis_s=1.415028603,
        projected_semimajor_axis_sigma_s=9.2e-8,
        gamma_obs_s=3.84045e-4,
        gamma_sigma_s=9.4e-8,
        inference_kind="mass_ratio_plus_s",
        source_label="Kramer et al. (2021), with R quoted from Kramer et al. (2006)",
        source_url="https://arxiv.org/abs/2112.06795",
        inference_note="Leave-one-out masses inferred from theory-independent mass ratio R and non-gamma PK information s and r; gamma excluded from the mass inference.",
        shapiro_range_us=6.162,
        shapiro_range_sigma_us=0.021,
        shapiro_log_shape=9.65,
        shapiro_log_shape_sigma=0.15,
        mass_ratio=1.0714,
        mass_ratio_sigma=0.0011,
    ),
]


def gamma_gr(pb_days: np.ndarray, e: np.ndarray, mp_msun: np.ndarray, mc_msun: np.ndarray) -> np.ndarray:
    pb = pb_days * SECONDS_PER_DAY
    mtot = mp_msun + mc_msun
    return (
        T_SUN ** (2.0 / 3.0)
        * (pb / (2.0 * math.pi)) ** (1.0 / 3.0)
        * e
        * mc_msun
        * (mp_msun + 2.0 * mc_msun)
        / mtot ** (4.0 / 3.0)
    )


def gamma_clock(
    pb_days: np.ndarray,
    e: np.ndarray,
    mp_msun: np.ndarray,
    mc_msun: np.ndarray,
    self_gravity_fraction: np.ndarray,
    zeta1: np.ndarray,
    zeta2: np.ndarray,
    observable_model: str,
) -> np.ndarray:
    gamma0 = gamma_gr(pb_days, e, mp_msun, mc_msun)
    eta = zeta1[:, None, None] * self_gravity_fraction[None, None, :] + zeta2[None, :, None] * self_gravity_fraction[None, None, :] ** 2
    if observable_model == "potential_only":
        xc = mc_msun / (mp_msun + mc_msun)
        return gamma0[None, None, :] * (1.0 + eta / (1.0 + xc[None, None, :]))
    if observable_model == "full_rescale":
        return gamma0[None, None, :] * (1.0 + eta)
    raise ValueError(f"unknown observable_model={observable_model!r}")


def total_mass_from_omdot(pb_days: np.ndarray, e: np.ndarray, omdot_deg_yr: np.ndarray) -> np.ndarray:
    pb = pb_days * SECONDS_PER_DAY
    omdot = omdot_deg_yr * math.pi / 180.0 / SECONDS_PER_YEAR
    factor = omdot * (1.0 - e**2) / 3.0
    return (factor ** 1.5) * (pb / (2.0 * math.pi)) ** 2.5 / T_SUN


def mass_function(projected_semimajor_axis_s: np.ndarray, pb_days: np.ndarray) -> np.ndarray:
    pb = pb_days * SECONDS_PER_DAY
    return (4.0 * math.pi**2 * projected_semimajor_axis_s**3) / (T_SUN * pb**2)


def mass_conditioned_s_prior(mp_msun: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    center = 0.10 + 0.25 * (mp_msun - 1.20)
    lower = np.maximum(0.10, center - 0.02)
    upper = np.minimum(0.20, center + 0.02)
    return lower, upper


def _draw_valid_samples(system: PulsarClockSystem, config: ClockConfig, rng: np.random.Generator) -> dict[str, np.ndarray]:
    remaining = config.mass_samples
    collected: dict[str, list[np.ndarray]] = {"pb_days": [], "eccentricity": [], "mp_msun": [], "mc_msun": [], "self_gravity_fraction": []}

    while remaining > 0:
        draw_count = max(512, 2 * remaining)
        pb = rng.normal(system.pb_days, system.pb_sigma_days, draw_count)
        eccentricity = rng.normal(system.eccentricity, system.eccentricity_sigma, draw_count)
        x = rng.normal(system.projected_semimajor_axis_s, system.projected_semimajor_axis_sigma_s, draw_count)

        if system.inference_kind == "omdot_plus_s":
            assert system.periastron_advance_deg_yr is not None
            assert system.periastron_advance_sigma_deg_yr is not None
            assert system.shapiro_shape is not None
            assert system.shapiro_shape_sigma is not None
            omdot = rng.normal(system.periastron_advance_deg_yr, system.periastron_advance_sigma_deg_yr, draw_count)
            shapiro_shape = rng.normal(system.shapiro_shape, system.shapiro_shape_sigma, draw_count)
            total_mass = total_mass_from_omdot(pb, eccentricity, omdot)
            mf = mass_function(x, pb)
            companion_mass = np.cbrt(mf * total_mass**2 / np.clip(shapiro_shape, 1.0e-12, None) ** 3)
            pulsar_mass = total_mass - companion_mass
            valid = (
                np.isfinite(total_mass)
                & np.isfinite(companion_mass)
                & (pb > 0.0)
                & (eccentricity > 0.0)
                & (eccentricity < 1.0)
                & (shapiro_shape > 0.0)
                & (shapiro_shape < 1.0)
                & (pulsar_mass > 0.0)
                & (companion_mass > 0.0)
            )
        elif system.inference_kind == "mass_ratio_plus_s":
            assert system.mass_ratio is not None
            assert system.mass_ratio_sigma is not None
            assert system.shapiro_log_shape is not None
            assert system.shapiro_log_shape_sigma is not None
            mass_ratio = rng.normal(system.mass_ratio, system.mass_ratio_sigma, draw_count)
            shapiro_log_shape = rng.normal(system.shapiro_log_shape, system.shapiro_log_shape_sigma, draw_count)
            shapiro_shape = 1.0 - np.exp(-shapiro_log_shape)
            mf = mass_function(x, pb)
            companion_mass = mf * (1.0 + mass_ratio) ** 2 / np.clip(shapiro_shape, 1.0e-12, None) ** 3
            pulsar_mass = mass_ratio * companion_mass
            valid = (
                np.isfinite(companion_mass)
                & (pb > 0.0)
                & (eccentricity > 0.0)
                & (eccentricity < 1.0)
                & (mass_ratio > 0.0)
                & (shapiro_shape > 0.0)
                & (shapiro_shape < 1.0)
                & (pulsar_mass > 0.0)
                & (companion_mass > 0.0)
            )
        else:
            raise ValueError(f"unsupported inference_kind={system.inference_kind!r}")

        valid_count = int(np.sum(valid))
        if valid_count == 0:
            continue

        take = min(remaining, valid_count)
        valid_indices = np.flatnonzero(valid)[:take]
        pulsar_mass = pulsar_mass[valid_indices]
        companion_mass = companion_mass[valid_indices]
        pb = pb[valid_indices]
        eccentricity = eccentricity[valid_indices]
        sg_lo, sg_hi = mass_conditioned_s_prior(pulsar_mass)
        sg = rng.uniform(sg_lo, sg_hi)

        collected["pb_days"].append(pb)
        collected["eccentricity"].append(eccentricity)
        collected["mp_msun"].append(pulsar_mass)
        collected["mc_msun"].append(companion_mass)
        collected["self_gravity_fraction"].append(sg)
        remaining -= take

    return {key: np.concatenate(value) for key, value in collected.items()}


def zeta_grids(config: ClockConfig) -> tuple[np.ndarray, np.ndarray]:
    z1 = np.linspace(config.zeta1_min, config.zeta1_max, config.zeta1_points)
    z2 = np.linspace(config.zeta2_min, config.zeta2_max, config.zeta2_points)
    return z1, z2


def normalize(grid: np.ndarray) -> np.ndarray:
    total = float(np.sum(grid))
    return grid / total


def marginal_stats(grid: np.ndarray, marginal: np.ndarray) -> dict[str, float]:
    marginal = marginal / np.sum(marginal)
    mean = float(np.sum(grid * marginal))
    std = float(np.sqrt(np.sum((grid - mean) ** 2 * marginal)))
    mode = float(grid[int(np.argmax(marginal))])
    order = np.argsort(np.abs(grid))
    cumulative = np.cumsum(marginal[order])
    idx_95 = int(np.searchsorted(cumulative, 0.95, side="left"))
    bound_95 = float(np.abs(grid[order[min(idx_95, grid.size - 1)]]))
    return {"mean": mean, "std": std, "mode": mode, "abs_95": bound_95}


def summarize_posterior(z1: np.ndarray, z2: np.ndarray, posterior: np.ndarray) -> dict[str, dict[str, float]]:
    p1 = np.sum(posterior, axis=1)
    p2 = np.sum(posterior, axis=0)
    return {"zeta1": marginal_stats(z1, p1), "zeta2": marginal_stats(z2, p2)}


def _weighted_abs95(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(np.abs(values))
    cumulative = np.cumsum(weights[order])
    idx = int(np.searchsorted(cumulative, 0.95, side="left"))
    return float(np.abs(values[order[min(idx, values.size - 1)]]))


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    idx = int(np.searchsorted(cumulative, quantile, side="left"))
    return float(sorted_values[min(idx, sorted_values.size - 1)])


def eta_stats_at_sbar(z1: np.ndarray, z2: np.ndarray, posterior: np.ndarray, sbar: float) -> dict[str, float]:
    eta_grid = z1[:, None] * sbar + z2[None, :] * sbar**2
    weights = posterior.reshape(-1)
    values = eta_grid.reshape(-1)
    mean = float(np.sum(weights * values))
    std = float(np.sqrt(np.sum(weights * (values - mean) ** 2)))
    mode = float(values[int(np.argmax(weights))])
    return {
        "sbar": float(sbar),
        "mean": mean,
        "std": std,
        "mode": mode,
        "abs_95": _weighted_abs95(values, weights),
        "q16": _weighted_quantile(values, weights, 0.16),
        "q84": _weighted_quantile(values, weights, 0.84),
        "q025": _weighted_quantile(values, weights, 0.025),
        "q975": _weighted_quantile(values, weights, 0.975),
    }


def basis_stats_at_sstar(z1: np.ndarray, z2: np.ndarray, posterior: np.ndarray, sstar: float) -> dict[str, dict[str, float] | float]:
    weights = posterior.reshape(-1)
    z1_flat = np.repeat(z1, z2.size)
    z2_flat = np.tile(z2, z1.size)
    eta_star = z1_flat * sstar + z2_flat * sstar**2
    kappa_star = z1_flat + 2.0 * sstar * z2_flat

    def pack(values: np.ndarray) -> dict[str, float]:
        mean = float(np.sum(weights * values))
        std = float(np.sqrt(np.sum(weights * (values - mean) ** 2)))
        mode = float(values[int(np.argmax(weights))])
        return {
            "mean": mean,
            "std": std,
            "mode": mode,
            "abs_95": _weighted_abs95(values, weights),
            "q16": _weighted_quantile(values, weights, 0.16),
            "q84": _weighted_quantile(values, weights, 0.84),
            "q025": _weighted_quantile(values, weights, 0.025),
            "q975": _weighted_quantile(values, weights, 0.975),
        }

    return {"sstar": float(sstar), "eta_star": pack(eta_star), "kappa_star": pack(kappa_star)}


def delta_gamma_stats_from_eta(eta_stats: dict[str, float], xc: float) -> dict[str, float]:
    scale = 1.0 / (1.0 + xc)
    return {
        "xc": float(xc),
        "mean": eta_stats["mean"] * scale,
        "std": eta_stats["std"] * scale,
        "mode": eta_stats["mode"] * scale,
        "abs_95": eta_stats["abs_95"] * scale,
        "q16": eta_stats["q16"] * scale,
        "q84": eta_stats["q84"] * scale,
        "q025": eta_stats["q025"] * scale,
        "q975": eta_stats["q975"] * scale,
    }


def system_likelihood_grid(
    system: PulsarClockSystem,
    z1: np.ndarray,
    z2: np.ndarray,
    config: ClockConfig,
    rng: np.random.Generator,
) -> tuple[np.ndarray, dict[str, float], dict[str, np.ndarray]]:
    samples = _draw_valid_samples(system, config, rng)
    pb = samples["pb_days"]
    eccentricity = samples["eccentricity"]
    mp = samples["mp_msun"]
    mc = samples["mc_msun"]
    sg = samples["self_gravity_fraction"]
    gamma0 = gamma_gr(pb, eccentricity, mp, mc)
    weights = np.ones_like(gamma0)
    if system.shapiro_range_us is not None and system.shapiro_range_sigma_us is not None:
        shapiro_range_model = 1.0e6 * T_SUN * mc
        weights *= np.exp(-0.5 * ((system.shapiro_range_us - shapiro_range_model) / system.shapiro_range_sigma_us) ** 2)
    weights /= np.sum(weights)

    like_sum = np.zeros((z1.size, z2.size), dtype=float)
    for start in range(0, gamma0.size, config.batch_size):
        stop = min(start + config.batch_size, gamma0.size)
        model = gamma_clock(
            pb[start:stop],
            eccentricity[start:stop],
            mp[start:stop],
            mc[start:stop],
            sg[start:stop],
            z1,
            z2,
            config.observable_model,
        )
        resid = system.gamma_obs_s - model
        like_sum += (np.exp(-0.5 * (resid / system.gamma_sigma_s) ** 2) * weights[None, None, start:stop]).sum(axis=2)
    likelihood = like_sum

    meta = {
        "mp_mean_msun": float(np.sum(weights * mp)),
        "mp_std_msun": float(np.sqrt(np.sum(weights * (mp - np.sum(weights * mp)) ** 2))),
        "mc_mean_msun": float(np.sum(weights * mc)),
        "mc_std_msun": float(np.sqrt(np.sum(weights * (mc - np.sum(weights * mc)) ** 2))),
        "xc_mean": float(np.sum(weights * (mc / (mp + mc)))),
        "self_gravity_mean": float(np.sum(weights * sg)),
        "self_gravity_std": float(np.sqrt(np.sum(weights * (sg - np.sum(weights * sg)) ** 2))),
        "gamma_pred_gr_s": float(np.sum(weights * gamma0)),
        "gamma_pred_gr_std_s": float(np.sqrt(np.sum(weights * (gamma0 - np.sum(weights * gamma0)) ** 2))),
        "gamma_obs_over_gr_minus_one": float(system.gamma_obs_s / np.sum(weights * gamma0) - 1.0),
    }
    samples["weights"] = weights
    return likelihood, meta, samples


def joint_clock_likelihood(
    systems: list[PulsarClockSystem],
    z1: np.ndarray,
    z2: np.ndarray,
    config: ClockConfig,
) -> tuple[np.ndarray, dict[str, dict[str, float]], dict[str, dict[str, np.ndarray]]]:
    rng = np.random.default_rng(config.rng_seed)
    like = np.ones((z1.size, z2.size), dtype=float)
    meta: dict[str, dict[str, float]] = {}
    cached_samples: dict[str, dict[str, np.ndarray]] = {}
    for system in systems:
        sys_like, sys_meta, sys_samples = system_likelihood_grid(system, z1, z2, config, rng)
        like *= sys_like
        meta[system.name] = sys_meta
        cached_samples[system.name] = sys_samples
    return like, meta, cached_samples


def load_tied_prior(summary_path: str | Path, bound_name: str = "optimistic") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    summary = json.loads(Path(summary_path).read_text())
    z1 = np.array(summary["sigma1_grid"], dtype=float)
    z2 = np.array(summary["sigma2_grid"], dtype=float)
    prior = np.array(summary["bounds"][bound_name]["wide"]["posterior"], dtype=float)
    return z1, z2, normalize(prior)


def evidence_under_prior(prior: np.ndarray, likelihood: np.ndarray) -> float:
    return float(np.sum(prior * likelihood))


def posterior_from_prior(prior: np.ndarray, likelihood: np.ndarray) -> np.ndarray:
    return normalize(prior * likelihood)


def posterior_mean_model_pull(
    system: PulsarClockSystem,
    samples: dict[str, np.ndarray],
    zeta1: float,
    zeta2: float,
    observable_model: str,
) -> dict[str, float]:
    pb = samples["pb_days"]
    eccentricity = samples["eccentricity"]
    mp = samples["mp_msun"]
    mc = samples["mc_msun"]
    sg = samples["self_gravity_fraction"]
    weights = samples["weights"]
    gamma0 = gamma_gr(pb, eccentricity, mp, mc)
    eta = zeta1 * sg + zeta2 * sg**2
    if observable_model == "potential_only":
        xc = mc / (mp + mc)
        gamma_model = gamma0 * (1.0 + eta / (1.0 + xc))
    elif observable_model == "full_rescale":
        gamma_model = gamma0 * (1.0 + eta)
    else:
        raise ValueError(f"unknown observable_model={observable_model!r}")

    gamma_model_mean = float(np.sum(weights * gamma_model))
    gamma_model_std = float(np.sqrt(np.sum(weights * (gamma_model - gamma_model_mean) ** 2)))
    combined_sigma = math.sqrt(system.gamma_sigma_s**2 + gamma_model_std**2)
    pull = (system.gamma_obs_s - gamma_model_mean) / combined_sigma
    return {
        "name": system.name,
        "gamma_obs_s": system.gamma_obs_s,
        "gamma_model_s": gamma_model_mean,
        "gamma_model_std_s": gamma_model_std,
        "gamma_gr_s": float(np.sum(weights * gamma0)),
        "pull": float(pull),
        "mp_mean_msun": float(np.sum(weights * mp)),
        "mc_mean_msun": float(np.sum(weights * mc)),
        "self_gravity_mean": float(np.sum(weights * sg)),
    }


def write_posterior_table(path: str | Path, z1: np.ndarray, z2: np.ndarray, posterior: np.ndarray) -> Path:
    output_path = Path(path)
    rows = []
    for i, a in enumerate(z1):
        for j, b in enumerate(z2):
            rows.append([a, b, posterior[i, j]])
    np.savetxt(output_path, np.array(rows), header="zeta1\tzeta2\tposterior", comments="", delimiter="\t")
    return output_path


def _svg_heatmap(parts: list[str], x0: float, y0: float, w: float, h: float, xgrid: np.ndarray, ygrid: np.ndarray, posterior: np.ndarray, title: str) -> None:
    parts.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" fill="none" stroke="#111827" stroke-width="1.2"/>')
    parts.append(f'<text x="{x0:.1f}" y="{y0 - 10:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">{title}</text>')
    pnorm = posterior / np.max(posterior)
    for i in range(xgrid.size - 1):
        for j in range(ygrid.size - 1):
            frac = float(pnorm[i, j])
            color = f"rgb({255 - int(110*frac)},{246 - int(90*frac)},{255 - int(170*frac)})"
            x = x0 + w * i / (xgrid.size - 1)
            y = y0 + h * (1.0 - (j + 1) / (ygrid.size - 1))
            rw = w / (xgrid.size - 1)
            rh = h / (ygrid.size - 1)
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{rw:.2f}" height="{rh:.2f}" fill="{color}" stroke="none"/>')
    for frac, label in [(0.0, f"{xgrid[0]:.1e}"), (0.5, "0"), (1.0, f"{xgrid[-1]:.1e}")]:
        x = x0 + frac * w
        parts.append(f'<line x1="{x:.1f}" y1="{y0+h:.1f}" x2="{x:.1f}" y2="{y0+h+6:.1f}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x-18:.1f}" y="{y0+h+20:.1f}" font-family="monospace" font-size="11" fill="#374151">{label}</text>')
    for frac, label in [(0.0, f"{ygrid[0]:.1e}"), (0.5, "0"), (1.0, f"{ygrid[-1]:.1e}")]:
        y = y0 + h * (1.0 - frac)
        parts.append(f'<line x1="{x0-6:.1f}" y1="{y:.1f}" x2="{x0:.1f}" y2="{y:.1f}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x0-72:.1f}" y="{y+4:.1f}" font-family="monospace" font-size="11" fill="#374151">{label}</text>')
    parts.append(f'<text x="{x0+w/2-18:.1f}" y="{y0+h+38:.1f}" font-family="sans-serif" font-size="12" fill="#111827">zeta1</text>')
    parts.append(f'<text x="{x0-60:.1f}" y="{y0+h/2:.1f}" font-family="sans-serif" font-size="12" fill="#111827" transform="rotate(-90 {x0-60:.1f},{y0+h/2:.1f})">zeta2</text>')


def _svg_pull_panel(parts: list[str], x0: float, y0: float, w: float, h: float, pulls: list[dict[str, float]], title: str) -> None:
    parts.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" fill="none" stroke="#111827" stroke-width="1.2"/>')
    parts.append(f'<text x="{x0:.1f}" y="{y0 - 10:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">{title}</text>')

    y_mid = y0 + h / 2.0

    def map_y(v: float) -> float:
        return y_mid - (h / 6.0) * v

    for val in [-2, -1, 0, 1, 2]:
        y = map_y(val)
        parts.append(f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x0+w:.1f}" y2="{y:.1f}" stroke="#d1d5db" stroke-width="1"/>')
        parts.append(f'<text x="{x0-28:.1f}" y="{y+4:.1f}" font-family="monospace" font-size="11" fill="#374151">{val}</text>')

    n = len(pulls)
    for i, item in enumerate(pulls):
        x = x0 + w * (i + 0.5) / n
        y = map_y(float(item["pull"]))
        parts.append(f'<line x1="{x:.1f}" y1="{y_mid:.1f}" x2="{x:.1f}" y2="{y:.1f}" stroke="#2563eb" stroke-width="8" stroke-linecap="round"/>')
        parts.append(f'<text x="{x-48:.1f}" y="{y0+h+18:.1f}" font-family="sans-serif" font-size="11" fill="#111827">{item["name"]}</text>')


def write_summary_svg(path: str | Path, summary: dict[str, object]) -> Path:
    output_path = Path(path)
    dec = summary["models"]["clock_only"]
    tied = summary["models"]["tied_optimistic"]
    eta_global = summary["effective_combinations"]["global_sbar"]["clock_only"]
    width, height = 1320, 920
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 6: Corrected Clock-Sector Leave-One-Out Gamma Fit</text>',
        '<text x="48" y="72" font-family="sans-serif" font-size="14" fill="#4b5563">The observable now matches the potential-only clock EFT, and masses are inferred from non-gamma PK information only.</text>',
    ]
    _svg_heatmap(parts, 90, 130, 430, 330, np.array(summary["decoupled_grid"]["zeta1"]), np.array(summary["decoupled_grid"]["zeta2"]), np.array(dec["posterior"]), "Clock-only posterior")
    _svg_heatmap(parts, 770, 130, 430, 330, np.array(summary["tied_grid"]["zeta1"]), np.array(summary["tied_grid"]["zeta2"]), np.array(tied["posterior"]), "Tied posterior (Request 5 prior)")
    _svg_pull_panel(parts, 140, 560, 1040, 220, dec["pulls"], "System-by-system pulls at the clock-only posterior mean")
    parts.append(f'<text x="90" y="830" font-family="monospace" font-size="13" fill="#111827">eta(sbar={eta_global["sbar"]:.3f}) = {eta_global["mean"]:.3e} +- {eta_global["std"]:.3e}, |eta|_95 = {eta_global["abs_95"]:.3e}</text>')
    parts.append(f'<text x="90" y="852" font-family="monospace" font-size="13" fill="#111827">Bayes factors: B(clock/tied_opt) = {summary["bayes_factors"]["clock_over_tied_optimistic"]:.3e}, B(clock/tied_cons) = {summary["bayes_factors"]["clock_over_tied_conservative"]:.3e}</text>')
    parts.append("</svg>")
    output_path.write_text("\n".join(parts))
    return output_path


def run_request6(output_dir: str | Path = ".") -> dict[str, object]:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    config = ClockConfig()
    z1_dec, z2_dec = zeta_grids(config)
    like_dec, system_meta, dec_samples = joint_clock_likelihood(SYSTEMS, z1_dec, z2_dec, config)
    prior_dec = np.ones_like(like_dec) / like_dec.size
    evidence_dec = evidence_under_prior(prior_dec, like_dec)
    posterior_dec = posterior_from_prior(prior_dec, like_dec)
    stats_dec = summarize_posterior(z1_dec, z2_dec, posterior_dec)

    req5_summary = output_root / "request5_j0337_phaseA_summary.json"
    z1_tied_opt, z2_tied_opt, prior_tied_opt = load_tied_prior(req5_summary, "optimistic")
    like_tied_opt, _, tied_opt_samples = joint_clock_likelihood(SYSTEMS, z1_tied_opt, z2_tied_opt, config)
    evidence_tied_opt = evidence_under_prior(prior_tied_opt, like_tied_opt)
    posterior_tied_opt = posterior_from_prior(prior_tied_opt, like_tied_opt)
    stats_tied_opt = summarize_posterior(z1_tied_opt, z2_tied_opt, posterior_tied_opt)

    z1_tied_cons, z2_tied_cons, prior_tied_cons = load_tied_prior(req5_summary, "conservative")
    like_tied_cons, _, _ = joint_clock_likelihood(SYSTEMS, z1_tied_cons, z2_tied_cons, config)
    evidence_tied_cons = evidence_under_prior(prior_tied_cons, like_tied_cons)
    posterior_tied_cons = posterior_from_prior(prior_tied_cons, like_tied_cons)
    stats_tied_cons = summarize_posterior(z1_tied_cons, z2_tied_cons, posterior_tied_cons)

    pulls_dec = [
        posterior_mean_model_pull(system, dec_samples[system.name], stats_dec["zeta1"]["mean"], stats_dec["zeta2"]["mean"], config.observable_model)
        for system in SYSTEMS
    ]
    pulls_tied = [
        posterior_mean_model_pull(system, tied_opt_samples[system.name], stats_tied_opt["zeta1"]["mean"], stats_tied_opt["zeta2"]["mean"], config.observable_model)
        for system in SYSTEMS
    ]

    system_sbars = {system.name: system_meta[system.name]["self_gravity_mean"] for system in SYSTEMS}
    system_xc = {system.name: system_meta[system.name]["xc_mean"] for system in SYSTEMS}
    global_sbar = float(np.mean(list(system_sbars.values())))
    sbar_spread = float(max(system_sbars.values()) - min(system_sbars.values()))
    basis_clock = basis_stats_at_sstar(z1_dec, z2_dec, posterior_dec, global_sbar)
    basis_tied_opt = basis_stats_at_sstar(z1_tied_opt, z2_tied_opt, posterior_tied_opt, global_sbar)
    basis_tied_cons = basis_stats_at_sstar(z1_tied_cons, z2_tied_cons, posterior_tied_cons, global_sbar)
    delta_gamma_per_system = {
        system.name: {
            "clock_only": delta_gamma_stats_from_eta(
                eta_stats_at_sbar(z1_dec, z2_dec, posterior_dec, system_sbars[system.name]),
                system_xc[system.name],
            ),
            "tied_optimistic": delta_gamma_stats_from_eta(
                eta_stats_at_sbar(z1_tied_opt, z2_tied_opt, posterior_tied_opt, system_sbars[system.name]),
                system_xc[system.name],
            ),
            "tied_conservative": delta_gamma_stats_from_eta(
                eta_stats_at_sbar(z1_tied_cons, z2_tied_cons, posterior_tied_cons, system_sbars[system.name]),
                system_xc[system.name],
            ),
        }
        for system in SYSTEMS
    }
    gamma_weights = np.array([1.0 / system.gamma_sigma_s**2 for system in SYSTEMS], dtype=float)
    gamma_weights /= np.sum(gamma_weights)
    slope_information_proxy = float(
        np.sum(
            gamma_weights
            * np.array(
                [
                    (system_sbars[system.name] - global_sbar) ** 2 / (1.0 + system_xc[system.name]) ** 2
                    for system in SYSTEMS
                ]
            )
        )
    )
    effective_combinations = {
        "global_sbar": {
            "clock_only": eta_stats_at_sbar(z1_dec, z2_dec, posterior_dec, global_sbar),
            "tied_optimistic": eta_stats_at_sbar(z1_tied_opt, z2_tied_opt, posterior_tied_opt, global_sbar),
            "tied_conservative": eta_stats_at_sbar(z1_tied_cons, z2_tied_cons, posterior_tied_cons, global_sbar),
        },
        "linearized_basis": {
            "clock_only": basis_clock,
            "tied_optimistic": basis_tied_opt,
            "tied_conservative": basis_tied_cons,
        },
        "delta_gamma_per_system": delta_gamma_per_system,
        "per_system_sbar": {
            system.name: {
                "clock_only": eta_stats_at_sbar(z1_dec, z2_dec, posterior_dec, system_sbars[system.name]),
                "tied_optimistic": eta_stats_at_sbar(z1_tied_opt, z2_tied_opt, posterior_tied_opt, system_sbars[system.name]),
                "tied_conservative": eta_stats_at_sbar(z1_tied_cons, z2_tied_cons, posterior_tied_cons, system_sbars[system.name]),
            }
            for system in SYSTEMS
        },
        "sbar_diagnostics": {
            "global_sbar": global_sbar,
            "sbar_spread": sbar_spread,
            "relative_spread": float(sbar_spread / global_sbar) if global_sbar != 0.0 else 0.0,
            "slope_information_proxy": slope_information_proxy,
        },
    }

    write_posterior_table(output_root / "request6_clock_sector_posterior_clock_only.tsv", z1_dec, z2_dec, posterior_dec)
    write_posterior_table(output_root / "request6_clock_sector_posterior_tied_optimistic.tsv", z1_tied_opt, z2_tied_opt, posterior_tied_opt)
    write_posterior_table(output_root / "request6_clock_sector_posterior_tied_conservative.tsv", z1_tied_cons, z2_tied_cons, posterior_tied_cons)

    summary = {
        "config": asdict(config),
        "systems": [asdict(system) for system in SYSTEMS],
        "system_meta": system_meta,
        "decoupled_grid": {"zeta1": z1_dec.tolist(), "zeta2": z2_dec.tolist()},
        "tied_grid": {"zeta1": z1_tied_opt.tolist(), "zeta2": z2_tied_opt.tolist()},
        "effective_combinations": effective_combinations,
        "models": {
            "clock_only": {
                "stats": stats_dec,
                "evidence": evidence_dec,
                "posterior": posterior_dec.tolist(),
                "pulls": pulls_dec,
            },
            "tied_optimistic": {
                "stats": stats_tied_opt,
                "evidence": evidence_tied_opt,
                "posterior": posterior_tied_opt.tolist(),
                "pulls": pulls_tied,
            },
            "tied_conservative": {
                "stats": stats_tied_cons,
                "evidence": evidence_tied_cons,
                "posterior": posterior_tied_cons.tolist(),
            },
        },
        "bayes_factors": {
            "clock_over_tied_optimistic": evidence_dec / evidence_tied_opt,
            "clock_over_tied_conservative": evidence_dec / evidence_tied_cons,
        },
    }

    summary_path = output_root / "request6_clock_sector_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    svg_path = write_summary_svg(output_root / "request6_clock_sector_summary.svg", summary)
    return {"summary": summary, "summary_path": summary_path, "svg_path": svg_path}


def print_summary(summary: dict[str, object]) -> None:
    print("=== Observable model ===")
    print(summary["config"]["observable_model"])
    print()
    print("=== Clock-only leave-one-out model ===")
    stats = summary["models"]["clock_only"]["stats"]
    print(
        f"zeta1 = {stats['zeta1']['mean']:.3e} ± {stats['zeta1']['std']:.3e}, "
        f"|zeta1|_95 = {stats['zeta1']['abs_95']:.3e}"
    )
    print(
        f"zeta2 = {stats['zeta2']['mean']:.3e} ± {stats['zeta2']['std']:.3e}, "
        f"|zeta2|_95 = {stats['zeta2']['abs_95']:.3e}"
    )
    print()
    eta_global = summary["effective_combinations"]["global_sbar"]["clock_only"]
    sbar_diag = summary["effective_combinations"]["sbar_diagnostics"]
    print("=== Effective combination actually seen by the current sample ===")
    print(
        f"sbar(global) = {sbar_diag['global_sbar']:.3f}, "
        f"spread = {sbar_diag['sbar_spread']:.3e} "
        f"(relative {sbar_diag['relative_spread']:.3%})"
    )
    print(
        f"eta(sbar) = {eta_global['mean']:.3e} ± {eta_global['std']:.3e}, "
        f"|eta(sbar)|_95 = {eta_global['abs_95']:.3e}"
    )
    basis = summary["effective_combinations"]["linearized_basis"]["clock_only"]
    print(
        f"eta_* = {basis['eta_star']['mean']:.3e} ± {basis['eta_star']['std']:.3e}, "
        f"|eta_*|_95 = {basis['eta_star']['abs_95']:.3e}"
    )
    print(
        f"kappa_* = {basis['kappa_star']['mean']:.3e} ± {basis['kappa_star']['std']:.3e}, "
        f"|kappa_*|_95 = {basis['kappa_star']['abs_95']:.3e}"
    )
    print()
    print("=== System GR leave-one-out diagnostics ===")
    for name, meta in summary["system_meta"].items():
        print(
            f"{name}: "
            f"mp = {meta['mp_mean_msun']:.6f} ± {meta['mp_std_msun']:.6f} Msun, "
            f"mc = {meta['mc_mean_msun']:.6f} ± {meta['mc_std_msun']:.6f} Msun, "
            f"Xc = {meta['xc_mean']:.3f}, "
            f"gamma_obs/gamma_GR - 1 = {meta['gamma_obs_over_gr_minus_one']:.3e}"
        )
    print()
    print("=== Tied model (zeta = sigma) ===")
    for name in ["tied_optimistic", "tied_conservative"]:
        stats = summary["models"][name]["stats"]
        print(
            f"{name}: "
            f"|zeta1|_95 = {stats['zeta1']['abs_95']:.3e}, "
            f"|zeta2|_95 = {stats['zeta2']['abs_95']:.3e}"
        )
    print()
    print("=== Bayes factors ===")
    print(
        f"B(clock/tied_opt) = {summary['bayes_factors']['clock_over_tied_optimistic']:.3e}, "
        f"B(clock/tied_cons) = {summary['bayes_factors']['clock_over_tied_conservative']:.3e}"
    )
    print()
    print("=== Approximate compressed delta_gamma by system ===")
    for system_name, stats in summary["effective_combinations"]["delta_gamma_per_system"].items():
        delta = stats["clock_only"]
        print(
            f"{system_name}: delta_gamma ~= {delta['mean']:.3e} ± {delta['std']:.3e}, "
            f"|delta_gamma|_95 = {delta['abs_95']:.3e}"
        )
    print()
    print("=== System pulls (clock-only posterior mean) ===")
    for item in summary["models"]["clock_only"]["pulls"]:
        print(
            f"{item['name']}: pull = {item['pull']:.3f}, "
            f"gamma_obs/gamma_GR - 1 = {(item['gamma_obs_s']/item['gamma_gr_s'] - 1.0):.3e}, "
            f"self-gravity mean = {item['self_gravity_mean']:.3f}"
        )


def main(output_dir: str | Path = ".") -> None:
    result = run_request6(output_dir)
    print_summary(result["summary"])
    print()
    print(f"wrote {result['summary_path']}")
    print(f"wrote {result['svg_path']}")


if __name__ == "__main__":
    main()
