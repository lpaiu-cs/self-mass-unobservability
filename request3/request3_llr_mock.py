from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


G = 6.67430e-11
C = 299_792_458.0
DAY = 86_400.0
YEAR = 365.25 * DAY
AU = 149_597_870_700.0

M_SUN = 1.98847e30
M_EARTH = 5.9722e24
M_MOON = 7.34767309e22
M_EM = M_EARTH + M_MOON

A_EMB = AU
A_MOON = 384_400_000.0

SIDEREAL_MONTH = 27.321661 * DAY
SYNODIC_MONTH = 29.530588853 * DAY

S_SUN = 3.52e-6
S_EARTH = 4.64e-10
S_MOON = 1.88e-11

SIGMA1_TEMPLATE_SCALE = 1.0
SIGMA2_TEMPLATE_SCALE = 1.0e7

ANALYTIC_SIGMA1_COEFF_M = 13.1
ANALYTIC_SIGMA2_COEFF_M = ANALYTIC_SIGMA1_COEFF_M * (S_EARTH**2 - S_MOON**2) / (S_EARTH - S_MOON)


@dataclass(frozen=True)
class LLRConfig:
    dt: float = 6.0 * 3600.0
    duration_days: int = 720
    output_step_days: int = 1
    range_noise_sigma_m: float = 0.005
    rng_seed: int = 20260422


@dataclass(frozen=True)
class InjectionCase:
    name: str
    sigma1: float
    sigma2: float
    bias_m: float
    srp_m: float
    thermal_m: float
    love_m: float
    init_dx_m: float
    init_dy_m: float
    init_dvx_mps: float
    init_dvy_mps: float


@dataclass
class SimulationResult:
    times: np.ndarray
    ranges: np.ndarray
    positions: np.ndarray
    velocities: np.ndarray


@dataclass
class ScenarioResult:
    name: str
    injected_sigma1: float
    injected_sigma2: float
    recovered_sigma1_mean: float
    recovered_sigma1_std: float
    recovered_sigma2_mean: float
    recovered_sigma2_std: float
    map_sigma1: float
    map_sigma2: float
    analytic_cosD_amplitude_m: float
    injected_cosD_amplitude_m: float
    fit_cosD_amplitude_m: float
    max_abs_residual_mm: float


def passive_factors(sigma1: float, sigma2: float) -> np.ndarray:
    s = np.array([S_SUN, S_EARTH, S_MOON], dtype=float)
    return 1.0 + sigma1 * s + sigma2 * s**2


def initial_state() -> tuple[np.ndarray, np.ndarray]:
    positions = np.zeros((3, 3), dtype=float)
    velocities = np.zeros((3, 3), dtype=float)

    n_emb = math.sqrt(G * (M_SUN + M_EM) / A_EMB**3)
    v_emb = n_emb * A_EMB
    r_emb = np.array([A_EMB, 0.0, 0.0], dtype=float)
    v_emb_vec = np.array([0.0, v_emb, 0.0], dtype=float)

    n_moon = math.sqrt(G * M_EM / A_MOON**3)
    v_rel = n_moon * A_MOON
    r_earth_rel = np.array([-M_MOON / M_EM * A_MOON, 0.0, 0.0], dtype=float)
    r_moon_rel = np.array([M_EARTH / M_EM * A_MOON, 0.0, 0.0], dtype=float)
    v_earth_rel = np.array([0.0, -M_MOON / M_EM * v_rel, 0.0], dtype=float)
    v_moon_rel = np.array([0.0, M_EARTH / M_EM * v_rel, 0.0], dtype=float)

    positions[1] = r_emb + r_earth_rel
    positions[2] = r_emb + r_moon_rel
    velocities[1] = v_emb_vec + v_earth_rel
    velocities[2] = v_emb_vec + v_moon_rel

    positions[0] = -(M_EARTH * positions[1] + M_MOON * positions[2]) / M_SUN
    velocities[0] = -(M_EARTH * velocities[1] + M_MOON * velocities[2]) / M_SUN
    return positions, velocities


def apply_moon_relative_perturbation(
    positions: np.ndarray,
    velocities: np.ndarray,
    dx: float = 0.0,
    dy: float = 0.0,
    dvx: float = 0.0,
    dvy: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    pos = positions.copy()
    vel = velocities.copy()

    pos_delta = np.array([dx, dy, 0.0], dtype=float)
    vel_delta = np.array([dvx, dvy, 0.0], dtype=float)

    pos[2] += M_EARTH / M_EM * pos_delta
    pos[1] -= M_MOON / M_EM * pos_delta
    vel[2] += M_EARTH / M_EM * vel_delta
    vel[1] -= M_MOON / M_EM * vel_delta

    pos[0] = -(M_EARTH * pos[1] + M_MOON * pos[2]) / M_SUN
    vel[0] = -(M_EARTH * vel[1] + M_MOON * vel[2]) / M_SUN
    return pos, vel


def accelerations(positions: np.ndarray, velocities: np.ndarray, passive: np.ndarray) -> np.ndarray:
    masses = np.array([M_SUN, M_EARTH, M_MOON], dtype=float)
    acc = np.zeros_like(positions)
    for i in range(3):
        for j in range(i + 1, 3):
            rij = positions[i] - positions[j]
            vij = velocities[i] - velocities[j]
            dist = np.linalg.norm(rij)
            dist3 = dist**3

            acc[i] += -passive[i] * G * masses[j] * rij / dist3
            acc[j] += passive[j] * G * masses[i] * rij / dist3

            common = 4.0 * G * (masses[i] + masses[j]) / dist - np.dot(vij, vij)
            rv = np.dot(rij, vij)
            pn_vec = common * rij + 4.0 * rv * vij
            pn_scale = 1.0 / (C**2 * dist3)

            acc[i] += passive[i] * G * masses[j] * pn_scale * pn_vec
            acc[j] += -passive[j] * G * masses[i] * pn_scale * pn_vec
    return acc


def rk4_step(positions: np.ndarray, velocities: np.ndarray, dt: float, passive: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    def rhs(pos: np.ndarray, vel: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return vel, accelerations(pos, vel, passive)

    k1_r, k1_v = rhs(positions, velocities)
    k2_r, k2_v = rhs(positions + 0.5 * dt * k1_r, velocities + 0.5 * dt * k1_v)
    k3_r, k3_v = rhs(positions + 0.5 * dt * k2_r, velocities + 0.5 * dt * k2_v)
    k4_r, k4_v = rhs(positions + dt * k3_r, velocities + dt * k3_v)

    positions_next = positions + (dt / 6.0) * (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r)
    velocities_next = velocities + (dt / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v)
    return positions_next, velocities_next


def simulate_system(
    config: LLRConfig,
    sigma1: float = 0.0,
    sigma2: float = 0.0,
    init_dx: float = 0.0,
    init_dy: float = 0.0,
    init_dvx: float = 0.0,
    init_dvy: float = 0.0,
) -> SimulationResult:
    passive = passive_factors(sigma1, sigma2)
    positions, velocities = initial_state()
    positions, velocities = apply_moon_relative_perturbation(
        positions,
        velocities,
        dx=init_dx,
        dy=init_dy,
        dvx=init_dvx,
        dvy=init_dvy,
    )

    total_steps = int(config.duration_days * DAY / config.dt)
    output_every = int(config.output_step_days * DAY / config.dt)
    num_outputs = total_steps // output_every + 1

    times = np.zeros(num_outputs, dtype=float)
    ranges = np.zeros(num_outputs, dtype=float)
    pos_out = np.zeros((num_outputs, 3, 3), dtype=float)
    vel_out = np.zeros((num_outputs, 3, 3), dtype=float)

    out = 0
    times[out] = 0.0
    ranges[out] = np.linalg.norm(positions[2] - positions[1])
    pos_out[out] = positions
    vel_out[out] = velocities
    out += 1

    for step in range(1, total_steps + 1):
        positions, velocities = rk4_step(positions, velocities, config.dt, passive)
        if step % output_every == 0:
            times[out] = step * config.dt
            ranges[out] = np.linalg.norm(positions[2] - positions[1])
            pos_out[out] = positions
            vel_out[out] = velocities
            out += 1

    return SimulationResult(times=times, ranges=ranges, positions=pos_out, velocities=vel_out)


def synodic_phase(times: np.ndarray) -> np.ndarray:
    return 2.0 * np.pi * times / SYNODIC_MONTH


def sidereal_phase(times: np.ndarray) -> np.ndarray:
    return 2.0 * np.pi * times / SIDEREAL_MONTH


def nuisance_basis(times: np.ndarray) -> dict[str, np.ndarray]:
    d = synodic_phase(times)
    l = sidereal_phase(times)
    y = 2.0 * np.pi * times / YEAR
    return {
        "bias": np.ones_like(times),
        "srp": np.cos(d + 0.35),
        "thermal": np.cos(y + 0.7),
        "love": np.cos(2.0 * l + 0.2),
    }


def fit_synodic_amplitude(times: np.ndarray, series: np.ndarray) -> dict[str, float]:
    d = synodic_phase(times)
    design = np.column_stack([np.ones_like(times), np.cos(d), np.sin(d)])
    coeff, _, _, _ = np.linalg.lstsq(design, series, rcond=None)
    return {
        "offset": float(coeff[0]),
        "cos_coeff": float(coeff[1]),
        "sin_coeff": float(coeff[2]),
        "amplitude": float(np.hypot(coeff[1], coeff[2])),
    }


def linear_nuisance_templates(config: LLRConfig, baseline: SimulationResult) -> dict[str, np.ndarray]:
    templates = {}
    perturbations = {
        "init_dx": {"init_dx": 100.0},
        "init_dy": {"init_dy": 100.0},
        "init_dvx": {"init_dvx": 1.0e-3},
        "init_dvy": {"init_dvy": 1.0e-3},
    }
    for name, kwargs in perturbations.items():
        shifted = simulate_system(config, **kwargs)
        scale = next(iter(kwargs.values()))
        templates[name] = (shifted.ranges - baseline.ranges) / scale
    return templates


def build_mock_pipeline(config: LLRConfig) -> dict[str, np.ndarray | SimulationResult | dict[str, np.ndarray]]:
    baseline = simulate_system(config)
    sigma1_template_run = simulate_system(config, sigma1=SIGMA1_TEMPLATE_SCALE)
    sigma2_template_run = simulate_system(config, sigma2=SIGMA2_TEMPLATE_SCALE)

    sigma1_template = (sigma1_template_run.ranges - baseline.ranges) / SIGMA1_TEMPLATE_SCALE
    sigma2_template = (sigma2_template_run.ranges - baseline.ranges) / SIGMA2_TEMPLATE_SCALE

    obs_basis = nuisance_basis(baseline.times)
    init_templates = linear_nuisance_templates(config, baseline)

    return {
        "baseline": baseline,
        "sigma1_template": sigma1_template,
        "sigma2_template": sigma2_template,
        "obs_basis": obs_basis,
        "init_templates": init_templates,
    }


def build_observed_range(
    config: LLRConfig,
    case: InjectionCase,
    baseline: SimulationResult,
    init_templates: dict[str, np.ndarray],
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    injected = simulate_system(config, sigma1=case.sigma1, sigma2=case.sigma2)
    basis = nuisance_basis(baseline.times)
    obs = injected.ranges.copy()
    obs += case.bias_m * basis["bias"]
    obs += case.srp_m * basis["srp"]
    obs += case.thermal_m * basis["thermal"]
    obs += case.love_m * basis["love"]
    obs += case.init_dx_m * init_templates["init_dx"]
    obs += case.init_dy_m * init_templates["init_dy"]
    obs += case.init_dvx_mps * init_templates["init_dvx"]
    obs += case.init_dvy_mps * init_templates["init_dvy"]
    obs += rng.normal(0.0, config.range_noise_sigma_m, size=obs.shape)
    return obs, injected.ranges


def nuisance_design_matrix(obs_basis: dict[str, np.ndarray], init_templates: dict[str, np.ndarray]) -> tuple[np.ndarray, list[str]]:
    names = ["bias", "srp", "thermal", "love", "init_dx", "init_dy", "init_dvx", "init_dvy"]
    matrix = np.column_stack(
        [
            obs_basis["bias"],
            obs_basis["srp"],
            obs_basis["thermal"],
            obs_basis["love"],
            init_templates["init_dx"],
            init_templates["init_dy"],
            init_templates["init_dvx"],
            init_templates["init_dvy"],
        ]
    )
    return matrix, names


def profile_chi2(
    data_residual: np.ndarray,
    sigma1_grid: np.ndarray,
    sigma2_grid: np.ndarray,
    sigma1_template: np.ndarray,
    sigma2_template: np.ndarray,
    nuisance_matrix: np.ndarray,
    noise_sigma: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    weight = 1.0 / noise_sigma
    y = data_residual * weight
    nmat = nuisance_matrix * weight
    chi2 = np.zeros((sigma1_grid.size, sigma2_grid.size), dtype=float)
    beta_store = np.zeros((sigma1_grid.size, sigma2_grid.size, nuisance_matrix.shape[1]), dtype=float)
    for i, s1 in enumerate(sigma1_grid):
        for j, s2 in enumerate(sigma2_grid):
            model_sep = (s1 * sigma1_template + s2 * sigma2_template) * weight
            rhs = y - model_sep
            beta, _, _, _ = np.linalg.lstsq(nmat, rhs, rcond=None)
            resid = rhs - nmat @ beta
            chi2[i, j] = np.dot(resid, resid)
            beta_store[i, j] = beta

    log_post = -0.5 * (chi2 - np.min(chi2))
    posterior = np.exp(log_post)
    posterior /= np.sum(posterior)
    return chi2, posterior, beta_store


def summarize_posterior(
    sigma1_grid: np.ndarray,
    sigma2_grid: np.ndarray,
    posterior: np.ndarray,
) -> dict[str, float]:
    p1 = np.sum(posterior, axis=1)
    p2 = np.sum(posterior, axis=0)
    mean1 = float(np.sum(sigma1_grid * p1))
    mean2 = float(np.sum(sigma2_grid * p2))
    std1 = float(np.sqrt(np.sum((sigma1_grid - mean1) ** 2 * p1)))
    std2 = float(np.sqrt(np.sum((sigma2_grid - mean2) ** 2 * p2)))
    mode_sigma1 = float(sigma1_grid[int(np.argmax(p1))])
    mode_sigma2 = float(sigma2_grid[int(np.argmax(p2))])
    return {
        "sigma1_mean": mean1,
        "sigma1_std": std1,
        "sigma2_mean": mean2,
        "sigma2_std": std2,
        "mode_sigma1": mode_sigma1,
        "mode_sigma2": mode_sigma2,
    }


def fit_scenario(
    config: LLRConfig,
    pipeline: dict[str, np.ndarray | SimulationResult | dict[str, np.ndarray]],
    case: InjectionCase,
    sigma1_grid: np.ndarray,
    sigma2_grid: np.ndarray,
    rng: np.random.Generator,
) -> tuple[ScenarioResult, dict[str, np.ndarray]]:
    baseline = pipeline["baseline"]
    sigma1_template = pipeline["sigma1_template"]
    sigma2_template = pipeline["sigma2_template"]
    obs_basis = pipeline["obs_basis"]
    init_templates = pipeline["init_templates"]
    assert isinstance(baseline, SimulationResult)
    assert isinstance(obs_basis, dict)
    assert isinstance(init_templates, dict)

    observed_range, injected_dynamic = build_observed_range(config, case, baseline, init_templates, rng)
    residual_data = observed_range - baseline.ranges

    nuisance_matrix, nuisance_names = nuisance_design_matrix(obs_basis, init_templates)
    chi2, posterior, beta_store = profile_chi2(
        residual_data,
        sigma1_grid,
        sigma2_grid,
        sigma1_template,
        sigma2_template,
        nuisance_matrix,
        config.range_noise_sigma_m,
    )
    summary = summarize_posterior(sigma1_grid, sigma2_grid, posterior)

    map_i = int(np.argmin(np.abs(sigma1_grid - summary["sigma1_mean"])))
    map_j = int(np.argmin(np.abs(sigma2_grid - summary["sigma2_mean"])))
    beta_map = beta_store[map_i, map_j]

    fitted_sep = summary["sigma1_mean"] * sigma1_template + summary["sigma2_mean"] * sigma2_template
    fitted_nuisance = nuisance_matrix @ beta_map
    fitted_total = baseline.ranges + fitted_sep + fitted_nuisance
    fit_residual = observed_range - fitted_total

    injected_sep_like = injected_dynamic - baseline.ranges
    injected_amp = fit_synodic_amplitude(baseline.times, injected_sep_like)
    fit_amp = fit_synodic_amplitude(baseline.times, fitted_sep)
    analytic_amp = ANALYTIC_SIGMA1_COEFF_M * case.sigma1 + ANALYTIC_SIGMA2_COEFF_M * case.sigma2

    scenario = ScenarioResult(
        name=case.name,
        injected_sigma1=case.sigma1,
        injected_sigma2=case.sigma2,
        recovered_sigma1_mean=summary["sigma1_mean"],
        recovered_sigma1_std=summary["sigma1_std"],
        recovered_sigma2_mean=summary["sigma2_mean"],
        recovered_sigma2_std=summary["sigma2_std"],
        map_sigma1=summary["mode_sigma1"],
        map_sigma2=summary["mode_sigma2"],
        analytic_cosD_amplitude_m=analytic_amp,
        injected_cosD_amplitude_m=injected_amp["amplitude"],
        fit_cosD_amplitude_m=fit_amp["amplitude"],
        max_abs_residual_mm=1000.0 * float(np.max(np.abs(fit_residual))),
    )

    arrays = {
        "times": baseline.times,
        "observed_range": observed_range,
        "baseline_range": baseline.ranges,
        "injected_dynamic": injected_dynamic,
        "fitted_range": fitted_total,
        "fit_residual": fit_residual,
        "sigma1_grid": sigma1_grid,
        "sigma2_grid": sigma2_grid,
        "posterior": posterior,
        "chi2": chi2,
        "sigma1_template": sigma1_template,
        "sigma2_template": sigma2_template,
        "fitted_sep": fitted_sep,
        "fitted_nuisance": fitted_nuisance,
        "map_nuisance": beta_map,
        "nuisance_names": np.array(nuisance_names, dtype=object),
    }
    return scenario, arrays


def write_residual_table(path: str | Path, arrays: dict[str, np.ndarray]) -> Path:
    output_path = Path(path)
    data = np.column_stack(
        [
            arrays["times"] / DAY,
            arrays["observed_range"],
            arrays["baseline_range"],
            arrays["fitted_range"],
            arrays["fit_residual"],
        ]
    )
    header = "time_days\tobserved_range_m\tbaseline_range_m\tfitted_range_m\tfit_residual_m"
    np.savetxt(output_path, data, header=header, comments="", delimiter="\t")
    return output_path


def write_posterior_table(path: str | Path, sigma1_grid: np.ndarray, sigma2_grid: np.ndarray, posterior: np.ndarray) -> Path:
    output_path = Path(path)
    rows = []
    for i, s1 in enumerate(sigma1_grid):
        for j, s2 in enumerate(sigma2_grid):
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
    sigma1_grid: np.ndarray,
    sigma2_grid: np.ndarray,
    posterior: np.ndarray,
    title: str,
) -> None:
    parts.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" fill="none" stroke="#111827" stroke-width="1.2"/>')
    parts.append(f'<text x="{x0:.1f}" y="{y0 - 10:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">{title}</text>')
    pnorm = posterior / np.max(posterior)
    for i in range(sigma1_grid.size - 1):
        for j in range(sigma2_grid.size - 1):
            frac = float(pnorm[i, j])
            shade = int(255 - 170 * frac)
            color = f"rgb({shade},{245 - int(80*frac)},{255 - int(140*frac)})"
            x = x0 + w * i / (sigma1_grid.size - 1)
            y = y0 + h * (1.0 - (j + 1) / (sigma2_grid.size - 1))
            rw = w / (sigma1_grid.size - 1)
            rh = h / (sigma2_grid.size - 1)
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{rw:.2f}" height="{rh:.2f}" fill="{color}" stroke="none"/>')
    for frac, label in [(0.0, f"{sigma1_grid[0]:.1e}"), (0.5, "0"), (1.0, f"{sigma1_grid[-1]:.1e}")]:
        x = x0 + frac * w
        parts.append(f'<line x1="{x:.1f}" y1="{y0+h:.1f}" x2="{x:.1f}" y2="{y0+h+6:.1f}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x - 20:.1f}" y="{y0+h+20:.1f}" font-family="monospace" font-size="11" fill="#374151">{label}</text>')
    for frac, label in [(0.0, f"{sigma2_grid[0]:.1e}"), (0.5, "0"), (1.0, f"{sigma2_grid[-1]:.1e}")]:
        y = y0 + h * (1.0 - frac)
        parts.append(f'<line x1="{x0-6:.1f}" y1="{y:.1f}" x2="{x0:.1f}" y2="{y:.1f}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x0-76:.1f}" y="{y+4:.1f}" font-family="monospace" font-size="11" fill="#374151">{label}</text>')
    parts.append(f'<text x="{x0 + w/2 - 20:.1f}" y="{y0+h+40:.1f}" font-family="sans-serif" font-size="12" fill="#111827">sigma1</text>')
    parts.append(f'<text x="{x0-62:.1f}" y="{y0+h/2:.1f}" font-family="sans-serif" font-size="12" fill="#111827" transform="rotate(-90 {x0-62:.1f},{y0+h/2:.1f})">sigma2</text>')


def _svg_line_panel(
    parts: list[str],
    x0: float,
    y0: float,
    w: float,
    h: float,
    times_days: np.ndarray,
    residual_mm: np.ndarray,
    title: str,
) -> None:
    parts.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" fill="none" stroke="#111827" stroke-width="1.2"/>')
    parts.append(f'<text x="{x0:.1f}" y="{y0 - 10:.1f}" font-family="sans-serif" font-size="14" font-weight="700" fill="#111827">{title}</text>')
    x_min, x_max = float(times_days[0]), float(times_days[-1])
    y_max = max(1.0, float(np.max(np.abs(residual_mm))) * 1.1)

    def map_x(x: float) -> float:
        return x0 + w * (x - x_min) / (x_max - x_min)

    def map_y(y: float) -> float:
        return y0 + h * (0.5 - 0.5 * y / y_max)

    for frac in [-1.0, 0.0, 1.0]:
        y = map_y(frac * y_max)
        parts.append(f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x0+w:.1f}" y2="{y:.1f}" stroke="#d1d5db" stroke-width="1"/>')
        parts.append(f'<text x="{x0-46:.1f}" y="{y+4:.1f}" font-family="monospace" font-size="11" fill="#374151">{frac*y_max:.1f}</text>')
    points = " ".join(f"{map_x(float(t)):.2f},{map_y(float(r)):.2f}" for t, r in zip(times_days, residual_mm))
    parts.append(f'<polyline fill="none" stroke="#dc2626" stroke-width="1.6" points="{points}"/>')
    parts.append(f'<text x="{x0 + w/2 - 25:.1f}" y="{y0+h+38:.1f}" font-family="sans-serif" font-size="12" fill="#111827">time [days]</text>')
    parts.append(f'<text x="{x0-62:.1f}" y="{y0+h/2:.1f}" font-family="sans-serif" font-size="12" fill="#111827" transform="rotate(-90 {x0-62:.1f},{y0+h/2:.1f})">residual [mm]</text>')


def write_summary_svg(
    path: str | Path,
    injected_arrays: dict[str, np.ndarray],
    null_arrays: dict[str, np.ndarray],
) -> Path:
    width, height = 1280, 900
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="44" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 3: LLR Mock Injection-Recovery</text>',
        '<text x="48" y="70" font-family="sans-serif" font-size="14" fill="#4b5563">Toy Sun-Earth-Moon 1PN integration with linearized sigma1/sigma2 templates and nuisance marginalization.</text>',
    ]
    _svg_line_panel(
        parts,
        70,
        120,
        1140,
        230,
        injected_arrays["times"] / DAY,
        1000.0 * injected_arrays["fit_residual"],
        "Injected-case fit residual after nuisance marginalization",
    )
    _svg_heatmap(
        parts,
        120,
        450,
        420,
        320,
        injected_arrays["sigma1_grid"],
        injected_arrays["sigma2_grid"],
        injected_arrays["posterior"],
        "Injected-case posterior",
    )
    _svg_heatmap(
        parts,
        720,
        450,
        420,
        320,
        null_arrays["sigma1_grid"],
        null_arrays["sigma2_grid"],
        null_arrays["posterior"],
        "Null-case posterior",
    )
    parts.append("</svg>")
    output_path = Path(path)
    output_path.write_text("\n".join(parts))
    return output_path


def run_request3(output_dir: str | Path = ".") -> dict[str, object]:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    config = LLRConfig()
    rng = np.random.default_rng(config.rng_seed)
    pipeline = build_mock_pipeline(config)

    sigma1_grid = np.linspace(-8.0e-4, 8.0e-4, 81)
    sigma2_grid = np.linspace(-2.0e5, 2.0e5, 81)

    null_case = InjectionCase(
        name="null",
        sigma1=0.0,
        sigma2=0.0,
        bias_m=0.012,
        srp_m=0.010,
        thermal_m=-0.008,
        love_m=0.006,
        init_dx_m=1.8,
        init_dy_m=-1.2,
        init_dvx_mps=2.5e-5,
        init_dvy_mps=-1.5e-5,
    )
    injected_case = InjectionCase(
        name="injected",
        sigma1=5.0e-4,
        sigma2=0.0,
        bias_m=0.018,
        srp_m=0.012,
        thermal_m=-0.010,
        love_m=0.008,
        init_dx_m=2.6,
        init_dy_m=-1.5,
        init_dvx_mps=3.0e-5,
        init_dvy_mps=-2.0e-5,
    )

    null_result, null_arrays = fit_scenario(config, pipeline, null_case, sigma1_grid, sigma2_grid, rng)
    injected_result, injected_arrays = fit_scenario(config, pipeline, injected_case, sigma1_grid, sigma2_grid, rng)

    write_residual_table(output_root / "request3_llr_mock_residuals_null.tsv", null_arrays)
    write_residual_table(output_root / "request3_llr_mock_residuals_injected.tsv", injected_arrays)
    write_posterior_table(output_root / "request3_llr_mock_posterior_null.tsv", sigma1_grid, sigma2_grid, null_arrays["posterior"])
    write_posterior_table(output_root / "request3_llr_mock_posterior_injected.tsv", sigma1_grid, sigma2_grid, injected_arrays["posterior"])
    write_summary_svg(output_root / "request3_llr_mock_summary.svg", injected_arrays, null_arrays)

    summary = {
        "config": asdict(config),
        "analytic_coefficients_m": {
            "sigma1": ANALYTIC_SIGMA1_COEFF_M,
            "sigma2": ANALYTIC_SIGMA2_COEFF_M,
        },
        "scenarios": [asdict(null_result), asdict(injected_result)],
    }
    summary_path = output_root / "request3_llr_mock_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    return {
        "summary": summary,
        "summary_path": summary_path,
    }


def print_summary(summary: dict[str, object]) -> None:
    analytic = summary["analytic_coefficients_m"]
    print("=== Analytic weak-field benchmark ===")
    print(f"A_cosD ≈ ({analytic['sigma1']:.6g}) sigma1 + ({analytic['sigma2']:.6g}) sigma2 [m]")
    print()
    print("=== Scenario summary ===")
    for scenario in summary["scenarios"]:
        print(
            f"{scenario['name']}: "
            f"inj=({scenario['injected_sigma1']:.3e}, {scenario['injected_sigma2']:.3e}), "
            f"mean=({scenario['recovered_sigma1_mean']:.3e} ± {scenario['recovered_sigma1_std']:.3e}, "
            f"{scenario['recovered_sigma2_mean']:.3e} ± {scenario['recovered_sigma2_std']:.3e}), "
            f"mode=({scenario['map_sigma1']:.3e}, {scenario['map_sigma2']:.3e}), "
            f"AcosD[inj/fit/analytic]=({scenario['injected_cosD_amplitude_m']:.3e}, "
            f"{scenario['fit_cosD_amplitude_m']:.3e}, {scenario['analytic_cosD_amplitude_m']:.3e}) m, "
            f"max|resid|={scenario['max_abs_residual_mm']:.2f} mm"
        )


def main(output_dir: str | Path = ".") -> None:
    result = run_request3(output_dir)
    print_summary(result["summary"])
    print()
    print(f"wrote {result['summary_path']}")


if __name__ == "__main__":
    main()
