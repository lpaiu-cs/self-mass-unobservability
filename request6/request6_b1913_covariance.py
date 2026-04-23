from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

import request6_clock_sector as r6
import request6_lever_arm_audit as audit


@dataclass(frozen=True)
class B1913Config:
    observable_model: str = "potential_only"
    mass_samples: int = 4096
    batch_size: int = 128
    rng_seed: int = 191316
    branch_prior_heavy_pulsar: float = 0.5
    current_summary_path: str = "request6_clock_sector_summary.json"
    current_posterior_path: str = "request6_clock_sector_posterior_clock_only.tsv"
    current_source_scout_path: str = "request6_source_scout_summary.json"


@dataclass(frozen=True)
class B1913System:
    name: str = "PSR B1913+16"
    pb_days: float = 0.322997448911
    pb_sigma_days: float = 4.0e-12
    eccentricity: float = 0.6171334
    eccentricity_sigma: float = 5.0e-7
    projected_semimajor_axis_s: float = 2.341782
    projected_semimajor_axis_sigma_s: float = 3.0e-6
    periastron_advance_deg_yr: float = 4.226598
    periastron_advance_sigma_deg_yr: float = 5.0e-6
    pbdot_obs: float = -2.423e-12
    pbdot_obs_sigma: float = 1.0e-15
    pbdot_gal_mean: float = -0.027e-12
    pbdot_gal_sigma: float = 0.005e-12
    gamma_obs_s: float = 4.2992e-3
    gamma_sigma_s: float = 8.0e-7
    source_label: str = "Weisberg, Nice, Taylor (2010)"
    source_url: str = "https://arxiv.org/abs/1011.0718"
    notes: str = (
        "Covariance-aware here means a hierarchical nuisance-aware PK fit over "
        "(Pb, e, omdot, Pbdot_obs, DeltaPbdot_gal) with mass-root branch "
        "marginalization. No published full timing covariance matrix is used."
    )


def pbdot_gr(pb_days: np.ndarray, e: np.ndarray, mp_msun: np.ndarray, mc_msun: np.ndarray) -> np.ndarray:
    pb = pb_days * r6.SECONDS_PER_DAY
    mtot = mp_msun + mc_msun
    ecc_factor = (1.0 + (73.0 / 24.0) * e**2 + (37.0 / 96.0) * e**4) / np.clip((1.0 - e**2) ** 3.5, 1.0e-18, None)
    coeff = (-192.0 * math.pi / 5.0) * (r6.T_SUN ** (5.0 / 3.0)) * (pb / (2.0 * math.pi)) ** (-5.0 / 3.0) * ecc_factor
    return coeff * mp_msun * mc_msun / np.cbrt(mtot)


def mass_roots_from_total_and_pbdot(pb_days: np.ndarray, e: np.ndarray, total_mass: np.ndarray, pbdot_intrinsic: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pb = pb_days * r6.SECONDS_PER_DAY
    ecc_factor = (1.0 + (73.0 / 24.0) * e**2 + (37.0 / 96.0) * e**4) / np.clip((1.0 - e**2) ** 3.5, 1.0e-18, None)
    coeff = (-192.0 * math.pi / 5.0) * (r6.T_SUN ** (5.0 / 3.0)) * (pb / (2.0 * math.pi)) ** (-5.0 / 3.0) * ecc_factor
    pair_product = pbdot_intrinsic * np.cbrt(total_mass) / coeff
    discriminant = total_mass**2 - 4.0 * pair_product
    valid = (
        np.isfinite(pair_product)
        & np.isfinite(discriminant)
        & (pair_product > 0.0)
        & (discriminant >= 0.0)
        & (pbdot_intrinsic < 0.0)
        & (coeff < 0.0)
    )
    root = np.sqrt(np.clip(discriminant, 0.0, None))
    heavy = 0.5 * (total_mass + root)
    light = 0.5 * (total_mass - root)
    valid &= (light > 0.0) & (heavy > 0.0)
    return heavy, light, valid


def load_current_posterior_grid(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = np.loadtxt(path, skiprows=1)
    z1 = np.unique(table[:, 0])
    z2 = np.unique(table[:, 1])
    posterior = table[:, 2].reshape(z1.size, z2.size)
    return z1, z2, posterior


def sample_b1913(system: B1913System, config: B1913Config) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(config.rng_seed)
    remaining = config.mass_samples
    collected: dict[str, list[np.ndarray]] = {
        "pb_days": [],
        "eccentricity": [],
        "projected_semimajor_axis_s": [],
        "total_mass_msun": [],
        "pbdot_obs": [],
        "pbdot_gal": [],
        "pbdot_intrinsic": [],
        "m_heavy_msun": [],
        "m_light_msun": [],
        "s_heavy": [],
        "s_light": [],
        "sin_i_if_heavy_pulsar": [],
        "sin_i_if_light_pulsar": [],
    }

    while remaining > 0:
        draw_count = max(1024, 2 * remaining)
        pb = rng.normal(system.pb_days, system.pb_sigma_days, draw_count)
        e = rng.normal(system.eccentricity, system.eccentricity_sigma, draw_count)
        x = rng.normal(system.projected_semimajor_axis_s, system.projected_semimajor_axis_sigma_s, draw_count)
        omdot = rng.normal(system.periastron_advance_deg_yr, system.periastron_advance_sigma_deg_yr, draw_count)
        pbdot_obs = rng.normal(system.pbdot_obs, system.pbdot_obs_sigma, draw_count)
        pbdot_gal = rng.normal(system.pbdot_gal_mean, system.pbdot_gal_sigma, draw_count)

        total_mass = r6.total_mass_from_omdot(pb, e, omdot)
        pbdot_intrinsic = pbdot_obs - pbdot_gal
        heavy, light, valid = mass_roots_from_total_and_pbdot(pb, e, total_mass, pbdot_intrinsic)
        valid &= (pb > 0.0) & (e > 0.0) & (e < 1.0) & np.isfinite(total_mass)

        if not np.any(valid):
            continue

        take = min(remaining, int(np.sum(valid)))
        idx = np.flatnonzero(valid)[:take]
        pb = pb[idx]
        e = e[idx]
        x = x[idx]
        total_mass = total_mass[idx]
        pbdot_obs = pbdot_obs[idx]
        pbdot_gal = pbdot_gal[idx]
        pbdot_intrinsic = pbdot_intrinsic[idx]
        heavy = heavy[idx]
        light = light[idx]

        s_lo_heavy, s_hi_heavy = r6.mass_conditioned_s_prior(heavy)
        s_lo_light, s_hi_light = r6.mass_conditioned_s_prior(light)
        s_heavy = rng.uniform(s_lo_heavy, s_hi_heavy)
        s_light = rng.uniform(s_lo_light, s_hi_light)

        mf = r6.mass_function(x, pb)
        sin_i_heavy = np.cbrt(np.clip(mf * total_mass**2 / np.clip(light, 1.0e-18, None) ** 3, 0.0, None))
        sin_i_light = np.cbrt(np.clip(mf * total_mass**2 / np.clip(heavy, 1.0e-18, None) ** 3, 0.0, None))

        collected["pb_days"].append(pb)
        collected["eccentricity"].append(e)
        collected["projected_semimajor_axis_s"].append(x)
        collected["total_mass_msun"].append(total_mass)
        collected["pbdot_obs"].append(pbdot_obs)
        collected["pbdot_gal"].append(pbdot_gal)
        collected["pbdot_intrinsic"].append(pbdot_intrinsic)
        collected["m_heavy_msun"].append(heavy)
        collected["m_light_msun"].append(light)
        collected["s_heavy"].append(s_heavy)
        collected["s_light"].append(s_light)
        collected["sin_i_if_heavy_pulsar"].append(sin_i_heavy)
        collected["sin_i_if_light_pulsar"].append(sin_i_light)
        remaining -= take

    return {key: np.concatenate(value) for key, value in collected.items()}


def branch_gamma_models(
    samples: dict[str, np.ndarray],
    z1: np.ndarray,
    z2: np.ndarray,
    observable_model: str,
    start: int,
    stop: int,
) -> tuple[np.ndarray, np.ndarray]:
    pb = samples["pb_days"][start:stop]
    eccentricity = samples["eccentricity"][start:stop]

    gamma_heavy = r6.gamma_clock(
        pb,
        eccentricity,
        samples["m_heavy_msun"][start:stop],
        samples["m_light_msun"][start:stop],
        samples["s_heavy"][start:stop],
        z1,
        z2,
        observable_model,
    )
    gamma_light = r6.gamma_clock(
        pb,
        eccentricity,
        samples["m_light_msun"][start:stop],
        samples["m_heavy_msun"][start:stop],
        samples["s_light"][start:stop],
        z1,
        z2,
        observable_model,
    )
    return gamma_heavy, gamma_light


def branch_responsibilities(
    system: B1913System,
    samples: dict[str, np.ndarray],
    zeta1: float,
    zeta2: float,
    observable_model: str,
    branch_prior_heavy: float,
) -> dict[str, np.ndarray | float]:
    pb = samples["pb_days"]
    e = samples["eccentricity"]
    gamma_heavy = r6.gamma_clock(pb, e, samples["m_heavy_msun"], samples["m_light_msun"], samples["s_heavy"], np.array([zeta1]), np.array([zeta2]), observable_model)[0, 0]
    gamma_light = r6.gamma_clock(pb, e, samples["m_light_msun"], samples["m_heavy_msun"], samples["s_light"], np.array([zeta1]), np.array([zeta2]), observable_model)[0, 0]
    like_heavy = branch_prior_heavy * np.exp(-0.5 * ((system.gamma_obs_s - gamma_heavy) / system.gamma_sigma_s) ** 2)
    like_light = (1.0 - branch_prior_heavy) * np.exp(-0.5 * ((system.gamma_obs_s - gamma_light) / system.gamma_sigma_s) ** 2)
    norm = np.clip(like_heavy + like_light, 1.0e-300, None)
    resp_heavy = like_heavy / norm
    resp_light = like_light / norm
    total = float(np.sum(like_heavy) + np.sum(like_light))
    return {
        "resp_heavy": resp_heavy,
        "resp_light": resp_light,
        "like_heavy": like_heavy,
        "like_light": like_light,
        "mean_heavy": float(np.sum(like_heavy) / total),
        "mean_light": float(np.sum(like_light) / total),
        "gamma_heavy_mean_s": float(np.mean(gamma_heavy)),
        "gamma_light_mean_s": float(np.mean(gamma_light)),
    }


def b1913_likelihood_grid(
    system: B1913System,
    samples: dict[str, np.ndarray],
    z1: np.ndarray,
    z2: np.ndarray,
    config: B1913Config,
) -> np.ndarray:
    like_sum = np.zeros((z1.size, z2.size), dtype=float)
    prior_heavy = config.branch_prior_heavy_pulsar
    prior_light = 1.0 - prior_heavy
    for start in range(0, samples["pb_days"].size, config.batch_size):
        stop = min(start + config.batch_size, samples["pb_days"].size)
        model_heavy, model_light = branch_gamma_models(samples, z1, z2, config.observable_model, start, stop)
        resid_heavy = system.gamma_obs_s - model_heavy
        resid_light = system.gamma_obs_s - model_light
        like_heavy = np.exp(-0.5 * (resid_heavy / system.gamma_sigma_s) ** 2)
        like_light = np.exp(-0.5 * (resid_light / system.gamma_sigma_s) ** 2)
        like_sum += (prior_heavy * like_heavy + prior_light * like_light).sum(axis=2)
    return like_sum / samples["pb_days"].size


def summarize_branch_selected_meta(
    system: B1913System,
    samples: dict[str, np.ndarray],
    zeta1: float,
    zeta2: float,
    observable_model: str,
    branch_prior_heavy: float,
) -> dict[str, object]:
    branch = branch_responsibilities(system, samples, zeta1, zeta2, observable_model, branch_prior_heavy)
    like_heavy = np.asarray(branch["like_heavy"], dtype=float)
    like_light = np.asarray(branch["like_light"], dtype=float)
    combined_weights = np.concatenate([like_heavy, like_light])
    combined_weights /= np.sum(combined_weights)

    mp = np.concatenate([samples["m_heavy_msun"], samples["m_light_msun"]])
    mc = np.concatenate([samples["m_light_msun"], samples["m_heavy_msun"]])
    sg = np.concatenate([samples["s_heavy"], samples["s_light"]])
    xc = mc / (mp + mc)
    total_mass = np.concatenate([samples["total_mass_msun"], samples["total_mass_msun"]])
    pbdot_gal = np.concatenate([samples["pbdot_gal"], samples["pbdot_gal"]])
    pbdot_intrinsic = np.concatenate([samples["pbdot_intrinsic"], samples["pbdot_intrinsic"]])
    gamma_gr = np.concatenate(
        [
            r6.gamma_gr(samples["pb_days"], samples["eccentricity"], samples["m_heavy_msun"], samples["m_light_msun"]),
            r6.gamma_gr(samples["pb_days"], samples["eccentricity"], samples["m_light_msun"], samples["m_heavy_msun"]),
        ]
    )
    gamma_model = np.concatenate(
        [
            r6.gamma_clock(
                samples["pb_days"],
                samples["eccentricity"],
                samples["m_heavy_msun"],
                samples["m_light_msun"],
                samples["s_heavy"],
                np.array([zeta1]),
                np.array([zeta2]),
                observable_model,
            )[0, 0],
            r6.gamma_clock(
                samples["pb_days"],
                samples["eccentricity"],
                samples["m_light_msun"],
                samples["m_heavy_msun"],
                samples["s_light"],
                np.array([zeta1]),
                np.array([zeta2]),
                observable_model,
            )[0, 0],
        ]
    )

    def mean_std(values: np.ndarray) -> tuple[float, float]:
        mean = float(np.sum(combined_weights * values))
        std = float(np.sqrt(np.sum(combined_weights * (values - mean) ** 2)))
        return mean, std

    mp_mean, mp_std = mean_std(mp)
    mc_mean, mc_std = mean_std(mc)
    sg_mean, sg_std = mean_std(sg)
    xc_mean, xc_std = mean_std(xc)
    total_mass_mean, total_mass_std = mean_std(total_mass)
    pbdot_gal_mean, pbdot_gal_std = mean_std(pbdot_gal)
    pbdot_intrinsic_mean, pbdot_intrinsic_std = mean_std(pbdot_intrinsic)
    gamma_gr_mean, gamma_gr_std = mean_std(gamma_gr)
    gamma_model_mean, gamma_model_std = mean_std(gamma_model)

    cov_fields = {
        "mp_msun": mp,
        "mc_msun": mc,
        "self_gravity": sg,
        "gamma_gr_s": gamma_gr,
    }
    field_names = list(cov_fields)
    matrix = np.vstack([cov_fields[name] for name in field_names])
    covariance = np.cov(matrix, aweights=combined_weights, bias=True)

    return {
        "branch_probabilities": {
            "heavy_pulsar": float(np.sum(like_heavy) / (np.sum(like_heavy) + np.sum(like_light))),
            "light_pulsar": float(np.sum(like_light) / (np.sum(like_heavy) + np.sum(like_light))),
        },
        "mp_mean_msun": mp_mean,
        "mp_std_msun": mp_std,
        "mc_mean_msun": mc_mean,
        "mc_std_msun": mc_std,
        "self_gravity_mean": sg_mean,
        "self_gravity_std": sg_std,
        "xc_mean": xc_mean,
        "xc_std": xc_std,
        "total_mass_mean_msun": total_mass_mean,
        "total_mass_std_msun": total_mass_std,
        "pbdot_gal_mean": pbdot_gal_mean,
        "pbdot_gal_std": pbdot_gal_std,
        "pbdot_intrinsic_mean": pbdot_intrinsic_mean,
        "pbdot_intrinsic_std": pbdot_intrinsic_std,
        "gamma_pred_gr_s": gamma_gr_mean,
        "gamma_pred_gr_std_s": gamma_gr_std,
        "gamma_model_s": gamma_model_mean,
        "gamma_model_std_s": gamma_model_std,
        "gamma_obs_over_gr_minus_one": float(system.gamma_obs_s / gamma_gr_mean - 1.0),
        "covariance_fields": field_names,
        "covariance_matrix": covariance.tolist(),
    }


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(weights * values))


def _weighted_std(values: np.ndarray, weights: np.ndarray) -> float:
    mean = _weighted_mean(values, weights)
    return float(np.sqrt(np.sum(weights * (values - mean) ** 2)))


def precision_weighted_sref(rows: list[dict[str, float | str]]) -> float:
    weights = np.array([1.0 / float(row["sigma_delta"]) ** 2 for row in rows], dtype=float)
    weights /= np.sum(weights)
    return float(np.sum(weights * np.array([float(row["sbar"]) for row in rows], dtype=float)))


def make_audit_row_from_meta(name: str, meta: dict[str, object], gamma_obs_s: float, gamma_sigma_s: float) -> dict[str, float | str]:
    sigma_total = math.sqrt(gamma_sigma_s**2 + float(meta["gamma_model_std_s"]) ** 2)
    return {
        "name": name,
        "sbar": float(meta["self_gravity_mean"]),
        "xc": float(meta["xc_mean"]),
        "gamma_gr_s": float(meta["gamma_pred_gr_s"]),
        "sigma_gamma_s": sigma_total,
        "sigma_delta": float(sigma_total / float(meta["gamma_pred_gr_s"])),
        "delta_gamma_obs": float(gamma_obs_s / float(meta["gamma_pred_gr_s"]) - 1.0),
    }


def summarize_non_gamma_predictive(samples: dict[str, np.ndarray]) -> dict[str, float]:
    gamma_heavy = r6.gamma_gr(samples["pb_days"], samples["eccentricity"], samples["m_heavy_msun"], samples["m_light_msun"])
    gamma_light = r6.gamma_gr(samples["pb_days"], samples["eccentricity"], samples["m_light_msun"], samples["m_heavy_msun"])
    return {
        "heavy_branch_gamma_mean_s": float(np.mean(gamma_heavy)),
        "heavy_branch_gamma_std_s": float(np.std(gamma_heavy)),
        "light_branch_gamma_mean_s": float(np.mean(gamma_light)),
        "light_branch_gamma_std_s": float(np.std(gamma_light)),
        "total_mass_mean_msun": float(np.mean(samples["total_mass_msun"])),
        "total_mass_std_msun": float(np.std(samples["total_mass_msun"])),
        "pbdot_gal_mean": float(np.mean(samples["pbdot_gal"])),
        "pbdot_gal_std": float(np.std(samples["pbdot_gal"])),
        "pbdot_intrinsic_mean": float(np.mean(samples["pbdot_intrinsic"])),
        "pbdot_intrinsic_std": float(np.std(samples["pbdot_intrinsic"])),
    }


def write_summary_svg(path: str | Path, summary: dict[str, object]) -> Path:
    width, height = 1320, 900
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 6: B1913+16 covariance-aware clock-only extension</text>',
        '<text x="48" y="72" font-family="sans-serif" font-size="14" fill="#475569">Masses are inferred from omdot and Pbdot with a Galactic-bias nuisance, and the unordered mass roots are branch-marginalized before gamma is applied.</text>',
    ]
    r6._svg_heatmap(
        parts,
        90,
        130,
        430,
        330,
        np.array(summary["zeta_grid"]["zeta1"]),
        np.array(summary["zeta_grid"]["zeta2"]),
        np.array(summary["b1913_only"]["posterior"]),
        "B1913-only posterior",
    )
    r6._svg_heatmap(
        parts,
        770,
        130,
        430,
        330,
        np.array(summary["zeta_grid"]["zeta1"]),
        np.array(summary["zeta_grid"]["zeta2"]),
        np.array(summary["combined_with_request6"]["posterior"]),
        "Current Request 6 + B1913 posterior",
    )
    parts.append(
        f'<text x="90" y="520" font-family="monospace" font-size="13" fill="#111827">branch weight (heavy pulsar | GR gamma): {summary["gr_branch_selection"]["heavy_pulsar"]:.6f}</text>'
    )
    parts.append(
        f'<text x="90" y="544" font-family="monospace" font-size="13" fill="#111827">non-gamma predictive gamma_GR(heavy branch) = {summary["non_gamma_predictive"]["heavy_branch_gamma_mean_s"]*1.0e3:.4f} +- {summary["non_gamma_predictive"]["heavy_branch_gamma_std_s"]*1.0e3:.4f} ms</text>'
    )
    parts.append(
        f'<text x="90" y="568" font-family="monospace" font-size="13" fill="#111827">B1913-only eta(s_B1913={summary["b1913_only"]["eta_at_b1913"]["sbar"]:.3f}) = {summary["b1913_only"]["eta_at_b1913"]["mean"]:.3e} +- {summary["b1913_only"]["eta_at_b1913"]["std"]:.3e}</text>'
    )
    parts.append(
        f'<text x="90" y="592" font-family="monospace" font-size="13" fill="#111827">combined eta*(s_ref={summary["combined_with_request6"]["sref"]:.3f}) = {summary["combined_with_request6"]["basis"]["eta_star"]["mean"]:.3e} +- {summary["combined_with_request6"]["basis"]["eta_star"]["std"]:.3e}</text>'
    )
    parts.append(
        f'<text x="90" y="616" font-family="monospace" font-size="13" fill="#111827">combined |kappa_*|_95 = {summary["combined_with_request6"]["basis"]["kappa_star"]["abs_95"]:.3e}; audit update = {summary["lever_arm_update"]["with_b1913"]["abs95_kappa_star"]:.3e}</text>'
    )
    parts.append(
        f'<text x="90" y="640" font-family="monospace" font-size="13" fill="#111827">effective sigma_delta(B1913) = {summary["b1913_effective_row"]["sigma_delta"]:.3e}; source-scout target sigma_delta ~ 2e-4</text>'
    )
    parts.append('<rect x="90" y="684" width="1110" height="140" fill="white" stroke="#111827" stroke-width="1.2"/>')
    parts.append('<text x="108" y="712" font-family="sans-serif" font-size="15" font-weight="700" fill="#111827">Interpretation</text>')
    parts.append(
        '<text x="108" y="738" font-family="sans-serif" font-size="13" fill="#374151">B1913 sharply selects the published heavy-pulsar branch through gamma, and it tightens the decoupled clock-only posterior substantially.</text>'
    )
    parts.append(
        '<text x="108" y="760" font-family="sans-serif" font-size="13" fill="#374151">But this is still a clock-side strengthening step, not a clean final tied-model verdict, because the mass reconstruction itself relies on relativistic orbital dynamics plus a Galactic-bias nuisance.</text>'
    )
    parts.append(
        '<text x="108" y="782" font-family="sans-serif" font-size="13" fill="#374151">The next bottleneck remains the low-side source with explicit nuisance control, exactly as the source scout predicted.</text>'
    )
    parts.append("</svg>")
    output = Path(path)
    output.write_text("\n".join(parts))
    return output


def build_summary(output_dir: str | Path = ".") -> dict[str, object]:
    output_root = Path(output_dir)
    config = B1913Config()
    system = B1913System()

    z1, z2, current_posterior = load_current_posterior_grid(output_root / config.current_posterior_path)
    current_summary = json.loads((output_root / config.current_summary_path).read_text())
    current_rows = audit.current_rows(current_summary)
    current_sstar = float(current_summary["effective_combinations"]["linearized_basis"]["clock_only"]["sstar"])

    samples = sample_b1913(system, config)
    b1913_like = b1913_likelihood_grid(system, samples, z1, z2, config)
    prior_uniform = np.ones_like(b1913_like) / b1913_like.size
    b1913_only_posterior = r6.posterior_from_prior(prior_uniform, b1913_like)
    combined_posterior = r6.normalize(current_posterior * b1913_like)

    stats_b1913 = r6.summarize_posterior(z1, z2, b1913_only_posterior)
    stats_combined = r6.summarize_posterior(z1, z2, combined_posterior)

    non_gamma_predictive = summarize_non_gamma_predictive(samples)
    branch_gr = summarize_branch_selected_meta(system, samples, 0.0, 0.0, config.observable_model, config.branch_prior_heavy_pulsar)
    branch_b1913_mean = summarize_branch_selected_meta(
        system,
        samples,
        stats_b1913["zeta1"]["mean"],
        stats_b1913["zeta2"]["mean"],
        config.observable_model,
        config.branch_prior_heavy_pulsar,
    )
    branch_combined_mean = summarize_branch_selected_meta(
        system,
        samples,
        stats_combined["zeta1"]["mean"],
        stats_combined["zeta2"]["mean"],
        config.observable_model,
        config.branch_prior_heavy_pulsar,
    )

    b1913_only_eta = r6.eta_stats_at_sbar(z1, z2, b1913_only_posterior, float(branch_gr["self_gravity_mean"]))
    b1913_only_basis = r6.basis_stats_at_sstar(z1, z2, b1913_only_posterior, float(branch_gr["self_gravity_mean"]))

    effective_row = make_audit_row_from_meta(system.name, branch_gr, system.gamma_obs_s, system.gamma_sigma_s)
    combined_rows = list(current_rows) + [effective_row]
    combined_sref = precision_weighted_sref(combined_rows)
    combined_basis = r6.basis_stats_at_sstar(z1, z2, combined_posterior, combined_sref)
    combined_eta = r6.eta_stats_at_sbar(z1, z2, combined_posterior, combined_sref)
    current_fisher = audit.fisher_eta_kappa(current_rows, current_sstar)
    updated_fisher = audit.fisher_eta_kappa(combined_rows, current_sstar)

    current_source_scout = json.loads((output_root / config.current_source_scout_path).read_text())
    scout_b1913 = None
    for row in current_source_scout["candidate_rows_ranked"]:
        if row["name"] == system.name:
            scout_b1913 = row
            break
    assert scout_b1913 is not None

    r6.write_posterior_table(output_root / "request6_b1913_covariance_posterior_b1913_only.tsv", z1, z2, b1913_only_posterior)
    r6.write_posterior_table(output_root / "request6_b1913_covariance_posterior_combined.tsv", z1, z2, combined_posterior)

    summary = {
        "config": asdict(config),
        "system": asdict(system),
        "zeta_grid": {"zeta1": z1.tolist(), "zeta2": z2.tolist()},
        "sample_count": int(samples["pb_days"].size),
        "non_gamma_predictive": non_gamma_predictive,
        "current_request6_reference": {
            "sstar": current_sstar,
            "basis": current_summary["effective_combinations"]["linearized_basis"]["clock_only"],
        },
        "gr_branch_selection": branch_gr["branch_probabilities"],
        "source_scout_reference": scout_b1913,
        "b1913_effective_row": effective_row,
        "b1913_only": {
            "stats": stats_b1913,
            "posterior": b1913_only_posterior.tolist(),
            "eta_at_b1913": b1913_only_eta,
            "basis": b1913_only_basis,
            "branch_selected_meta_at_mean": branch_b1913_mean,
        },
        "combined_with_request6": {
            "stats": stats_combined,
            "posterior": combined_posterior.tolist(),
            "sref": combined_sref,
            "eta_at_sref": combined_eta,
            "basis": combined_basis,
            "branch_selected_meta_at_mean": branch_combined_mean,
        },
        "lever_arm_update": {
            "current_only": current_fisher,
            "with_b1913": updated_fisher,
            "improvement_factor": float(current_fisher["abs95_kappa_star"] / updated_fisher["abs95_kappa_star"]),
        },
    }
    (output_root / "request6_b1913_covariance_summary.json").write_text(json.dumps(summary, indent=2))
    write_summary_svg(output_root / "request6_b1913_covariance_summary.svg", summary)
    return summary


def print_summary(summary: dict[str, object]) -> None:
    print("=== Request 6: B1913 covariance-aware extension ===")
    print(f"samples = {summary['sample_count']}")
    print(
        f"GR branch selection: heavy-pulsar = {summary['gr_branch_selection']['heavy_pulsar']:.6f}, "
        f"light-pulsar = {summary['gr_branch_selection']['light_pulsar']:.6f}"
    )
    print()
    ng = summary["non_gamma_predictive"]
    print("=== Non-gamma predictive check before gamma weighting ===")
    print(
        f"heavy-branch gamma_GR = {1.0e3*ng['heavy_branch_gamma_mean_s']:.4f} ± {1.0e3*ng['heavy_branch_gamma_std_s']:.4f} ms, "
        f"light-branch gamma_GR = {1.0e3*ng['light_branch_gamma_mean_s']:.4f} ± {1.0e3*ng['light_branch_gamma_std_s']:.4f} ms"
    )
    print(
        f"DeltaPbdot_gal prior predictive = {ng['pbdot_gal_mean']:.3e} ± {ng['pbdot_gal_std']:.3e}, "
        f"Pbdot_intrinsic = {ng['pbdot_intrinsic_mean']:.3e} ± {ng['pbdot_intrinsic_std']:.3e}"
    )
    print()
    row = summary["b1913_effective_row"]
    print("=== Effective B1913 row after branch selection ===")
    print(
        f"sbar = {row['sbar']:.6f}, xc = {row['xc']:.6f}, "
        f"sigma_delta = {row['sigma_delta']:.3e}, "
        f"delta_gamma_obs = {row['delta_gamma_obs']:.3e}"
    )
    print()
    b1913 = summary["b1913_only"]
    print("=== B1913-only posterior ===")
    print(
        f"|zeta1|_95 = {b1913['stats']['zeta1']['abs_95']:.3e}, "
        f"|zeta2|_95 = {b1913['stats']['zeta2']['abs_95']:.3e}"
    )
    print(
        f"eta(s_B1913={b1913['eta_at_b1913']['sbar']:.3f}) = "
        f"{b1913['eta_at_b1913']['mean']:.3e} ± {b1913['eta_at_b1913']['std']:.3e}, "
        f"|eta|_95 = {b1913['eta_at_b1913']['abs_95']:.3e}"
    )
    print()
    combined = summary["combined_with_request6"]
    print("=== Current Request 6 + B1913 ===")
    print(
        f"|zeta1|_95 = {combined['stats']['zeta1']['abs_95']:.3e}, "
        f"|zeta2|_95 = {combined['stats']['zeta2']['abs_95']:.3e}"
    )
    print(
        f"s_ref = {combined['sref']:.6f}, "
        f"eta_* = {combined['basis']['eta_star']['mean']:.3e} ± {combined['basis']['eta_star']['std']:.3e}, "
        f"|eta_*|_95 = {combined['basis']['eta_star']['abs_95']:.3e}"
    )
    print(
        f"kappa_* = {combined['basis']['kappa_star']['mean']:.3e} ± {combined['basis']['kappa_star']['std']:.3e}, "
        f"|kappa_*|_95 = {combined['basis']['kappa_star']['abs_95']:.3e}"
    )
    print()
    update = summary["lever_arm_update"]
    print("=== Lever-arm audit update at the original Request 6 s* ===")
    print(
        f"current |kappa_*|_95 = {update['current_only']['abs95_kappa_star']:.3e}, "
        f"with B1913 = {update['with_b1913']['abs95_kappa_star']:.3e}, "
        f"improvement = {update['improvement_factor']:.2f}x"
    )
    print()
    meta = combined["branch_selected_meta_at_mean"]
    print("=== B1913 branch-selected meta at combined posterior mean ===")
    print(
        f"mp = {meta['mp_mean_msun']:.6f} ± {meta['mp_std_msun']:.6f} Msun, "
        f"mc = {meta['mc_mean_msun']:.6f} ± {meta['mc_std_msun']:.6f} Msun, "
        f"gamma_obs/gamma_GR - 1 = {meta['gamma_obs_over_gr_minus_one']:.3e}"
    )


def main(output_dir: str | Path = ".") -> None:
    summary = build_summary(output_dir)
    print_summary(summary)
    print()
    print(f"wrote {Path(output_dir) / 'request6_b1913_covariance_summary.json'}")
    print(f"wrote {Path(output_dir) / 'request6_b1913_covariance_summary.svg'}")


if __name__ == "__main__":
    main()
