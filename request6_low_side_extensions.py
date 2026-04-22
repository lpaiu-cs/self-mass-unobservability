from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

import request6_b1913_covariance as b1913
import request6_clock_sector as r6
import request6_lever_arm_audit as audit


@dataclass(frozen=True)
class InclinationBranch:
    label: str
    inclination_mean_deg: float
    inclination_sigma_deg: float
    prior_weight: float
    source_url: str
    note: str
    correlated_with_gamma: bool = False


@dataclass(frozen=True)
class LowSideSource:
    name: str
    pb_days: float
    pb_sigma_days: float
    eccentricity: float
    eccentricity_sigma: float
    projected_semimajor_axis_s: float
    projected_semimajor_axis_sigma_s: float
    periastron_advance_deg_yr: float
    periastron_advance_sigma_deg_yr: float
    gamma_obs_s: float
    gamma_sigma_s: float
    system_type: str
    source_label: str
    source_url: str
    branches: tuple[InclinationBranch, ...]
    nuisance_notes: tuple[str, ...]


@dataclass(frozen=True)
class LowSideConfig:
    observable_model: str = "potential_only"
    mass_samples_per_branch: int = 4096
    batch_size: int = 128
    rng_seed: int = 654500
    baseline_summary_path: str = "request6_b1913_covariance_summary.json"
    baseline_posterior_path: str = "request6_b1913_covariance_posterior_combined.tsv"
    request6_summary_path: str = "request6_clock_sector_summary.json"


J1141 = LowSideSource(
    name="PSR J1141-6545",
    pb_days=0.1976509593,
    pb_sigma_days=1.0e-10,
    eccentricity=0.171884,
    eccentricity_sigma=2.0e-6,
    projected_semimajor_axis_s=1.858922,
    projected_semimajor_axis_sigma_s=6.0e-6,
    periastron_advance_deg_yr=5.3096,
    periastron_advance_sigma_deg_yr=4.0e-4,
    gamma_obs_s=7.73e-4,
    gamma_sigma_s=1.1e-5,
    system_type="NS-WD",
    source_label="Bhat, Bailes, Verbiest (2008) with xdot caveat from Venkatraman Krishnan et al. (2020)",
    source_url="https://arxiv.org/abs/0804.0956",
    branches=(
        InclinationBranch(
            label="scintillation_geometry",
            inclination_mean_deg=76.0,
            inclination_sigma_deg=2.5,
            prior_weight=1.0,
            source_url="https://arxiv.org/abs/0804.0956",
            note="Independent scintillation-based inclination prior used for non-gamma mass inference.",
            correlated_with_gamma=False,
        ),
    ),
    nuisance_notes=(
        "2020 timing reports xdot_obs = (1.7 ± 0.3)e-13 s s^-1.",
        "Profile-precession arguments limit aberration to < 21% of xdot_obs at 99% confidence, so >79% is attributed to WD-driven spin-orbit coupling.",
        "The same paper quotes P_WD < 900 s at 99% confidence, or < 200 s with the prograde-rotation prior.",
    ),
)


J1906 = LowSideSource(
    name="PSR J1906+0746",
    pb_days=0.16599304525,
    pb_sigma_days=4.0e-11,
    eccentricity=0.0852989,
    eccentricity_sigma=4.0e-7,
    projected_semimajor_axis_s=1.4199546,
    projected_semimajor_axis_sigma_s=8.0e-7,
    periastron_advance_deg_yr=7.5841,
    periastron_advance_sigma_deg_yr=2.0e-4,
    gamma_obs_s=4.59e-4,
    gamma_sigma_s=2.0e-6,
    system_type="compact companion (DNS or WD)",
    source_label="van Leeuwen et al. (2026)",
    source_url="https://arxiv.org/abs/2602.05947",
    branches=(
        InclinationBranch(
            label="external_i",
            inclination_mean_deg=45.0,
            inclination_sigma_deg=3.0,
            prior_weight=0.5,
            source_url="https://arxiv.org/abs/2602.05947",
            note="External inclination from Desvignes et al. (2019) as quoted in the 2026 timing paper.",
            correlated_with_gamma=False,
        ),
        InclinationBranch(
            label="xdot_conditioned_i",
            inclination_mean_deg=50.0,
            inclination_sigma_deg=1.0,
            prior_weight=0.5,
            source_url="https://arxiv.org/abs/2602.05947",
            note="DDGR+xdot branch from the 2026 timing paper; kept as an explicit correlated nuisance branch.",
            correlated_with_gamma=True,
        ),
    ),
    nuisance_notes=(
        "2026 timing reports xdot_obs = (-1.8 ± 0.6)e-13 s s^-1.",
        "The same paper states that fitting xdot shifts the individual masses by ~3.5 sigma because xdot and gamma are correlated.",
        "In the xdot-conditioned branch the pulsar mass moves to ~1.42 Msun, pushing the source from the low-side toward the positive-s side.",
    ),
)


def load_posterior_grid(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = np.loadtxt(path, skiprows=1)
    z1 = np.unique(table[:, 0])
    z2 = np.unique(table[:, 1])
    posterior = table[:, 2].reshape(z1.size, z2.size)
    return z1, z2, posterior


def sample_branch(source: LowSideSource, branch: InclinationBranch, config: LowSideConfig, rng: np.random.Generator) -> dict[str, np.ndarray]:
    remaining = config.mass_samples_per_branch
    collected: dict[str, list[np.ndarray]] = {
        "pb_days": [],
        "eccentricity": [],
        "projected_semimajor_axis_s": [],
        "total_mass_msun": [],
        "inclination_deg": [],
        "mp_msun": [],
        "mc_msun": [],
        "self_gravity_fraction": [],
    }
    while remaining > 0:
        draw_count = max(1024, 2 * remaining)
        pb = rng.normal(source.pb_days, source.pb_sigma_days, draw_count)
        eccentricity = rng.normal(source.eccentricity, source.eccentricity_sigma, draw_count)
        x = rng.normal(source.projected_semimajor_axis_s, source.projected_semimajor_axis_sigma_s, draw_count)
        omdot = rng.normal(source.periastron_advance_deg_yr, source.periastron_advance_sigma_deg_yr, draw_count)
        inclination_deg = rng.normal(branch.inclination_mean_deg, branch.inclination_sigma_deg, draw_count)

        total_mass = r6.total_mass_from_omdot(pb, eccentricity, omdot)
        sin_i = np.sin(np.deg2rad(inclination_deg))
        mf = r6.mass_function(x, pb)
        companion_mass = np.cbrt(mf * total_mass**2 / np.clip(sin_i, 1.0e-12, None) ** 3)
        pulsar_mass = total_mass - companion_mass

        valid = (
            (pb > 0.0)
            & (eccentricity > 0.0)
            & (eccentricity < 1.0)
            & (inclination_deg > 0.0)
            & (inclination_deg < 90.0)
            & np.isfinite(total_mass)
            & np.isfinite(companion_mass)
            & (pulsar_mass > 0.0)
            & (companion_mass > 0.0)
        )
        if not np.any(valid):
            continue
        take = min(remaining, int(np.sum(valid)))
        idx = np.flatnonzero(valid)[:take]
        pb = pb[idx]
        eccentricity = eccentricity[idx]
        x = x[idx]
        total_mass = total_mass[idx]
        inclination_deg = inclination_deg[idx]
        pulsar_mass = pulsar_mass[idx]
        companion_mass = companion_mass[idx]
        sg_lo, sg_hi = r6.mass_conditioned_s_prior(pulsar_mass)
        sg_hi = np.maximum(sg_hi, sg_lo)
        self_gravity = rng.uniform(sg_lo, sg_hi)

        collected["pb_days"].append(pb)
        collected["eccentricity"].append(eccentricity)
        collected["projected_semimajor_axis_s"].append(x)
        collected["total_mass_msun"].append(total_mass)
        collected["inclination_deg"].append(inclination_deg)
        collected["mp_msun"].append(pulsar_mass)
        collected["mc_msun"].append(companion_mass)
        collected["self_gravity_fraction"].append(self_gravity)
        remaining -= take
    return {key: np.concatenate(value) for key, value in collected.items()}


def sample_source(source: LowSideSource, config: LowSideConfig) -> dict[str, dict[str, np.ndarray]]:
    rng = np.random.default_rng(config.rng_seed + abs(hash(source.name)) % 10_000)
    return {branch.label: sample_branch(source, branch, config, rng) for branch in source.branches}


def source_likelihood_grid(
    source: LowSideSource,
    branch_samples: dict[str, dict[str, np.ndarray]],
    z1: np.ndarray,
    z2: np.ndarray,
    config: LowSideConfig,
) -> np.ndarray:
    likelihood = np.zeros((z1.size, z2.size), dtype=float)
    for branch in source.branches:
        samples = branch_samples[branch.label]
        branch_like = np.zeros_like(likelihood)
        for start in range(0, samples["pb_days"].size, config.batch_size):
            stop = min(start + config.batch_size, samples["pb_days"].size)
            model = r6.gamma_clock(
                samples["pb_days"][start:stop],
                samples["eccentricity"][start:stop],
                samples["mp_msun"][start:stop],
                samples["mc_msun"][start:stop],
                samples["self_gravity_fraction"][start:stop],
                z1,
                z2,
                config.observable_model,
            )
            resid = source.gamma_obs_s - model
            branch_like += np.exp(-0.5 * (resid / source.gamma_sigma_s) ** 2).sum(axis=2)
        likelihood += branch.prior_weight * (branch_like / samples["pb_days"].size)
    return likelihood


def source_branch_selected_meta(
    source: LowSideSource,
    branch_samples: dict[str, dict[str, np.ndarray]],
    zeta1: float,
    zeta2: float,
    observable_model: str,
) -> dict[str, object]:
    weights_list = []
    branch_probabilities: dict[str, float] = {}
    arrays: dict[str, list[np.ndarray]] = {
        "mp_msun": [],
        "mc_msun": [],
        "self_gravity": [],
        "gamma_gr_s": [],
        "gamma_model_s": [],
        "xc": [],
        "inclination_deg": [],
        "total_mass_msun": [],
    }
    total_like = 0.0
    branch_like_sums: dict[str, float] = {}
    for branch in source.branches:
        samples = branch_samples[branch.label]
        gamma_model = r6.gamma_clock(
            samples["pb_days"],
            samples["eccentricity"],
            samples["mp_msun"],
            samples["mc_msun"],
            samples["self_gravity_fraction"],
            np.array([zeta1]),
            np.array([zeta2]),
            observable_model,
        )[0, 0]
        like = branch.prior_weight * np.exp(-0.5 * ((source.gamma_obs_s - gamma_model) / source.gamma_sigma_s) ** 2)
        weights_list.append(like)
        total_like += float(np.sum(like))
        branch_like_sums[branch.label] = float(np.sum(like))

        gamma_gr = r6.gamma_gr(samples["pb_days"], samples["eccentricity"], samples["mp_msun"], samples["mc_msun"])
        arrays["mp_msun"].append(samples["mp_msun"])
        arrays["mc_msun"].append(samples["mc_msun"])
        arrays["self_gravity"].append(samples["self_gravity_fraction"])
        arrays["gamma_gr_s"].append(gamma_gr)
        arrays["gamma_model_s"].append(gamma_model)
        arrays["xc"].append(samples["mc_msun"] / (samples["mp_msun"] + samples["mc_msun"]))
        arrays["inclination_deg"].append(samples["inclination_deg"])
        arrays["total_mass_msun"].append(samples["total_mass_msun"])

    for branch in source.branches:
        branch_probabilities[branch.label] = branch_like_sums[branch.label] / total_like if total_like > 0.0 else 0.0

    combined_weights = np.concatenate(weights_list)
    combined_weights /= np.sum(combined_weights)
    merged = {key: np.concatenate(value) for key, value in arrays.items()}

    def mean_std(values: np.ndarray) -> tuple[float, float]:
        mean = float(np.sum(combined_weights * values))
        std = float(np.sqrt(np.sum(combined_weights * (values - mean) ** 2)))
        return mean, std

    mp_mean, mp_std = mean_std(merged["mp_msun"])
    mc_mean, mc_std = mean_std(merged["mc_msun"])
    sg_mean, sg_std = mean_std(merged["self_gravity"])
    gamma_gr_mean, gamma_gr_std = mean_std(merged["gamma_gr_s"])
    gamma_model_mean, gamma_model_std = mean_std(merged["gamma_model_s"])
    xc_mean, xc_std = mean_std(merged["xc"])
    inc_mean, inc_std = mean_std(merged["inclination_deg"])
    mtot_mean, mtot_std = mean_std(merged["total_mass_msun"])

    return {
        "branch_probabilities": branch_probabilities,
        "mp_mean_msun": mp_mean,
        "mp_std_msun": mp_std,
        "mc_mean_msun": mc_mean,
        "mc_std_msun": mc_std,
        "self_gravity_mean": sg_mean,
        "self_gravity_std": sg_std,
        "gamma_pred_gr_s": gamma_gr_mean,
        "gamma_pred_gr_std_s": gamma_gr_std,
        "gamma_model_s": gamma_model_mean,
        "gamma_model_std_s": gamma_model_std,
        "xc_mean": xc_mean,
        "xc_std": xc_std,
        "inclination_mean_deg": inc_mean,
        "inclination_std_deg": inc_std,
        "total_mass_mean_msun": mtot_mean,
        "total_mass_std_msun": mtot_std,
        "gamma_obs_over_gr_minus_one": float(source.gamma_obs_s / gamma_gr_mean - 1.0),
    }


def make_effective_row(source: LowSideSource, meta: dict[str, object]) -> dict[str, float | str]:
    sigma_total = math.sqrt(source.gamma_sigma_s**2 + float(meta["gamma_model_std_s"]) ** 2)
    return {
        "name": source.name,
        "sbar": float(meta["self_gravity_mean"]),
        "xc": float(meta["xc_mean"]),
        "gamma_gr_s": float(meta["gamma_pred_gr_s"]),
        "sigma_gamma_s": sigma_total,
        "sigma_delta": float(sigma_total / float(meta["gamma_pred_gr_s"])),
        "delta_gamma_obs": float(source.gamma_obs_s / float(meta["gamma_pred_gr_s"]) - 1.0),
    }


def precision_weighted_sref(rows: list[dict[str, float | str]]) -> float:
    weights = np.array([1.0 / float(row["sigma_delta"]) ** 2 for row in rows], dtype=float)
    weights /= np.sum(weights)
    return float(np.sum(weights * np.array([float(row["sbar"]) for row in rows], dtype=float)))


def summarize_source(
    source: LowSideSource,
    branch_samples: dict[str, dict[str, np.ndarray]],
    z1: np.ndarray,
    z2: np.ndarray,
    baseline_posterior: np.ndarray,
    config: LowSideConfig,
) -> dict[str, object]:
    like = source_likelihood_grid(source, branch_samples, z1, z2, config)
    source_only_posterior = r6.posterior_from_prior(np.ones_like(like) / like.size, like)
    combined_posterior = r6.normalize(baseline_posterior * like)

    source_stats = r6.summarize_posterior(z1, z2, source_only_posterior)
    combined_stats = r6.summarize_posterior(z1, z2, combined_posterior)
    meta_gr = source_branch_selected_meta(source, branch_samples, 0.0, 0.0, config.observable_model)
    meta_mean = source_branch_selected_meta(
        source,
        branch_samples,
        combined_stats["zeta1"]["mean"],
        combined_stats["zeta2"]["mean"],
        config.observable_model,
    )

    return {
        "likelihood": like,
        "source_only_posterior": source_only_posterior,
        "combined_posterior": combined_posterior,
        "source_only_stats": source_stats,
        "combined_stats": combined_stats,
        "meta_at_gr": meta_gr,
        "meta_at_combined_mean": meta_mean,
        "effective_row": make_effective_row(source, meta_gr),
    }


def write_summary_svg(path: str | Path, summary: dict[str, object]) -> Path:
    width, height = 1400, 1020
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 6: Low-side source extensions</text>',
        '<text x="48" y="72" font-family="sans-serif" font-size="14" fill="#475569">Baseline is current Request 6 plus the B1913 covariance-aware extension. J1141 and J1906 are added as source-specific surrogate likelihoods with their nuisance structures left explicit.</text>',
    ]
    z1 = np.array(summary["zeta_grid"]["zeta1"])
    z2 = np.array(summary["zeta_grid"]["zeta2"])
    r6._svg_heatmap(parts, 70, 120, 380, 280, z1, z2, np.array(summary["baseline"]["posterior"]), "Baseline: current + B1913")
    r6._svg_heatmap(parts, 510, 120, 380, 280, z1, z2, np.array(summary["j1141"]["combined"]["posterior"]), "Baseline + J1141")
    r6._svg_heatmap(parts, 950, 120, 380, 280, z1, z2, np.array(summary["j1906"]["combined"]["posterior"]), "Baseline + J1906")
    parts.append('<rect x="70" y="450" width="1260" height="470" fill="white" stroke="#111827" stroke-width="1.2"/>')
    lines = [
        f"baseline |kappa_*|_95 = {summary['baseline']['lever_update']['abs95_kappa_star']:.3e}",
        f"+ J1141  -> |kappa_*|_95 = {summary['j1141']['lever_update']['abs95_kappa_star']:.3e}",
        f"+ J1906  -> |kappa_*|_95 = {summary['j1906']['lever_update']['abs95_kappa_star']:.3e}",
        f"+ both    -> |kappa_*|_95 = {summary['both']['lever_update']['abs95_kappa_star']:.3e}",
        f"J1141 effective row: sbar = {summary['j1141']['effective_row']['sbar']:.3f}, sigma_delta = {summary['j1141']['effective_row']['sigma_delta']:.3e}",
        f"J1906 effective row: sbar = {summary['j1906']['effective_row']['sbar']:.3f}, sigma_delta = {summary['j1906']['effective_row']['sigma_delta']:.3e}",
        f"J1906 branch weights at GR: ext_i = {summary['j1906']['meta_at_gr']['branch_probabilities']['external_i']:.3f}, xdot_i = {summary['j1906']['meta_at_gr']['branch_probabilities']['xdot_conditioned_i']:.3f}",
        "Interpretation: J1141 stays a genuine low-side source. J1906 becomes a nuisance-sensitive mixed source because its xdot-conditioned branch moves toward positive s.",
    ]
    parts.append('<text x="88" y="482" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Audit update</text>')
    for idx, line in enumerate(lines):
        parts.append(f'<text x="88" y="{512 + 28*idx}" font-family="monospace" font-size="14" fill="#111827">{line}</text>')
    parts.append("</svg>")
    output = Path(path)
    output.write_text("\n".join(parts))
    return output


def build_summary(output_dir: str | Path = ".") -> dict[str, object]:
    output_root = Path(output_dir)
    config = LowSideConfig()
    z1, z2, baseline_posterior = load_posterior_grid(output_root / config.baseline_posterior_path)
    baseline_summary = json.loads((output_root / config.baseline_summary_path).read_text())
    request6_summary = json.loads((output_root / config.request6_summary_path).read_text())
    current_sstar = float(request6_summary["effective_combinations"]["linearized_basis"]["clock_only"]["sstar"])

    baseline_rows = audit.current_rows(request6_summary) + [baseline_summary["b1913_effective_row"]]
    baseline_sref = precision_weighted_sref(baseline_rows)
    baseline_basis = r6.basis_stats_at_sstar(z1, z2, baseline_posterior, baseline_sref)
    baseline_update = audit.fisher_eta_kappa(baseline_rows, current_sstar)

    j1141_samples = sample_source(J1141, config)
    j1141 = summarize_source(J1141, j1141_samples, z1, z2, baseline_posterior, config)
    j1906_samples = sample_source(J1906, config)
    j1906 = summarize_source(J1906, j1906_samples, z1, z2, baseline_posterior, config)

    posterior_both = r6.normalize(baseline_posterior * j1141["likelihood"] * j1906["likelihood"])
    stats_both = r6.summarize_posterior(z1, z2, posterior_both)
    rows_j1141 = baseline_rows + [j1141["effective_row"]]
    rows_j1906 = baseline_rows + [j1906["effective_row"]]
    rows_both = baseline_rows + [j1141["effective_row"], j1906["effective_row"]]
    sref_both = precision_weighted_sref(rows_both)
    basis_both = r6.basis_stats_at_sstar(z1, z2, posterior_both, sref_both)
    update_j1141 = audit.fisher_eta_kappa(rows_j1141, current_sstar)
    update_j1906 = audit.fisher_eta_kappa(rows_j1906, current_sstar)
    update_both = audit.fisher_eta_kappa(rows_both, current_sstar)

    r6.write_posterior_table(output_root / "request6_low_side_posterior_j1141.tsv", z1, z2, j1141["combined_posterior"])
    r6.write_posterior_table(output_root / "request6_low_side_posterior_j1906.tsv", z1, z2, j1906["combined_posterior"])
    r6.write_posterior_table(output_root / "request6_low_side_posterior_both.tsv", z1, z2, posterior_both)

    summary = {
        "config": asdict(config),
        "zeta_grid": {"zeta1": z1.tolist(), "zeta2": z2.tolist()},
        "current_request6_sstar": current_sstar,
        "baseline": {
            "posterior": baseline_posterior.tolist(),
            "sref": baseline_sref,
            "basis": baseline_basis,
            "lever_update": baseline_update,
        },
        "j1141": {
            "system": asdict(J1141),
            "effective_row": j1141["effective_row"],
            "meta_at_gr": j1141["meta_at_gr"],
            "meta_at_combined_mean": j1141["meta_at_combined_mean"],
            "source_only": {
                "stats": j1141["source_only_stats"],
                "posterior": j1141["source_only_posterior"].tolist(),
            },
            "combined": {
                "stats": j1141["combined_stats"],
                "posterior": j1141["combined_posterior"].tolist(),
            },
            "lever_update": update_j1141,
        },
        "j1906": {
            "system": asdict(J1906),
            "effective_row": j1906["effective_row"],
            "meta_at_gr": j1906["meta_at_gr"],
            "meta_at_combined_mean": j1906["meta_at_combined_mean"],
            "source_only": {
                "stats": j1906["source_only_stats"],
                "posterior": j1906["source_only_posterior"].tolist(),
            },
            "combined": {
                "stats": j1906["combined_stats"],
                "posterior": j1906["combined_posterior"].tolist(),
            },
            "lever_update": update_j1906,
        },
        "both": {
            "stats": stats_both,
            "posterior": posterior_both.tolist(),
            "sref": sref_both,
            "basis": basis_both,
            "lever_update": update_both,
        },
    }
    (output_root / "request6_low_side_extensions_summary.json").write_text(json.dumps(summary, indent=2))
    write_summary_svg(output_root / "request6_low_side_extensions_summary.svg", summary)
    return summary


def print_summary(summary: dict[str, object]) -> None:
    print("=== Request 6 low-side extensions ===")
    print(
        f"baseline |kappa_*|_95 = {summary['baseline']['lever_update']['abs95_kappa_star']:.3e} "
        f"at original s* = {summary['current_request6_sstar']:.6f}"
    )
    print()
    for key in ["j1141", "j1906"]:
        item = summary[key]
        row = item["effective_row"]
        print(f"=== {item['system']['name']} ===")
        print(
            f"sbar = {row['sbar']:.6f}, xc = {row['xc']:.6f}, "
            f"sigma_delta = {row['sigma_delta']:.3e}, delta_gamma_obs = {row['delta_gamma_obs']:.3e}"
        )
        print(
            f"meta at GR: mp = {item['meta_at_gr']['mp_mean_msun']:.6f} ± {item['meta_at_gr']['mp_std_msun']:.6f}, "
            f"mc = {item['meta_at_gr']['mc_mean_msun']:.6f} ± {item['meta_at_gr']['mc_std_msun']:.6f}, "
            f"i = {item['meta_at_gr']['inclination_mean_deg']:.3f} ± {item['meta_at_gr']['inclination_std_deg']:.3f} deg"
        )
        if key == "j1906":
            probs = item["meta_at_gr"]["branch_probabilities"]
            print(
                f"branch weights at GR: external_i = {probs['external_i']:.6f}, "
                f"xdot_conditioned_i = {probs['xdot_conditioned_i']:.6f}"
            )
        print(
            f"with source |kappa_*|_95 = {item['lever_update']['abs95_kappa_star']:.3e}, "
            f"improvement vs baseline = {summary['baseline']['lever_update']['abs95_kappa_star'] / item['lever_update']['abs95_kappa_star']:.2f}x"
        )
        print()
    print("=== Baseline + both low-side sources ===")
    print(
        f"|zeta1|_95 = {summary['both']['stats']['zeta1']['abs_95']:.3e}, "
        f"|zeta2|_95 = {summary['both']['stats']['zeta2']['abs_95']:.3e}"
    )
    print(
        f"s_ref = {summary['both']['sref']:.6f}, "
        f"|eta_*|_95 = {summary['both']['basis']['eta_star']['abs_95']:.3e}, "
        f"|kappa_*|_95 = {summary['both']['basis']['kappa_star']['abs_95']:.3e}"
    )
    print(
        f"audit update at original s*: {summary['both']['lever_update']['abs95_kappa_star']:.3e}"
    )


def main(output_dir: str | Path = ".") -> None:
    summary = build_summary(output_dir)
    print_summary(summary)
    print()
    print(f"wrote {Path(output_dir) / 'request6_low_side_extensions_summary.json'}")
    print(f"wrote {Path(output_dir) / 'request6_low_side_extensions_summary.svg'}")


if __name__ == "__main__":
    main()
