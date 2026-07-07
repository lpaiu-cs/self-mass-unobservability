from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent
REQ5_SUMMARY = REPO_ROOT / "request5" / "request5_j0337_phaseA_summary.json"
REQ6_SUMMARY = REPO_ROOT / "request6" / "request6_clock_sector_summary.json"
REQ6_LEVER_ARM_SUMMARY = REPO_ROOT / "request6" / "request6_lever_arm_audit_summary.json"
REQ3_SUMMARY = REPO_ROOT / "request3" / "request3_llr_mock_summary.json"

SUMMARY_JSON = ROOT / "request7_joint_consistency_summary.json"
SUMMARY_SVG = ROOT / "request7_joint_consistency_summary.svg"
SCENARIOS_TSV = ROOT / "request7_joint_consistency_scenarios.tsv"
REFERENCE_TIED_TSV = ROOT / "request7_joint_consistency_reference_tied_posterior.tsv"


@dataclass(frozen=True)
class PhaseABound:
    name: str
    delta_95: float
    hypothesis: str = ""
    source_label: str = ""
    source_url: str = ""


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


@dataclass(frozen=True)
class ClockPrior:
    name: str
    eta_abs_max: float
    kappa_abs_max: float
    grid_points: int


def repo_relative(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def sigma_grids(config: PhaseAConfig) -> tuple[np.ndarray, np.ndarray]:
    sigma1 = np.linspace(config.sigma1_min, config.sigma1_max, config.sigma1_points)
    sigma2 = np.linspace(config.sigma2_min, config.sigma2_max, config.sigma2_points)
    return sigma1, sigma2


def gaussian_sigma_from_95(delta_95: float) -> float:
    return delta_95 / 1.96


def normalize(weights: np.ndarray) -> np.ndarray:
    total = float(np.sum(weights))
    return weights / total


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    idx = int(np.searchsorted(cumulative, quantile, side="left"))
    return float(sorted_values[min(idx, sorted_values.size - 1)])


def weighted_abs95(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(np.abs(values))
    cumulative = np.cumsum(weights[order])
    idx = int(np.searchsorted(cumulative, 0.95, side="left"))
    return float(np.abs(values[order[min(idx, values.size - 1)]]))


def weighted_stats(values: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    mean = float(np.sum(weights * values))
    std = float(np.sqrt(np.sum(weights * (values - mean) ** 2)))
    mode = float(values[int(np.argmax(weights))])
    return {
        "mean": mean,
        "std": std,
        "mode": mode,
        "abs_95": weighted_abs95(values, weights),
        "q16": weighted_quantile(values, weights, 0.16),
        "q84": weighted_quantile(values, weights, 0.84),
        "q025": weighted_quantile(values, weights, 0.025),
        "q975": weighted_quantile(values, weights, 0.975),
    }


def summarize_sigma_posterior(sigma1: np.ndarray, sigma2: np.ndarray, posterior: np.ndarray) -> dict[str, dict[str, float]]:
    p1 = np.sum(posterior, axis=1)
    p2 = np.sum(posterior, axis=0)
    return {
        "sigma1": weighted_stats(sigma1, p1 / np.sum(p1)),
        "sigma2": weighted_stats(sigma2, p2 / np.sum(p2)),
    }


def phase_a_likelihood(config: PhaseAConfig, bound: PhaseABound, prior: EOSPrior) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sigma1, sigma2 = sigma_grids(config)
    s_ns = np.linspace(prior.s_ns_min, prior.s_ns_max, config.s_ns_points)
    sigma_delta = gaussian_sigma_from_95(bound.delta_95)
    delta = (
        sigma1[:, None, None] * (s_ns[None, None, :] - config.s_wd)
        + sigma2[None, :, None] * (s_ns[None, None, :] ** 2 - config.s_wd**2)
    )
    like = np.mean(np.exp(-0.5 * (delta / sigma_delta) ** 2), axis=2)
    return sigma1, sigma2, like


def load_clock_surrogate(summary_path: Path) -> dict[str, object]:
    summary = json.loads(summary_path.read_text())
    sstar = float(summary["effective_combinations"]["linearized_basis"]["clock_only"]["sstar"])
    z1 = np.array(summary["decoupled_grid"]["zeta1"], dtype=float)
    z2 = np.array(summary["decoupled_grid"]["zeta2"], dtype=float)
    posterior = normalize(np.array(summary["models"]["clock_only"]["posterior"], dtype=float))

    eta = z1[:, None] * sstar + z2[None, :] * sstar**2
    kappa = z1[:, None] + 2.0 * sstar * z2[None, :]

    w = posterior.reshape(-1)
    eta_flat = eta.reshape(-1)
    kappa_flat = kappa.reshape(-1)
    mean_eta = float(np.sum(w * eta_flat))
    mean_kappa = float(np.sum(w * kappa_flat))
    cov_eta_eta = float(np.sum(w * (eta_flat - mean_eta) ** 2))
    cov_kappa_kappa = float(np.sum(w * (kappa_flat - mean_kappa) ** 2))
    cov_eta_kappa = float(np.sum(w * (eta_flat - mean_eta) * (kappa_flat - mean_kappa)))
    cov = np.array([[cov_eta_eta, cov_eta_kappa], [cov_eta_kappa, cov_kappa_kappa]], dtype=float)
    cov += 1.0e-18 * np.eye(2)
    inv_cov = np.linalg.inv(cov)

    return {
        "summary": summary,
        "sstar": sstar,
        "mean": np.array([mean_eta, mean_kappa], dtype=float),
        "cov": cov,
        "inv_cov": inv_cov,
        "eta_stats": weighted_stats(eta_flat, w),
        "kappa_stats": weighted_stats(kappa_flat, w),
        "correlation": float(cov_eta_kappa / math.sqrt(cov_eta_eta * cov_kappa_kappa)),
    }


def build_clock_surrogate_variants(
    req6_summary_path: Path,
    lever_arm_path: Path,
) -> dict[str, dict[str, object]]:
    """Return the clock-branch Gaussian surrogates used in the joint scaffold.

    grid_matched: moments matched to the Request 6 clock-only grid posterior.
    Its kappa_* width inherits the flat zeta prior box of that grid
    (box-conditional), so it overstates the clock slope information.

    eta_only_kappa_free: eta_* moments from the same grid posterior (the
    data-driven band constraint), while the kappa_* width is taken from the
    Request 6 Fisher lever-arm audit of the same systems and the eta-kappa
    correlation is dropped. This removes the zeta-prior-box artifact and is
    the reference variant.
    """
    grid_matched = load_clock_surrogate(req6_summary_path)
    grid_matched["note"] = (
        "Moment-matched to the Request 6 clock-only grid posterior; kappa_* width is "
        "zeta-prior-box-conditional and overstates the clock slope information."
    )
    variants: dict[str, dict[str, object]] = {"grid_matched": grid_matched}
    if lever_arm_path.exists():
        lever = json.loads(lever_arm_path.read_text())
        kappa95 = float(lever["current_basis_gaussian"]["abs95_kappa_star"])
        eta_mean = float(grid_matched["mean"][0])
        eta_std = float(math.sqrt(grid_matched["cov"][0, 0]))
        kappa_std = kappa95 / 1.96
        cov = np.array([[eta_std**2, 0.0], [0.0, kappa_std**2]], dtype=float)
        variants["eta_only_kappa_free"] = {
            "sstar": grid_matched["sstar"],
            "mean": np.array([eta_mean, 0.0], dtype=float),
            "cov": cov,
            "inv_cov": np.linalg.inv(cov),
            "eta_stats": grid_matched["eta_stats"],
            "kappa_stats": {
                "mean": 0.0,
                "std": kappa_std,
                "abs_95": kappa95,
                "source": "request6_lever_arm_audit_summary.json:current_basis_gaussian.abs95_kappa_star",
            },
            "correlation": 0.0,
            "note": (
                "Data-driven variant: eta_* from the Request 6 grid posterior (band "
                f"constraint); kappa_* width from the Fisher lever-arm audit (|kappa_*|_95 = "
                f"{kappa95:.3e}); correlation dropped. Removes the zeta-prior-box artifact."
            ),
        }
    return variants


def local_clock_likelihood(
    eta: np.ndarray,
    kappa: np.ndarray,
    mean: np.ndarray,
    inv_cov: np.ndarray,
) -> np.ndarray:
    d0 = eta - mean[0]
    d1 = kappa - mean[1]
    mahal = inv_cov[0, 0] * d0**2 + 2.0 * inv_cov[0, 1] * d0 * d1 + inv_cov[1, 1] * d1**2
    return np.exp(-0.5 * mahal)


def clock_prior_stats(clock_prior: ClockPrior, surrogate: dict[str, object]) -> dict[str, object]:
    eta_grid = np.linspace(-clock_prior.eta_abs_max, clock_prior.eta_abs_max, clock_prior.grid_points)
    kappa_grid = np.linspace(-clock_prior.kappa_abs_max, clock_prior.kappa_abs_max, clock_prior.grid_points)
    eta_mesh = eta_grid[:, None]
    kappa_mesh = kappa_grid[None, :]
    like = local_clock_likelihood(eta_mesh, kappa_mesh, surrogate["mean"], surrogate["inv_cov"])
    evidence = float(np.mean(like))
    posterior = normalize(like)
    eta_weights = np.sum(posterior, axis=1)
    kappa_weights = np.sum(posterior, axis=0)
    return {
        "prior": asdict(clock_prior),
        "evidence": evidence,
        "stats": {
            "eta_star": weighted_stats(eta_grid, eta_weights / np.sum(eta_weights)),
            "kappa_star": weighted_stats(kappa_grid, kappa_weights / np.sum(kappa_weights)),
        },
    }


def transformed_basis_stats(
    eta_map: np.ndarray,
    kappa_map: np.ndarray,
    weights_2d: np.ndarray,
) -> dict[str, object]:
    weights = normalize(weights_2d).reshape(-1)
    eta_flat = eta_map.reshape(-1)
    kappa_flat = kappa_map.reshape(-1)
    mean_eta = float(np.sum(weights * eta_flat))
    mean_kappa = float(np.sum(weights * kappa_flat))
    cov_eta_eta = float(np.sum(weights * (eta_flat - mean_eta) ** 2))
    cov_kappa_kappa = float(np.sum(weights * (kappa_flat - mean_kappa) ** 2))
    cov_eta_kappa = float(np.sum(weights * (eta_flat - mean_eta) * (kappa_flat - mean_kappa)))
    return {
        "eta_star": weighted_stats(eta_flat, weights),
        "kappa_star": weighted_stats(kappa_flat, weights),
        "covariance": {
            "eta_eta": cov_eta_eta,
            "eta_kappa": cov_eta_kappa,
            "kappa_kappa": cov_kappa_kappa,
        },
    }


def mahalanobis_distance(point: np.ndarray, mean: np.ndarray, inv_cov: np.ndarray) -> float:
    diff = point - mean
    return float(math.sqrt(np.dot(diff, inv_cov @ diff)))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_posterior_table(path: Path, sigma1: np.ndarray, sigma2: np.ndarray, posterior: np.ndarray) -> None:
    rows = []
    for i, s1 in enumerate(sigma1):
        for j, s2 in enumerate(sigma2):
            rows.append([s1, s2, posterior[i, j]])
    np.savetxt(path, np.array(rows), header="sigma1\tsigma2\tposterior", comments="", delimiter="\t")


def write_summary_svg(path: Path, summary: dict[str, object]) -> None:
    width, height = 1560, 1500
    surrogates = summary["clock_surrogates"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="46" font-family="sans-serif" font-size="24" font-weight="700" fill="#111827">Request 7: provisional joint consistency scaffold</text>',
        '<text x="48" y="74" font-family="sans-serif" font-size="14" fill="#475569">Decoupled vs tied comparison using J0337 Phase A plus the Request 6 local clock audit. Request 3 is reference-only and not used in the likelihood.</text>',
        '<rect x="48" y="108" width="1464" height="190" fill="white" stroke="#cbd5e1"/>',
        f'<text x="66" y="138" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Clock Surrogate Variants (reference: {summary["reference_variant"]})</text>',
        f'<text x="66" y="166" font-family="monospace" font-size="13" fill="#111827">s* = {summary["clock_surrogate"]["sstar"]:.5f}</text>',
    ]
    y = 190
    for name, payload in surrogates.items():
        parts.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="12" fill="#111827">'
            f'{name}: eta_* = {payload["eta_star"]["mean"]:.3e} ± {payload["eta_star"]["std"]:.3e}, '
            f'kappa_* std = {payload["kappa_star"]["std"]:.3e}, corr = {payload["correlation"]:.3f}</text>'
        )
        y += 22
    parts.append('<text x="66" y="258" font-family="sans-serif" font-size="12" fill="#475569">grid_matched inherits the Request 6 zeta prior box in kappa_* (box-conditional); eta_only_kappa_free uses the Fisher lever-arm kappa_* width and is the data-driven reference.</text>')
    parts.append('<text x="66" y="278" font-family="sans-serif" font-size="12" fill="#475569">The clock branch is a local Gaussian surrogate in (eta_*, kappa_*) derived from Request 6, not a new TOA fit.</text>')
    parts.append('<text x="66" y="330" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Scenario Table (Bayes factor decoupled over tied; BF is an Occam-volume ratio set by the clock prior box)</text>')

    parts.append('<text x="66" y="358" font-family="monospace" font-size="11" fill="#374151">surrogate            bound        eos    clock_prior   log10 BF(dec/tied)   pre-clock(Mah)   tied |sigma1|_95   tied |sigma2|_95   tied |eta_*|_95   tied |kappa_*|_95</text>')
    y = 358
    for row in summary["scenario_rows"]:
        y += 18
        parts.append(
            f'<text x="66" y="{y}" font-family="monospace" font-size="11" fill="#111827">'
            f'{row["clock_surrogate"]:<21}{row["bound"]:<11}{row["eos"]:<7}{row["clock_prior"]:<14}'
            f'{row["log10_bf_dec_over_tied"]:>9.3f}{"":8}'
            f'{row["pre_clock_mahalanobis"]:>9.3f}{"":8}'
            f'{row["tied_sigma1_abs95"]:>9.3e}{"":6}'
            f'{row["tied_sigma2_abs95"]:>9.3e}{"":6}'
            f'{row["tied_eta_abs95"]:>9.3e}{"":5}'
            f'{row["tied_kappa_abs95"]:>9.3e}</text>'
        )

    ref = summary["reference_scenario"]
    box_y = y + 26
    parts.extend(
        [
            f'<rect x="48" y="{box_y}" width="1464" height="330" fill="white" stroke="#cbd5e1"/>',
            f'<text x="66" y="{box_y + 30}" font-family="sans-serif" font-size="16" font-weight="700" fill="#111827">Reference Scenario</text>',
            f'<text x="66" y="{box_y + 60}" font-family="monospace" font-size="13" fill="#111827">{ref["name"]}</text>',
            f'<text x="66" y="{box_y + 88}" font-family="monospace" font-size="13" fill="#111827">log10 BF(dec/tied) = {ref["comparison"]["log10_bf_dec_over_tied"]:.3f} [prior-volume-driven]</text>',
            f'<text x="66" y="{box_y + 112}" font-family="monospace" font-size="13" fill="#111827">pre-clock projected tension = {ref["comparison"]["pre_clock_mahalanobis"]:.3f}</text>',
            f'<text x="66" y="{box_y + 136}" font-family="monospace" font-size="13" fill="#111827">tied |sigma1|_95 = {ref["tied"]["sigma_stats"]["sigma1"]["abs_95"]:.3e}, tied |sigma2|_95 = {ref["tied"]["sigma_stats"]["sigma2"]["abs_95"]:.3e} [box-conditional, open ridge]</text>',
            f'<text x="66" y="{box_y + 160}" font-family="monospace" font-size="13" fill="#111827">tied |eta_*|_95 = {ref["tied"]["basis_stats"]["eta_star"]["abs_95"]:.3e} [data-limited], tied |kappa_*|_95 = {ref["tied"]["basis_stats"]["kappa_star"]["abs_95"]:.3e} [box-conditional]</text>',
            f'<text x="66" y="{box_y + 184}" font-family="monospace" font-size="13" fill="#111827">clock prior box = eta in [-{ref["clock_prior"]["eta_abs_max"]:.2e}, +{ref["clock_prior"]["eta_abs_max"]:.2e}], kappa in [-{ref["clock_prior"]["kappa_abs_max"]:.2e}, +{ref["clock_prior"]["kappa_abs_max"]:.2e}]</text>',
            f'<text x="66" y="{box_y + 212}" font-family="monospace" font-size="12" fill="#111827">Delta inputs: Voisin et al. (2025), A&amp;A 693, A143 (arXiv:2411.10066), planet vs red-noise hypotheses.</text>',
            f'<text x="66" y="{box_y + 240}" font-family="sans-serif" font-size="12" fill="#475569">{summary["request3_reference"]["note"]}</text>',
            '</svg>',
        ]
    )
    path.write_text("\n".join(parts) + "\n")


def surrogate_payload(surrogate: dict[str, object]) -> dict[str, object]:
    cov = np.asarray(surrogate["cov"], dtype=float)
    return {
        "sstar": float(surrogate["sstar"]),
        "eta_star": surrogate["eta_stats"],
        "kappa_star": surrogate["kappa_stats"],
        "covariance": {
            "eta_eta": float(cov[0, 0]),
            "eta_kappa": float(cov[0, 1]),
            "kappa_kappa": float(cov[1, 1]),
        },
        "correlation": float(surrogate["correlation"]),
        "note": surrogate.get("note", ""),
    }


def main() -> None:
    phase_config = PhaseAConfig()
    bounds = {
        "optimistic": PhaseABound(
            "optimistic",
            1.5e-6,
            hypothesis="circum-ternary planet timing-noise model",
            source_label="Voisin et al. (2025), A&A 693, A143, |Delta| < 1.5e-6 (95% CL, planet hypothesis)",
            source_url="https://arxiv.org/abs/2411.10066",
        ),
        "conservative": PhaseABound(
            "conservative",
            2.3e-6,
            hypothesis="achromatic red-noise timing-noise model",
            source_label="Voisin et al. (2025), A&A 693, A143, |Delta| < 2.3e-6 (95% CL, red-noise hypothesis)",
            source_url="https://arxiv.org/abs/2411.10066",
        ),
    }
    eos_priors = {
        "low": EOSPrior("low", 0.10, 0.15),
        "wide": EOSPrior("wide", 0.10, 0.20),
        "high": EOSPrior("high", 0.15, 0.20),
    }
    clock_priors = {
        "tight": ClockPrior("tight", 3.0e-3, 5.0e-2, 401),
        "medium": ClockPrior("medium", 1.0e-2, 1.0e-1, 401),
        "wide": ClockPrior("wide", 3.0e-2, 2.0e-1, 401),
    }

    surrogate_variants = build_clock_surrogate_variants(REQ6_SUMMARY, REQ6_LEVER_ARM_SUMMARY)
    reference_variant = "eta_only_kappa_free" if "eta_only_kappa_free" in surrogate_variants else "grid_matched"
    clock_prior_results = {
        variant_name: {name: clock_prior_stats(prior, surrogate) for name, prior in clock_priors.items()}
        for variant_name, surrogate in surrogate_variants.items()
    }

    sigma1, sigma2 = sigma_grids(phase_config)
    scenario_rows: list[dict[str, object]] = []
    scenarios: list[dict[str, object]] = []
    reference_name = f"optimistic / wide / medium / {reference_variant}"
    reference_scenario: dict[str, object] | None = None

    sstar = float(surrogate_variants["grid_matched"]["sstar"])
    eta_map = sigma1[:, None] * sstar + sigma2[None, :] * sstar**2
    kappa_map = sigma1[:, None] + 2.0 * sstar * sigma2[None, :]

    for variant_name, surrogate in surrogate_variants.items():
        clock_like_map = local_clock_likelihood(eta_map, kappa_map, surrogate["mean"], surrogate["inv_cov"])
        for bound_name, bound in bounds.items():
            for eos_name, eos_prior in eos_priors.items():
                _, _, jlike = phase_a_likelihood(phase_config, bound, eos_prior)
                jposterior = normalize(jlike)
                sigma_stats_dec = summarize_sigma_posterior(sigma1, sigma2, jposterior)
                projected_basis = transformed_basis_stats(eta_map, kappa_map, jposterior)
                projected_mean = np.array(
                    [
                        projected_basis["eta_star"]["mean"],
                        projected_basis["kappa_star"]["mean"],
                    ],
                    dtype=float,
                )
                pre_clock_mahal = mahalanobis_distance(projected_mean, surrogate["mean"], surrogate["inv_cov"])
                clock_like_expectation = float(np.sum(jposterior * clock_like_map))

                for clock_prior_name, clock_prior_result in clock_prior_results[variant_name].items():
                    tied_like = jlike * clock_like_map
                    tied_posterior = normalize(tied_like)
                    tied_sigma_stats = summarize_sigma_posterior(sigma1, sigma2, tied_posterior)
                    tied_basis_stats = transformed_basis_stats(eta_map, kappa_map, tied_posterior)
                    dec_evidence = float(np.mean(jlike) * clock_prior_result["evidence"])
                    tied_evidence = float(np.mean(tied_like))
                    bf_dec_over_tied = dec_evidence / tied_evidence
                    scenario = {
                        "name": f"{bound_name} / {eos_name} / {clock_prior_name} / {variant_name}",
                        "bound": bound_name,
                        "eos": eos_name,
                        "clock_surrogate": variant_name,
                        "clock_prior": clock_prior_result["prior"],
                        "decoupled": {
                            "sigma_stats": sigma_stats_dec,
                            "clock_stats": clock_prior_result["stats"],
                            "evidence": dec_evidence,
                        },
                        "tied": {
                            "sigma_stats": tied_sigma_stats,
                            "basis_stats": tied_basis_stats,
                            "evidence": tied_evidence,
                        },
                        "comparison": {
                            "bf_dec_over_tied": bf_dec_over_tied,
                            "log10_bf_dec_over_tied": float(math.log10(bf_dec_over_tied)),
                            "pre_clock_mahalanobis": pre_clock_mahal,
                            "clock_like_expectation_under_j0337": clock_like_expectation,
                        },
                    }
                    scenarios.append(scenario)
                    scenario_rows.append(
                        {
                            "clock_surrogate": variant_name,
                            "bound": bound_name,
                            "eos": eos_name,
                            "clock_prior": clock_prior_name,
                            "bf_dec_over_tied": bf_dec_over_tied,
                            "log10_bf_dec_over_tied": float(math.log10(bf_dec_over_tied)),
                            "pre_clock_mahalanobis": pre_clock_mahal,
                            "clock_like_expectation_under_j0337": clock_like_expectation,
                            "tied_sigma1_abs95": tied_sigma_stats["sigma1"]["abs_95"],
                            "tied_sigma2_abs95": tied_sigma_stats["sigma2"]["abs_95"],
                            "tied_eta_abs95": tied_basis_stats["eta_star"]["abs_95"],
                            "tied_kappa_abs95": tied_basis_stats["kappa_star"]["abs_95"],
                        }
                    )
                    if scenario["name"] == reference_name:
                        reference_scenario = scenario
                        write_posterior_table(REFERENCE_TIED_TSV, sigma1, sigma2, tied_posterior)

    if reference_scenario is None:
        raise RuntimeError("reference scenario was not created")

    write_tsv(SCENARIOS_TSV, scenario_rows)

    req3 = json.loads(REQ3_SUMMARY.read_text())
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "methodology": {
            "decoupled_model": "L_J0337(sigma1,sigma2) times L_clock(eta_*,kappa_*) with independent parameters",
            "tied_model": "zeta_i = sigma_i, mapped to eta_* and kappa_* at the Request 6 local expansion point",
            "request3_usage": "weak-field scale reference only; not included in the joint likelihood or evidence",
            "clock_surrogate_note": (
                "Clock likelihood is a Gaussian surrogate in (eta_*,kappa_*). Two variants are "
                "swept: grid_matched (moments from the Request 6 grid posterior; kappa_* width "
                "is zeta-prior-box-conditional) and eta_only_kappa_free (eta_* from the grid "
                "posterior, kappa_* width from the Fisher lever-arm audit; reference variant)."
            ),
            "bounds_source": {name: asdict(bound) for name, bound in bounds.items()},
            "bayes_factor_note": (
                "The clock surrogate likelihood is nearly constant over the J0337 sigma grid, "
                "so BF(dec/tied) reduces to an Occam-volume ratio set almost entirely by the "
                "assumed clock prior box; it is a prior-sensitivity audit, not a data-driven "
                "discriminant. The pre_clock_mahalanobis is likewise structurally stable "
                "because the symmetric J0337 posterior always projects to a near-zero mean."
            ),
            "tied_marginals_note": (
                "tied_sigma1_abs95/tied_sigma2_abs95/tied_kappa_abs95 inherit the open "
                "degeneracy ridge of the Phase A translation and its prior box "
                "(see request5 box_scale_check); only the eta_* direction is data-limited."
            ),
        },
        "inputs": {
            "request5_phaseA_summary": repo_relative(REQ5_SUMMARY),
            "request6_clock_sector_summary": repo_relative(REQ6_SUMMARY),
            "request6_lever_arm_audit_summary": repo_relative(REQ6_LEVER_ARM_SUMMARY),
            "request3_mock_summary": repo_relative(REQ3_SUMMARY),
        },
        "reference_variant": reference_variant,
        "clock_surrogate": surrogate_payload(surrogate_variants[reference_variant]),
        "clock_surrogates": {
            name: surrogate_payload(surrogate) for name, surrogate in surrogate_variants.items()
        },
        "clock_prior_sweeps": clock_prior_results,
        "scenario_rows": scenario_rows,
        "reference_scenario": reference_scenario,
        "request3_reference": {
            "analytic_cosD_amplitude_m_per_sigma1": req3["analytic_coefficients_m"]["sigma1"],
            "analytic_sigma2_coeff_m": req3["analytic_coefficients_m"]["sigma2"],
            "note": "Request 3 remains a weak-field scale sanity check only; these coefficients were not used as evidence in Request 7.",
        },
        "scenarios": scenarios,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    write_summary_svg(SUMMARY_SVG, summary)


if __name__ == "__main__":
    main()
