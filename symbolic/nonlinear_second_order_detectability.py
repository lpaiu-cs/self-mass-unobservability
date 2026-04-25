from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite, sqrt

from amplitude_weighted_resonant_design import (
    forcing_cosine_coefficient,
    projection_magnitude,
)
from physical_detectability_map import harmonic_noise_sigma
from resonant_comparator_audit import ResonanceBudget, line_shape_budget


@dataclass(frozen=True)
class GeneratedLine:
    k: int
    forcing_squared_coefficient: float
    transfer_magnitude: float
    projection_magnitude: float
    observable_weight: float
    relative_weight: float
    observable_amplitude: float
    noise_sigma: float | None
    snr: float | None
    usable: bool


@dataclass(frozen=True)
class NonlinearDetectabilityDesign:
    label: str
    p: float
    eccentricity: float
    rho: Fraction
    damping: float
    channel: str
    generated_cutoff: int
    criterion: str
    relative_cutoff: float | None
    amplitude_scale: float | None
    base_noise_sigma: float | None
    snr_threshold: float | None
    sideband_order: int
    projection_nuisance: int
    usable_count: int
    lower_usable_count: int
    upper_usable_count: int
    budget: ResonanceBudget
    verdict: str


def _check_positive(name: str, value: float) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def _check_nonnegative(name: str, value: float) -> None:
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")


def generated_forcing_coefficient(
    p: float,
    eccentricity: float,
    harmonic: int,
    samples: int = 8192,
) -> float:
    return forcing_cosine_coefficient(
        2 * p,
        eccentricity,
        harmonic,
        samples=samples,
    )


def nonlinear_drive_transfer_magnitude(
    harmonic: int,
    rho: Fraction,
    damping: float,
) -> float:
    _check_positive("harmonic", harmonic)
    _check_positive("rho", float(rho))
    _check_nonnegative("damping", damping)
    rho_sq = float(rho) ** 2
    denominator = sqrt((rho_sq - harmonic**2) ** 2 + (damping * harmonic) ** 2)
    if denominator == 0:
        return float("inf")
    return 1 / denominator


def generated_line_rows(
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    generated_cutoff: int,
    criterion: str,
    relative_cutoff: float | None = None,
    amplitude_scale: float | None = None,
    base_noise_sigma: float | None = None,
    snr_threshold: float | None = None,
    harmonic_noise_slope: float = 0.0,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> tuple[GeneratedLine, ...]:
    _check_positive("generated_cutoff", generated_cutoff)
    if criterion not in {"relative", "snr"}:
        raise ValueError("criterion must be 'relative' or 'snr'")
    if criterion == "relative":
        if relative_cutoff is None:
            raise ValueError("relative_cutoff is required for relative criterion")
        _check_nonnegative("relative_cutoff", relative_cutoff)
    if criterion == "snr":
        if amplitude_scale is None or base_noise_sigma is None or snr_threshold is None:
            raise ValueError("amplitude_scale, base_noise_sigma, and snr_threshold are required")
        _check_positive("amplitude_scale", amplitude_scale)
        _check_positive("base_noise_sigma", base_noise_sigma)
        _check_positive("snr_threshold", snr_threshold)

    raw_rows: list[tuple[int, float, float, float, float]] = []
    for harmonic in range(1, generated_cutoff + 1):
        coefficient = generated_forcing_coefficient(
            p,
            eccentricity,
            harmonic,
            samples=samples,
        )
        transfer = nonlinear_drive_transfer_magnitude(harmonic, rho, damping)
        projection = projection_magnitude(channel, harmonic, kappa_ratio=kappa_ratio)
        weight = abs(coefficient) * transfer * projection
        raw_rows.append((harmonic, coefficient, transfer, projection, weight))
    finite_weights = [row[-1] for row in raw_rows if isfinite(row[-1])]
    scale = max(finite_weights, default=0.0)
    rows: list[GeneratedLine] = []
    for harmonic, coefficient, transfer, projection, weight in raw_rows:
        relative = 0.0 if scale == 0 else weight / scale
        observable = weight if amplitude_scale is None else amplitude_scale * weight
        sigma: float | None = None
        snr: float | None = None
        if criterion == "relative":
            usable = isfinite(weight) and weight > 0 and relative >= float(relative_cutoff)
        else:
            assert amplitude_scale is not None
            assert base_noise_sigma is not None
            assert snr_threshold is not None
            sigma = harmonic_noise_sigma(
                base_noise_sigma,
                harmonic,
                harmonic_noise_slope=harmonic_noise_slope,
            )
            snr = observable / sigma
            usable = snr >= snr_threshold
        rows.append(
            GeneratedLine(
                k=harmonic,
                forcing_squared_coefficient=coefficient,
                transfer_magnitude=transfer,
                projection_magnitude=projection,
                observable_weight=weight,
                relative_weight=relative,
                observable_amplitude=observable,
                noise_sigma=sigma,
                snr=snr,
                usable=usable,
            )
        )
    return tuple(rows)


def classify_nonlinear_detectability(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    generated_cutoff: int,
    criterion: str,
    sideband_order: int,
    projection_nuisance: int,
    relative_cutoff: float | None = None,
    amplitude_scale: float | None = None,
    base_noise_sigma: float | None = None,
    snr_threshold: float | None = None,
    harmonic_noise_slope: float = 0.0,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> NonlinearDetectabilityDesign:
    budget = line_shape_budget(sideband_order, projection_nuisance)
    rows = generated_line_rows(
        p,
        eccentricity,
        rho,
        damping,
        channel,
        generated_cutoff,
        criterion,
        relative_cutoff=relative_cutoff,
        amplitude_scale=amplitude_scale,
        base_noise_sigma=base_noise_sigma,
        snr_threshold=snr_threshold,
        harmonic_noise_slope=harmonic_noise_slope,
        kappa_ratio=kappa_ratio,
        samples=samples,
    )
    usable = tuple(row for row in rows if row.usable)
    lower = tuple(row for row in usable if row.k < rho)
    upper = tuple(row for row in usable if row.k > rho)
    if len(usable) < budget.minimum_frequency_samples:
        verdict = "generated-underbudget-no-go"
    elif not lower or not upper:
        verdict = "no-generated-resonance-bracket"
    else:
        verdict = "nonlinear-generated-budget-breaking"
    return NonlinearDetectabilityDesign(
        label=label,
        p=p,
        eccentricity=eccentricity,
        rho=rho,
        damping=damping,
        channel=channel,
        generated_cutoff=generated_cutoff,
        criterion=criterion,
        relative_cutoff=relative_cutoff,
        amplitude_scale=amplitude_scale,
        base_noise_sigma=base_noise_sigma,
        snr_threshold=snr_threshold,
        sideband_order=sideband_order,
        projection_nuisance=projection_nuisance,
        usable_count=len(usable),
        lower_usable_count=len(lower),
        upper_usable_count=len(upper),
        budget=budget,
        verdict=verdict,
    )


def minimum_generated_amplitude_scale(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    generated_cutoff: int,
    base_noise_sigma: float,
    snr_threshold: float,
    sideband_order: int,
    projection_nuisance: int,
    harmonic_noise_slope: float = 0.0,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> NonlinearDetectabilityDesign | None:
    unit_rows = generated_line_rows(
        p,
        eccentricity,
        rho,
        damping,
        channel,
        generated_cutoff,
        criterion="snr",
        amplitude_scale=1.0,
        base_noise_sigma=base_noise_sigma,
        snr_threshold=snr_threshold,
        harmonic_noise_slope=harmonic_noise_slope,
        kappa_ratio=kappa_ratio,
        samples=samples,
    )
    candidate_scales = sorted(
        snr_threshold / row.snr
        for row in unit_rows
        if row.snr is not None and row.snr > 0
    )
    for scale in candidate_scales:
        test_scale = scale * (1.0 + 1.0e-12)
        design = classify_nonlinear_detectability(
            label,
            p,
            eccentricity,
            rho,
            damping,
            channel,
            generated_cutoff,
            criterion="snr",
            sideband_order=sideband_order,
            projection_nuisance=projection_nuisance,
            amplitude_scale=test_scale,
            base_noise_sigma=base_noise_sigma,
            snr_threshold=snr_threshold,
            harmonic_noise_slope=harmonic_noise_slope,
            kappa_ratio=kappa_ratio,
            samples=samples,
        )
        if design.verdict == "nonlinear-generated-budget-breaking":
            return design
    return None


def default_relative_rows() -> tuple[NonlinearDetectabilityDesign, ...]:
    rows: list[NonlinearDetectabilityDesign] = []
    for eccentricity in (0.1, 0.3, 0.6):
        for channel, nuisance in (("acceleration", 1), ("range", 2)):
            rows.append(
                classify_nonlinear_detectability(
                    f"relative p=2 e={eccentricity} rho=3/2 {channel}",
                    p=2.0,
                    eccentricity=eccentricity,
                    rho=Fraction(3, 2),
                    damping=0.2,
                    channel=channel,
                    generated_cutoff=6,
                    criterion="relative",
                    relative_cutoff=1.0e-3,
                    sideband_order=1,
                    projection_nuisance=nuisance,
                    samples=4096,
                )
            )
    return tuple(rows)


def default_snr_rows() -> tuple[NonlinearDetectabilityDesign, ...]:
    rows: list[NonlinearDetectabilityDesign] = []
    for amplitude_scale in (0.1, 0.3, 1.0):
        rows.append(
            classify_nonlinear_detectability(
                f"snr acceleration scale={amplitude_scale}",
                p=2.0,
                eccentricity=0.3,
                rho=Fraction(3, 2),
                damping=0.2,
                channel="acceleration",
                generated_cutoff=6,
                criterion="snr",
                amplitude_scale=amplitude_scale,
                base_noise_sigma=1.0e-3,
                snr_threshold=5.0,
                sideband_order=1,
                projection_nuisance=1,
                samples=4096,
            )
        )
        rows.append(
            classify_nonlinear_detectability(
                f"snr range scale={amplitude_scale}",
                p=2.0,
                eccentricity=0.3,
                rho=Fraction(3, 2),
                damping=0.2,
                channel="range",
                generated_cutoff=6,
                criterion="snr",
                amplitude_scale=amplitude_scale,
                base_noise_sigma=1.0e-6,
                snr_threshold=5.0,
                sideband_order=1,
                projection_nuisance=2,
                samples=4096,
            )
        )
    return tuple(rows)


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _format_design(design: NonlinearDetectabilityDesign) -> str:
    if design.criterion == "relative":
        threshold = f"eta={design.relative_cutoff:.3g}"
    else:
        threshold = (
            f"scale={design.amplitude_scale:.6g}, "
            f"sigma0={design.base_noise_sigma:.3g}, "
            f"snr_min={design.snr_threshold:.3g}"
        )
    return (
        f"{design.label} | rho={_format_fraction(design.rho)} | "
        f"H_gen={design.generated_cutoff} | criterion={design.criterion} | "
        f"{threshold} | "
        f"usable={design.usable_count} "
        f"(lower={design.lower_usable_count}, upper={design.upper_usable_count}) | "
        f"required={design.budget.minimum_frequency_samples} | verdict={design.verdict}"
    )


def audit_report() -> str:
    lines = [
        "Nonlinear second-order detectability audit",
        "",
        "Generated-line coefficient:",
        "- B_k(p,e) = A_k(2p,e), the exact-in-e coefficient of F^2",
        "",
        "Relative cutoff classifications:",
    ]
    for row in default_relative_rows():
        lines.append(f"- {_format_design(row)}")
    lines.append("")
    lines.append("SNR classifications:")
    for row in default_snr_rows():
        lines.append(f"- {_format_design(row)}")
    lines.append("")
    lines.append("Minimum generated-amplitude scale examples:")
    for channel, sigma, nuisance in (
        ("acceleration", 1.0e-3, 1),
        ("range", 1.0e-6, 2),
    ):
        minimum = minimum_generated_amplitude_scale(
            f"minimum nonlinear p=2 e=0.3 rho=3/2 {channel}",
            p=2.0,
            eccentricity=0.3,
            rho=Fraction(3, 2),
            damping=0.2,
            channel=channel,
            generated_cutoff=6,
            base_noise_sigma=sigma,
            snr_threshold=5.0,
            sideband_order=1,
            projection_nuisance=nuisance,
            samples=4096,
        )
        lines.append("- none" if minimum is None else f"- {_format_design(minimum)}")
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
