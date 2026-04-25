from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from amplitude_weighted_resonant_design import (
    forcing_cosine_coefficient,
    projection_magnitude,
    second_order_transfer_magnitude,
)
from resonant_comparator_audit import ResonanceBudget, line_shape_budget


@dataclass(frozen=True)
class DetectableHarmonic:
    k: int
    forcing_coefficient: float
    observable_amplitude: float
    noise_sigma: float
    snr: float
    detectable: bool


@dataclass(frozen=True)
class PhysicalDesign:
    label: str
    p: float
    eccentricity: float
    rho: Fraction
    damping: float
    channel: str
    harmonic_cutoff: int
    amplitude_scale: float
    base_noise_sigma: float
    snr_threshold: float
    polynomial_order: int
    projection_nuisance: int
    detectable_count: int
    lower_detectable_count: int
    upper_detectable_count: int
    budget: ResonanceBudget
    verdict: str


def _check_positive(name: str, value: float) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def harmonic_noise_sigma(
    base_noise_sigma: float,
    harmonic: int,
    harmonic_noise_slope: float = 0.0,
) -> float:
    _check_positive("base_noise_sigma", base_noise_sigma)
    _check_positive("harmonic", harmonic)
    return base_noise_sigma * harmonic**harmonic_noise_slope


def detectable_harmonics(
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    harmonic_cutoff: int,
    amplitude_scale: float,
    base_noise_sigma: float,
    snr_threshold: float,
    harmonic_noise_slope: float = 0.0,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> tuple[DetectableHarmonic, ...]:
    _check_positive("harmonic_cutoff", harmonic_cutoff)
    _check_positive("amplitude_scale", amplitude_scale)
    _check_positive("snr_threshold", snr_threshold)
    rows: list[DetectableHarmonic] = []
    for harmonic in range(1, harmonic_cutoff + 1):
        coefficient = forcing_cosine_coefficient(
            p,
            eccentricity,
            harmonic,
            samples=samples,
        )
        transfer = second_order_transfer_magnitude(harmonic, rho, damping)
        projection = projection_magnitude(channel, harmonic, kappa_ratio=kappa_ratio)
        observable = amplitude_scale * abs(coefficient) * transfer * projection
        sigma = harmonic_noise_sigma(
            base_noise_sigma,
            harmonic,
            harmonic_noise_slope=harmonic_noise_slope,
        )
        snr = observable / sigma
        rows.append(
            DetectableHarmonic(
                k=harmonic,
                forcing_coefficient=coefficient,
                observable_amplitude=observable,
                noise_sigma=sigma,
                snr=snr,
                detectable=snr >= snr_threshold,
            )
        )
    return tuple(rows)


def classify_physical_design(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    harmonic_cutoff: int,
    amplitude_scale: float,
    base_noise_sigma: float,
    snr_threshold: float,
    polynomial_order: int,
    projection_nuisance: int,
    harmonic_noise_slope: float = 0.0,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> PhysicalDesign:
    budget = line_shape_budget(polynomial_order, projection_nuisance)
    rows = detectable_harmonics(
        p,
        eccentricity,
        rho,
        damping,
        channel,
        harmonic_cutoff,
        amplitude_scale,
        base_noise_sigma,
        snr_threshold,
        harmonic_noise_slope=harmonic_noise_slope,
        kappa_ratio=kappa_ratio,
        samples=samples,
    )
    detectable = tuple(row for row in rows if row.detectable)
    lower = tuple(row for row in detectable if row.k < rho)
    upper = tuple(row for row in detectable if row.k > rho)
    if len(detectable) < budget.minimum_frequency_samples:
        verdict = "snr-underbudget-no-go"
    elif not lower or not upper:
        verdict = "no-detectable-resonance-bracket"
    else:
        verdict = "physically-budget-breaking"
    return PhysicalDesign(
        label=label,
        p=p,
        eccentricity=eccentricity,
        rho=rho,
        damping=damping,
        channel=channel,
        harmonic_cutoff=harmonic_cutoff,
        amplitude_scale=amplitude_scale,
        base_noise_sigma=base_noise_sigma,
        snr_threshold=snr_threshold,
        polynomial_order=polynomial_order,
        projection_nuisance=projection_nuisance,
        detectable_count=len(detectable),
        lower_detectable_count=len(lower),
        upper_detectable_count=len(upper),
        budget=budget,
        verdict=verdict,
    )


def minimum_physical_amplitude_scale(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    harmonic_cutoff: int,
    base_noise_sigma: float,
    snr_threshold: float,
    polynomial_order: int,
    projection_nuisance: int,
    harmonic_noise_slope: float = 0.0,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> PhysicalDesign | None:
    unit_rows = detectable_harmonics(
        p,
        eccentricity,
        rho,
        damping,
        channel,
        harmonic_cutoff,
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
        if row.snr > 0
    )
    for scale in candidate_scales:
        test_scale = scale * (1.0 + 1.0e-12)
        design = classify_physical_design(
            label,
            p,
            eccentricity,
            rho,
            damping,
            channel,
            harmonic_cutoff,
            test_scale,
            base_noise_sigma,
            snr_threshold,
            polynomial_order,
            projection_nuisance,
            harmonic_noise_slope=harmonic_noise_slope,
            kappa_ratio=kappa_ratio,
            samples=samples,
        )
        if design.verdict == "physically-budget-breaking":
            return design
    return None


def default_physical_design_rows() -> tuple[PhysicalDesign, ...]:
    rows: list[PhysicalDesign] = []
    for amplitude_scale in (0.1, 0.3, 1.0):
        rows.append(
            classify_physical_design(
                f"acceleration scale={amplitude_scale}",
                p=2.0,
                eccentricity=0.3,
                rho=Fraction(3, 2),
                damping=0.2,
                channel="acceleration",
                harmonic_cutoff=6,
                amplitude_scale=amplitude_scale,
                base_noise_sigma=1.0e-3,
                snr_threshold=5.0,
                polynomial_order=1,
                projection_nuisance=1,
                samples=4096,
            )
        )
        rows.append(
            classify_physical_design(
                f"range scale={amplitude_scale}",
                p=2.0,
                eccentricity=0.3,
                rho=Fraction(3, 2),
                damping=0.2,
                channel="range",
                harmonic_cutoff=6,
                amplitude_scale=amplitude_scale,
                base_noise_sigma=1.0e-6,
                snr_threshold=5.0,
                polynomial_order=1,
                projection_nuisance=2,
                samples=4096,
            )
        )
    return tuple(rows)


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _format_design(design: PhysicalDesign) -> str:
    return (
        f"{design.label} | rho={_format_fraction(design.rho)} | "
        f"H={design.harmonic_cutoff} | scale={design.amplitude_scale:.6g} | "
        f"sigma0={design.base_noise_sigma:.3g} | detectable={design.detectable_count} "
        f"(lower={design.lower_detectable_count}, upper={design.upper_detectable_count}) | "
        f"required={design.budget.minimum_frequency_samples} | "
        f"snr_min={design.snr_threshold} | verdict={design.verdict}"
    )


def audit_report() -> str:
    lines = [
        "Physical detectability map",
        "",
        "Parameterized noise model:",
        "- SNR_k = O_k / sigma_k",
        "- sigma_k = sigma_0 k^q",
        "- detectable iff SNR_k >= SNR_min",
        "",
        "Default classifications:",
    ]
    for row in default_physical_design_rows():
        lines.append(f"- {_format_design(row)}")
    lines.append("")
    lines.append("Minimum physical scale examples:")
    for channel, sigma, nuisance in (
        ("acceleration", 1.0e-3, 1),
        ("range", 1.0e-6, 2),
    ):
        minimum = minimum_physical_amplitude_scale(
            f"minimum p=2 e=0.3 rho=3/2 {channel}",
            p=2.0,
            eccentricity=0.3,
            rho=Fraction(3, 2),
            damping=0.2,
            channel=channel,
            harmonic_cutoff=6,
            base_noise_sigma=sigma,
            snr_threshold=5.0,
            polynomial_order=1,
            projection_nuisance=nuisance,
            samples=4096,
        )
        lines.append("- none" if minimum is None else f"- {_format_design(minimum)}")
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
