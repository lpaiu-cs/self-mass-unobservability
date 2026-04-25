from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import cos, isfinite, pi, sin, sqrt

from resonant_comparator_audit import ResonanceBudget, line_shape_budget


@dataclass(frozen=True)
class WeightedHarmonic:
    k: int
    forcing_coefficient: float
    transfer_magnitude: float
    projection_magnitude: float
    observable_amplitude: float
    relative_amplitude: float
    usable: bool


@dataclass(frozen=True)
class AmplitudeDesign:
    label: str
    p: float
    eccentricity: float
    rho: Fraction
    damping: float
    channel: str
    harmonic_cutoff: int
    relative_cutoff: float
    polynomial_order: int
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


def exact_in_e_coefficient_formula() -> str:
    return (
        "A_k(p,e) = pi^-1 int_0^(2pi) "
        "(1-e cos E)^(1-p) cos(k(E-e sin E)) dE"
    )


def forcing_cosine_coefficient(
    p: float,
    eccentricity: float,
    harmonic: int,
    samples: int = 8192,
) -> float:
    _check_positive("p", p)
    _check_positive("harmonic", harmonic)
    if not 0 <= eccentricity < 1:
        raise ValueError("eccentricity must satisfy 0 <= e < 1")
    _check_positive("samples", samples)
    total = 0.0
    for item in range(samples):
        eccentric_anomaly = 2 * pi * (item + 0.5) / samples
        mean_anomaly = eccentric_anomaly - eccentricity * sin(eccentric_anomaly)
        weight = (1 - eccentricity * cos(eccentric_anomaly)) ** (1 - p)
        total += weight * cos(harmonic * mean_anomaly)
    return 2 * total / samples


def forcing_coefficients(
    p: float,
    eccentricity: float,
    harmonic_cutoff: int,
    samples: int = 8192,
) -> tuple[float, ...]:
    _check_positive("harmonic_cutoff", harmonic_cutoff)
    return tuple(
        forcing_cosine_coefficient(p, eccentricity, harmonic, samples=samples)
        for harmonic in range(1, harmonic_cutoff + 1)
    )


def second_order_transfer_magnitude(harmonic: int, rho: Fraction, damping: float) -> float:
    _check_positive("harmonic", harmonic)
    _check_positive("rho", float(rho))
    _check_nonnegative("damping", damping)
    rho_sq = float(rho) ** 2
    denominator = sqrt((rho_sq - harmonic**2) ** 2 + (damping * harmonic) ** 2)
    if denominator == 0:
        return float("inf")
    return rho_sq / denominator


def projection_magnitude(
    channel: str,
    harmonic: int,
    kappa_ratio: float = 10.0,
) -> float:
    _check_positive("harmonic", harmonic)
    if channel == "acceleration":
        return 1.0
    if channel == "range":
        denominator = kappa_ratio**2 - harmonic**2
        if denominator == 0:
            return float("inf")
        return abs(1 / denominator)
    raise ValueError(f"unknown channel: {channel}")


def weighted_harmonics(
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    harmonic_cutoff: int,
    relative_cutoff: float,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> tuple[WeightedHarmonic, ...]:
    _check_nonnegative("relative_cutoff", relative_cutoff)
    coefficients = forcing_coefficients(p, eccentricity, harmonic_cutoff, samples=samples)
    raw_rows: list[tuple[int, float, float, float, float]] = []
    for index, coefficient in enumerate(coefficients, start=1):
        transfer = second_order_transfer_magnitude(index, rho, damping)
        projection = projection_magnitude(channel, index, kappa_ratio=kappa_ratio)
        amplitude = abs(coefficient) * transfer * projection
        raw_rows.append((index, coefficient, transfer, projection, amplitude))
    finite_amplitudes = [row[-1] for row in raw_rows if isfinite(row[-1])]
    scale = max(finite_amplitudes, default=0.0)
    rows: list[WeightedHarmonic] = []
    for harmonic, coefficient, transfer, projection, amplitude in raw_rows:
        relative = 0.0 if scale == 0 else amplitude / scale
        rows.append(
            WeightedHarmonic(
                k=harmonic,
                forcing_coefficient=coefficient,
                transfer_magnitude=transfer,
                projection_magnitude=projection,
                observable_amplitude=amplitude,
                relative_amplitude=relative,
                usable=isfinite(amplitude) and relative >= relative_cutoff and amplitude > 0,
            )
        )
    return tuple(rows)


def classify_amplitude_design(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    harmonic_cutoff: int,
    relative_cutoff: float,
    polynomial_order: int,
    projection_nuisance: int,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> AmplitudeDesign:
    budget = line_shape_budget(polynomial_order, projection_nuisance)
    rows = weighted_harmonics(
        p,
        eccentricity,
        rho,
        damping,
        channel,
        harmonic_cutoff,
        relative_cutoff,
        kappa_ratio=kappa_ratio,
        samples=samples,
    )
    usable = tuple(row for row in rows if row.usable)
    lower = tuple(row for row in usable if row.k < rho)
    upper = tuple(row for row in usable if row.k > rho)
    if len(usable) < budget.minimum_frequency_samples:
        verdict = "amplitude-underbudget-no-go"
    elif not lower or not upper:
        verdict = "no-usable-resonance-bracket"
    else:
        verdict = "amplitude-budget-breaking"
    return AmplitudeDesign(
        label=label,
        p=p,
        eccentricity=eccentricity,
        rho=rho,
        damping=damping,
        channel=channel,
        harmonic_cutoff=harmonic_cutoff,
        relative_cutoff=relative_cutoff,
        polynomial_order=polynomial_order,
        projection_nuisance=projection_nuisance,
        usable_count=len(usable),
        lower_usable_count=len(lower),
        upper_usable_count=len(upper),
        budget=budget,
        verdict=verdict,
    )


def minimum_amplitude_cutoff_design(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    relative_cutoff: float,
    polynomial_order: int,
    projection_nuisance: int,
    max_harmonic_cutoff: int = 20,
    kappa_ratio: float = 10.0,
    samples: int = 8192,
) -> AmplitudeDesign | None:
    for cutoff in range(1, max_harmonic_cutoff + 1):
        design = classify_amplitude_design(
            label,
            p,
            eccentricity,
            rho,
            damping,
            channel,
            cutoff,
            relative_cutoff,
            polynomial_order,
            projection_nuisance,
            kappa_ratio=kappa_ratio,
            samples=samples,
        )
        if design.verdict == "amplitude-budget-breaking":
            return design
    return None


def default_design_rows() -> tuple[AmplitudeDesign, ...]:
    rows: list[AmplitudeDesign] = []
    for eccentricity in (0.1, 0.3, 0.6):
        for channel, nuisance in (("acceleration", 1), ("range", 2)):
            rows.append(
                classify_amplitude_design(
                    f"p=2 e={eccentricity} rho=3/2 {channel}",
                    p=2.0,
                    eccentricity=eccentricity,
                    rho=Fraction(3, 2),
                    damping=0.2,
                    channel=channel,
                    harmonic_cutoff=6,
                    relative_cutoff=1.0e-3,
                    polynomial_order=1,
                    projection_nuisance=nuisance,
                    samples=4096,
                )
            )
    return tuple(rows)


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _format_design(design: AmplitudeDesign) -> str:
    return (
        f"{design.label} | rho={_format_fraction(design.rho)} | "
        f"H={design.harmonic_cutoff} | usable={design.usable_count} "
        f"(lower={design.lower_usable_count}, upper={design.upper_usable_count}) | "
        f"required={design.budget.minimum_frequency_samples} | "
        f"cutoff={design.relative_cutoff} | verdict={design.verdict}"
    )


def audit_report() -> str:
    lines = [
        "Amplitude-weighted resonant design audit",
        "",
        "Exact-in-e coefficient formula:",
        f"- {exact_in_e_coefficient_formula()}",
        "",
        "Default classifications:",
    ]
    for row in default_design_rows():
        lines.append(f"- {_format_design(row)}")
    minimum = minimum_amplitude_cutoff_design(
        "minimum p=2 e=0.3 rho=3/2 acceleration",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        relative_cutoff=1.0e-3,
        polynomial_order=1,
        projection_nuisance=1,
        max_harmonic_cutoff=12,
        samples=4096,
    )
    lines.append("")
    lines.append("Minimum example:")
    lines.append("- none" if minimum is None else f"- {_format_design(minimum)}")
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
