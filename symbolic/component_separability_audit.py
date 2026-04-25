from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from amplitude_weighted_resonant_design import forcing_cosine_coefficient
from nonlinear_second_order_detectability import (
    generated_forcing_coefficient,
    generated_line_rows,
)
from resonant_comparator_audit import line_shape_budget


@dataclass(frozen=True)
class JointModelParameters:
    q_static: float = 0.2
    q_linear: float = 1.0
    q_generated: float = 0.4
    rho: float = 1.5
    damping: float = 0.2
    kappa_ratio: float = 10.0


@dataclass(frozen=True)
class SeparabilityDesign:
    label: str
    p: float
    eccentricity: float
    rho: Fraction
    damping: float
    channel: str
    harmonic_cutoff: int
    criterion: str
    sideband_order: int
    projection_nuisance: int
    usable_harmonics: tuple[int, ...]
    generated_budget_required: int
    include_range_projection_nuisance: bool
    parameter_count: int
    full_rank: int
    rank_without_generated: int
    generated_adds_rank: bool
    verdict: str


def _check_positive(name: str, value: float) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def _denominator(harmonic: int, rho: float, damping: float) -> complex:
    return complex(rho**2 - harmonic**2, damping * harmonic)


def _projection(channel: str, harmonic: int, kappa_ratio: float) -> complex:
    if channel == "acceleration":
        return 1.0 + 0.0j
    if channel == "range":
        return 1 / (kappa_ratio**2 - harmonic**2)
    raise ValueError(f"unknown channel: {channel}")


def joint_observable(
    harmonic: int,
    p: float,
    eccentricity: float,
    channel: str,
    params: JointModelParameters,
    samples: int = 8192,
) -> complex:
    a_k = forcing_cosine_coefficient(
        p,
        eccentricity,
        harmonic,
        samples=samples,
    )
    b_k = generated_forcing_coefficient(
        p,
        eccentricity,
        harmonic,
        samples=samples,
    )
    denominator = _denominator(harmonic, params.rho, params.damping)
    projection = _projection(channel, harmonic, params.kappa_ratio)
    linear_internal = params.rho**2 / denominator
    generated_internal = 1 / denominator
    return projection * (
        params.q_static * a_k
        + params.q_linear * a_k * linear_internal
        + params.q_generated * b_k * generated_internal
    )


def _parameter_vector(
    params: JointModelParameters,
    include_range_projection_nuisance: bool,
) -> np.ndarray:
    values = [
        params.q_static,
        params.q_linear,
        params.q_generated,
        params.rho,
        params.damping,
    ]
    if include_range_projection_nuisance:
        values.append(params.kappa_ratio)
    return np.array(values, dtype=float)


def _params_from_vector(
    values: np.ndarray,
    include_range_projection_nuisance: bool,
) -> JointModelParameters:
    kappa_ratio = 10.0
    if include_range_projection_nuisance:
        kappa_ratio = float(values[5])
    return JointModelParameters(
        q_static=float(values[0]),
        q_linear=float(values[1]),
        q_generated=float(values[2]),
        rho=float(values[3]),
        damping=float(values[4]),
        kappa_ratio=kappa_ratio,
    )


def real_observation_vector(
    harmonics: tuple[int, ...],
    p: float,
    eccentricity: float,
    channel: str,
    params: JointModelParameters,
    samples: int = 8192,
) -> np.ndarray:
    values: list[float] = []
    for harmonic in harmonics:
        observation = joint_observable(
            harmonic,
            p,
            eccentricity,
            channel,
            params,
            samples=samples,
        )
        values.extend((observation.real, observation.imag))
    return np.array(values, dtype=float)


def jacobian_matrix(
    harmonics: tuple[int, ...],
    p: float,
    eccentricity: float,
    channel: str,
    params: JointModelParameters,
    include_range_projection_nuisance: bool = False,
    step: float = 1.0e-6,
    samples: int = 8192,
) -> np.ndarray:
    base = _parameter_vector(params, include_range_projection_nuisance)
    columns: list[np.ndarray] = []
    for index, value in enumerate(base):
        scale = step * max(1.0, abs(value))
        plus = base.copy()
        minus = base.copy()
        plus[index] += scale
        minus[index] -= scale
        plus_params = _params_from_vector(plus, include_range_projection_nuisance)
        minus_params = _params_from_vector(minus, include_range_projection_nuisance)
        plus_vector = real_observation_vector(
            harmonics,
            p,
            eccentricity,
            channel,
            plus_params,
            samples=samples,
        )
        minus_vector = real_observation_vector(
            harmonics,
            p,
            eccentricity,
            channel,
            minus_params,
            samples=samples,
        )
        columns.append((plus_vector - minus_vector) / (2 * scale))
    return np.column_stack(columns)


def usable_harmonics(
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    harmonic_cutoff: int,
    criterion: str = "relative",
    relative_cutoff: float = 1.0e-3,
    samples: int = 8192,
) -> tuple[int, ...]:
    _check_positive("harmonic_cutoff", harmonic_cutoff)
    rows = generated_line_rows(
        p,
        eccentricity,
        rho,
        damping,
        channel,
        harmonic_cutoff,
        criterion=criterion,
        relative_cutoff=relative_cutoff,
        samples=samples,
    )
    return tuple(row.k for row in rows if row.usable)


def classify_separability_design(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    harmonic_cutoff: int,
    criterion: str = "relative",
    relative_cutoff: float = 1.0e-3,
    sideband_order: int = 1,
    projection_nuisance: int = 0,
    include_range_projection_nuisance: bool = False,
    rank_tolerance: float = 1.0e-8,
    samples: int = 8192,
) -> SeparabilityDesign:
    harmonics = usable_harmonics(
        p,
        eccentricity,
        rho,
        damping,
        channel,
        harmonic_cutoff,
        criterion=criterion,
        relative_cutoff=relative_cutoff,
        samples=samples,
    )
    generated_budget = line_shape_budget(sideband_order, projection_nuisance)
    params = JointModelParameters(rho=float(rho), damping=damping)
    if not harmonics:
        full_rank = 0
        rank_without_generated = 0
        parameter_count = len(_parameter_vector(params, include_range_projection_nuisance))
    else:
        matrix = jacobian_matrix(
            harmonics,
            p,
            eccentricity,
            channel,
            params,
            include_range_projection_nuisance=include_range_projection_nuisance,
            samples=samples,
        )
        parameter_count = matrix.shape[1]
        full_rank = int(np.linalg.matrix_rank(matrix, tol=rank_tolerance))
        reduced = np.delete(matrix, 2, axis=1)
        rank_without_generated = int(np.linalg.matrix_rank(reduced, tol=rank_tolerance))
    generated_adds_rank = full_rank > rank_without_generated
    if not generated_adds_rank:
        verdict = "generated-component-degenerate"
    elif full_rank < parameter_count:
        verdict = "partial-separability-with-nuisance-degeneracy"
    elif len(harmonics) < generated_budget.minimum_frequency_samples:
        verdict = "component-separable-but-budget-underdesigned"
    else:
        verdict = "component-separable-and-budget-ready"
    return SeparabilityDesign(
        label=label,
        p=p,
        eccentricity=eccentricity,
        rho=rho,
        damping=damping,
        channel=channel,
        harmonic_cutoff=harmonic_cutoff,
        criterion=criterion,
        sideband_order=sideband_order,
        projection_nuisance=projection_nuisance,
        usable_harmonics=harmonics,
        generated_budget_required=generated_budget.minimum_frequency_samples,
        include_range_projection_nuisance=include_range_projection_nuisance,
        parameter_count=parameter_count,
        full_rank=full_rank,
        rank_without_generated=rank_without_generated,
        generated_adds_rank=generated_adds_rank,
        verdict=verdict,
    )


def minimum_separable_cutoff(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    max_harmonic_cutoff: int = 12,
    criterion: str = "relative",
    relative_cutoff: float = 1.0e-3,
    sideband_order: int = 1,
    projection_nuisance: int = 0,
    include_range_projection_nuisance: bool = False,
    samples: int = 8192,
) -> SeparabilityDesign | None:
    for cutoff in range(1, max_harmonic_cutoff + 1):
        design = classify_separability_design(
            label,
            p,
            eccentricity,
            rho,
            damping,
            channel,
            cutoff,
            criterion=criterion,
            relative_cutoff=relative_cutoff,
            sideband_order=sideband_order,
            projection_nuisance=projection_nuisance,
            include_range_projection_nuisance=include_range_projection_nuisance,
            samples=samples,
        )
        if design.verdict == "component-separable-and-budget-ready":
            return design
    return None


def default_design_rows() -> tuple[SeparabilityDesign, ...]:
    rows: list[SeparabilityDesign] = []
    for eccentricity in (0.1, 0.3, 0.6):
        rows.append(
            classify_separability_design(
                f"p=2 e={eccentricity} acceleration",
                p=2.0,
                eccentricity=eccentricity,
                rho=Fraction(3, 2),
                damping=0.2,
                channel="acceleration",
                harmonic_cutoff=6,
                relative_cutoff=1.0e-3,
                projection_nuisance=1,
                samples=4096,
            )
        )
        rows.append(
            classify_separability_design(
                f"p=2 e={eccentricity} range fixed-kappa",
                p=2.0,
                eccentricity=eccentricity,
                rho=Fraction(3, 2),
                damping=0.2,
                channel="range",
                harmonic_cutoff=6,
                relative_cutoff=1.0e-3,
                projection_nuisance=2,
                samples=4096,
            )
        )
        rows.append(
            classify_separability_design(
                f"p=2 e={eccentricity} range free-kappa",
                p=2.0,
                eccentricity=eccentricity,
                rho=Fraction(3, 2),
                damping=0.2,
                channel="range",
                harmonic_cutoff=6,
                relative_cutoff=1.0e-3,
                projection_nuisance=2,
                include_range_projection_nuisance=True,
                samples=4096,
            )
        )
    return tuple(rows)


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _format_design(design: SeparabilityDesign) -> str:
    harmonics = ",".join(str(item) for item in design.usable_harmonics) or "none"
    return (
        f"{design.label} | rho={_format_fraction(design.rho)} | "
        f"H={design.harmonic_cutoff} | usable={harmonics} | "
        f"generated_required={design.generated_budget_required} | "
        f"rank={design.full_rank}/{design.parameter_count} | "
        f"rank_without_generated={design.rank_without_generated} | "
        f"generated_adds_rank={design.generated_adds_rank} | verdict={design.verdict}"
    )


def audit_report() -> str:
    lines = [
        "Component separability audit",
        "",
        "Joint model:",
        "- O_k = Lambda_k [Q_Y A_k + Q_L A_k H_2(k) + Q_beta B_k/D_2(k)]",
        "- B_k = A_k(2p,e)",
        "",
        "Default classifications:",
    ]
    for row in default_design_rows():
        lines.append(f"- {_format_design(row)}")
    lines.append("")
    lines.append("Minimum separable cutoff examples:")
    for channel, projection_nuisance in (
        ("acceleration", False),
        ("range", False),
        ("range", True),
    ):
        minimum = minimum_separable_cutoff(
            f"minimum p=2 e=0.3 {channel}",
            p=2.0,
            eccentricity=0.3,
            rho=Fraction(3, 2),
            damping=0.2,
            channel=channel,
            projection_nuisance=1 if channel == "acceleration" else 2,
            include_range_projection_nuisance=projection_nuisance,
            samples=4096,
        )
        lines.append("- none" if minimum is None else f"- {_format_design(minimum)}")
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
