from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import sympy as sp

from component_separability_audit import classify_separability_design
from exact_in_e_resonant_forcing_audit import minimum_cutoff_for_bracket
from nonlinear_robustness_map import RobustnessRow, robustness_grid
from resonant_comparator_audit import line_shape_budget


@dataclass(frozen=True)
class BoundaryClassification:
    label: str
    rho: Fraction
    harmonic_cutoff: int
    count_minimum: int
    bracket_minimum: int | None
    analytic_minimum: int | None
    usable_count: int
    bracketed: bool
    component_verdict: str
    analytic_verdict: str


def generated_count_minimum(sideband_order: int, projection_nuisance: int) -> int:
    return line_shape_budget(sideband_order, projection_nuisance).minimum_frequency_samples


def analytic_harmonic_minimum(
    rho: Fraction,
    sideband_order: int,
    projection_nuisance: int,
) -> int | None:
    bracket_minimum = minimum_cutoff_for_bracket(rho)
    if bracket_minimum is None:
        return None
    return max(generated_count_minimum(sideband_order, projection_nuisance), bracket_minimum)


def qbeta_minor_numerator() -> sp.Expr:
    a1, a2, a3 = sp.symbols("A1 A2 A3")
    b1, b2, b3 = sp.symbols("B1 B2 B3")
    d1, d2, d3 = sp.symbols("D1 D2 D3")
    rho_sq = sp.symbols("rho_sq")
    matrix = sp.Matrix(
        [
            [a1, a1 * rho_sq / d1, b1 / d1],
            [a2, a2 * rho_sq / d2, b2 / d2],
            [a3, a3 * rho_sq / d3, b3 / d3],
        ]
    )
    return sp.factor(matrix.det() * d1 * d2 * d3)


def qbeta_minor_vanishes_for_static_proportional_generated() -> sp.Expr:
    a1, a2, a3 = sp.symbols("A1 A2 A3")
    d1, d2, d3 = sp.symbols("D1 D2 D3")
    c = sp.symbols("c")
    numerator = qbeta_minor_numerator()
    return sp.factor(
        numerator.subs(
            {
                sp.Symbol("B1"): c * a1,
                sp.Symbol("B2"): c * a2,
                sp.Symbol("B3"): c * a3,
            }
        )
    )


def _is_bracketed(harmonics: tuple[int, ...], rho: Fraction) -> bool:
    return any(item < rho for item in harmonics) and any(item > rho for item in harmonics)


def classify_boundary(
    label: str,
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    harmonic_cutoff: int,
    sideband_order: int,
    projection_nuisance: int,
    include_range_projection_nuisance: bool = False,
    relative_cutoff: float = 1.0e-3,
    samples: int = 1024,
) -> BoundaryClassification:
    count_minimum = generated_count_minimum(sideband_order, projection_nuisance)
    bracket_minimum = minimum_cutoff_for_bracket(rho)
    analytic_minimum = analytic_harmonic_minimum(rho, sideband_order, projection_nuisance)
    design = classify_separability_design(
        label,
        p=p,
        eccentricity=eccentricity,
        rho=rho,
        damping=damping,
        channel=channel,
        harmonic_cutoff=harmonic_cutoff,
        criterion="relative",
        relative_cutoff=relative_cutoff,
        sideband_order=sideband_order,
        projection_nuisance=projection_nuisance,
        include_range_projection_nuisance=include_range_projection_nuisance,
        samples=samples,
    )
    bracketed = _is_bracketed(design.usable_harmonics, rho)
    if analytic_minimum is None:
        analytic_verdict = "no-single-system-bracket"
    elif harmonic_cutoff < analytic_minimum:
        analytic_verdict = "harmonic-cutoff-under-minimum"
    elif design.verdict == "generated-component-degenerate":
        analytic_verdict = "rank-degenerate"
    elif design.verdict == "partial-separability-with-nuisance-degeneracy":
        analytic_verdict = "nuisance-rank-degenerate"
    elif len(design.usable_harmonics) < count_minimum:
        analytic_verdict = "usable-count-underbudget"
    elif not bracketed:
        analytic_verdict = "usable-bracket-failure"
    elif design.verdict == "component-separable-but-budget-underdesigned":
        analytic_verdict = "component-separable-but-budget-underdesigned"
    elif design.verdict == "component-separable-and-budget-ready":
        analytic_verdict = "analytic-positive"
    else:
        analytic_verdict = design.verdict
    return BoundaryClassification(
        label=label,
        rho=rho,
        harmonic_cutoff=harmonic_cutoff,
        count_minimum=count_minimum,
        bracket_minimum=bracket_minimum,
        analytic_minimum=analytic_minimum,
        usable_count=len(design.usable_harmonics),
        bracketed=bracketed,
        component_verdict=design.verdict,
        analytic_verdict=analytic_verdict,
    )


def classify_from_robustness_row(row: RobustnessRow) -> BoundaryClassification:
    projection_nuisance = 1 if row.channel == "acceleration" else 2
    include_projection = row.projection_mode == "range-free-kappa"
    return classify_boundary(
        row.projection_mode,
        p=row.p,
        eccentricity=row.eccentricity,
        rho=row.rho,
        damping=row.damping,
        channel=row.channel,
        harmonic_cutoff=row.harmonic_cutoff,
        sideband_order=1,
        projection_nuisance=projection_nuisance,
        include_range_projection_nuisance=include_projection,
        samples=512,
    )


def boundary_matches_robustness(row: RobustnessRow) -> bool:
    boundary = classify_from_robustness_row(row)
    if row.verdict == "robust-positive":
        return boundary.analytic_verdict == "analytic-positive"
    if row.verdict == "component-separable-but-budget-underdesigned":
        return boundary.analytic_verdict in {
            "usable-count-underbudget",
            "component-separable-but-budget-underdesigned",
        }
    if row.verdict == "generated-component-degenerate":
        return boundary.analytic_verdict == "rank-degenerate"
    return boundary.analytic_verdict == row.verdict


def default_boundary_rows() -> tuple[BoundaryClassification, ...]:
    return (
        classify_boundary(
            "low-e acceleration",
            p=2.0,
            eccentricity=0.1,
            rho=Fraction(3, 2),
            damping=0.2,
            channel="acceleration",
            harmonic_cutoff=6,
            sideband_order=1,
            projection_nuisance=1,
            samples=1024,
        ),
        classify_boundary(
            "moderate-e acceleration",
            p=2.0,
            eccentricity=0.3,
            rho=Fraction(3, 2),
            damping=0.2,
            channel="acceleration",
            harmonic_cutoff=6,
            sideband_order=1,
            projection_nuisance=1,
            samples=1024,
        ),
        classify_boundary(
            "moderate-e range free-kappa",
            p=2.0,
            eccentricity=0.3,
            rho=Fraction(3, 2),
            damping=0.2,
            channel="range",
            harmonic_cutoff=6,
            sideband_order=1,
            projection_nuisance=2,
            include_range_projection_nuisance=True,
            samples=1024,
        ),
    )


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _format_boundary(row: BoundaryClassification) -> str:
    return (
        f"{row.label} | rho={_format_fraction(row.rho)} | H={row.harmonic_cutoff} | "
        f"H_count={row.count_minimum} | H_bracket={row.bracket_minimum} | "
        f"H_min={row.analytic_minimum} | usable={row.usable_count} | "
        f"bracketed={row.bracketed} | component={row.component_verdict} | "
        f"verdict={row.analytic_verdict}"
    )


def audit_report() -> str:
    rows = robustness_grid(samples=512)
    matches = sum(1 for row in rows if boundary_matches_robustness(row))
    lines = [
        "Analytic phase-boundary audit",
        "",
        "Boundary rules:",
        "- H_count = M + 2 + K_Lambda",
        "- H_bracket(rho) exists iff rho > 1",
        "- H_min = max(H_count, H_bracket)",
        "- usable generated harmonics must satisfy the amplitude cutoff and bracket rho",
        "- Q_beta must add Jacobian rank after shared nuisance parameters are included",
        "",
        "Q_beta amplitude minor numerator:",
        f"- {sp.sstr(qbeta_minor_numerator())}",
        "Static-proportional generated collapse check:",
        f"- {sp.sstr(qbeta_minor_vanishes_for_static_proportional_generated())}",
        "",
        "Example boundary rows:",
    ]
    for row in default_boundary_rows():
        lines.append(f"- {_format_boundary(row)}")
    lines.append("")
    lines.append("Grid consistency:")
    lines.append(f"- boundary matches robustness verdicts: {matches}/{len(rows)}")
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
