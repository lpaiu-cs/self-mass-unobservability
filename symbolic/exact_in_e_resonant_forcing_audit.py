from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from resonant_comparator_audit import ResonanceBudget, line_shape_budget


@dataclass(frozen=True)
class HarmonicBracket:
    resonance_ratio: Fraction
    harmonic_cutoff: int
    lower_count: int
    upper_count: int
    exact_hit: bool
    bracketed: bool


@dataclass(frozen=True)
class ResonantForcingDesign:
    label: str
    polynomial_order: int
    projection_nuisance: int
    resonance_ratio: Fraction
    minimum_harmonic_cutoff: int | None
    budget: ResonanceBudget
    verdict: str


def _check_positive(name: str, value: int) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def harmonics(harmonic_cutoff: int) -> tuple[int, ...]:
    _check_positive("harmonic_cutoff", harmonic_cutoff)
    return tuple(range(1, harmonic_cutoff + 1))


def harmonic_bracket(
    resonance_ratio: Fraction,
    harmonic_cutoff: int,
) -> HarmonicBracket:
    _check_positive("harmonic_cutoff", harmonic_cutoff)
    if resonance_ratio <= 0:
        raise ValueError("resonance_ratio must be positive")
    lower_count = sum(1 for item in harmonics(harmonic_cutoff) if item < resonance_ratio)
    upper_count = sum(1 for item in harmonics(harmonic_cutoff) if item > resonance_ratio)
    exact_hit = any(item == resonance_ratio for item in harmonics(harmonic_cutoff))
    return HarmonicBracket(
        resonance_ratio=resonance_ratio,
        harmonic_cutoff=harmonic_cutoff,
        lower_count=lower_count,
        upper_count=upper_count,
        exact_hit=exact_hit,
        bracketed=lower_count > 0 and upper_count > 0,
    )


def minimum_cutoff_for_bracket(resonance_ratio: Fraction) -> int | None:
    if resonance_ratio <= 1:
        return None
    if resonance_ratio.denominator == 1:
        return int(resonance_ratio) + 1
    return int(resonance_ratio) + 1


def minimum_cutoff_for_budget_and_bracket(
    resonance_ratio: Fraction,
    polynomial_order: int,
    projection_nuisance: int,
) -> int | None:
    budget = line_shape_budget(polynomial_order, projection_nuisance)
    bracket_cutoff = minimum_cutoff_for_bracket(resonance_ratio)
    if bracket_cutoff is None:
        return None
    return max(budget.minimum_frequency_samples, bracket_cutoff)


def classify_resonant_forcing_design(
    label: str,
    resonance_ratio: Fraction,
    polynomial_order: int,
    projection_nuisance: int,
) -> ResonantForcingDesign:
    budget = line_shape_budget(polynomial_order, projection_nuisance)
    minimum_cutoff = minimum_cutoff_for_budget_and_bracket(
        resonance_ratio,
        polynomial_order,
        projection_nuisance,
    )
    verdict = "no single-system harmonic bracket" if minimum_cutoff is None else "budget-breaking-if-usable"
    return ResonantForcingDesign(
        label=label,
        polynomial_order=polynomial_order,
        projection_nuisance=projection_nuisance,
        resonance_ratio=resonance_ratio,
        minimum_harmonic_cutoff=minimum_cutoff,
        budget=budget,
        verdict=verdict,
    )


def default_design_rows() -> tuple[ResonantForcingDesign, ...]:
    ratios = (
        ("near 3/2 wrap", Fraction(3, 2)),
        ("near 5/2 wrap", Fraction(5, 2)),
        ("integer 4 wrap", Fraction(4, 1)),
        ("sub-fundamental wrap", Fraction(3, 4)),
    )
    comparator_cases = (
        ("acceleration N=1", 1, 1),
        ("range N=1", 1, 2),
        ("range N=2", 2, 2),
    )
    rows: list[ResonantForcingDesign] = []
    for ratio_label, ratio in ratios:
        for comparator_label, order, nuisance in comparator_cases:
            rows.append(
                classify_resonant_forcing_design(
                    f"{ratio_label}; {comparator_label}",
                    ratio,
                    order,
                    nuisance,
                )
            )
    return tuple(rows)


def _format_ratio(ratio: Fraction) -> str:
    if ratio.denominator == 1:
        return str(ratio.numerator)
    return f"{ratio.numerator}/{ratio.denominator}"


def _format_design(row: ResonantForcingDesign) -> str:
    cutoff = "none" if row.minimum_harmonic_cutoff is None else str(row.minimum_harmonic_cutoff)
    return (
        f"{row.label} | rho={_format_ratio(row.resonance_ratio)} | "
        f"N={row.polynomial_order} K={row.projection_nuisance} | "
        f"min_samples={row.budget.minimum_frequency_samples} | "
        f"H_min={cutoff} | verdict={row.verdict}"
    )


def audit_report() -> str:
    lines = [
        "Exact-in-e resonant forcing audit",
        "",
        "Bracket examples:",
        f"- rho=3/2, H=4 -> {harmonic_bracket(Fraction(3, 2), 4)}",
        f"- rho=4, H=4 -> {harmonic_bracket(Fraction(4, 1), 4)}",
        f"- rho=4, H=5 -> {harmonic_bracket(Fraction(4, 1), 5)}",
        "",
        "Design rows:",
    ]
    for row in default_design_rows():
        lines.append(f"- {_format_design(row)}")
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
