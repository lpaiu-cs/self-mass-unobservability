from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SampleBudget:
    linear_order: int
    sideband_degree: int
    projection_nuisance: int
    linear_parameter_budget: int
    sideband_parameter_budget: int
    total_shared_budget: int
    minimum_linear_samples: int
    minimum_sideband_pairs: int
    minimum_total_samples: int


@dataclass(frozen=True)
class DesignClassification:
    name: str
    linear_samples: int
    sideband_pairs: int
    budget: SampleBudget
    beats_linear_budget: bool
    beats_sideband_budget: bool
    beats_total_shared_budget: bool
    verdict: str


def _check_nonnegative(name: str, value: int) -> None:
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")


def bivariate_monomial_count(total_degree: int) -> int:
    _check_nonnegative("total_degree", total_degree)
    return (total_degree + 1) * (total_degree + 2) // 2


def compute_sample_budget(
    linear_order: int,
    sideband_degree: int,
    projection_nuisance: int = 0,
) -> SampleBudget:
    _check_nonnegative("linear_order", linear_order)
    _check_nonnegative("sideband_degree", sideband_degree)
    _check_nonnegative("projection_nuisance", projection_nuisance)

    linear_budget = linear_order + 1
    sideband_budget = bivariate_monomial_count(sideband_degree)
    total_budget = linear_budget + sideband_budget + projection_nuisance

    minimum_linear = linear_budget + 1
    individual_minimum_sideband = sideband_budget + 1
    minimum_total = max(
        linear_budget + sideband_budget + 2,
        total_budget + 1,
    )
    canonical_sideband = minimum_total - minimum_linear
    minimum_sideband = max(individual_minimum_sideband, canonical_sideband)
    minimum_total = minimum_linear + minimum_sideband

    return SampleBudget(
        linear_order=linear_order,
        sideband_degree=sideband_degree,
        projection_nuisance=projection_nuisance,
        linear_parameter_budget=linear_budget,
        sideband_parameter_budget=sideband_budget,
        total_shared_budget=total_budget,
        minimum_linear_samples=minimum_linear,
        minimum_sideband_pairs=minimum_sideband,
        minimum_total_samples=minimum_total,
    )


def classify_design(
    name: str,
    linear_samples: int,
    sideband_pairs: int,
    linear_order: int,
    sideband_degree: int,
    projection_nuisance: int = 0,
) -> DesignClassification:
    _check_nonnegative("linear_samples", linear_samples)
    _check_nonnegative("sideband_pairs", sideband_pairs)
    budget = compute_sample_budget(linear_order, sideband_degree, projection_nuisance)
    beats_linear = linear_samples > budget.linear_parameter_budget
    beats_sideband = sideband_pairs > budget.sideband_parameter_budget
    beats_total = (
        linear_samples + sideband_pairs
    ) > budget.total_shared_budget
    verdict = (
        "budget-breaking"
        if beats_linear and beats_sideband and beats_total
        else "underbudget-no-go"
    )
    return DesignClassification(
        name=name,
        linear_samples=linear_samples,
        sideband_pairs=sideband_pairs,
        budget=budget,
        beats_linear_budget=beats_linear,
        beats_sideband_budget=beats_sideband,
        beats_total_shared_budget=beats_total,
        verdict=verdict,
    )


def orbital_design() -> tuple[int, int]:
    return 2, 1


def two_tone_design(include_difference: bool = False) -> tuple[int, int]:
    linear_samples = 2
    sum_sideband_pairs = 3
    difference_pair = 1 if include_difference else 0
    return linear_samples, sum_sideband_pairs + difference_pair


def default_comparator_grid() -> tuple[tuple[str, int, int, int], ...]:
    return (
        ("constant linear + constant sideband, no projection nuisance", 0, 0, 0),
        ("first-derivative linear + constant sideband, no projection nuisance", 1, 0, 0),
        ("first-derivative linear + linear sideband, no projection nuisance", 1, 1, 0),
        ("first-derivative linear + quadratic sideband, no projection nuisance", 1, 2, 0),
        ("first-derivative linear + linear sideband, two projection nuisances", 1, 1, 2),
    )


def classification_rows() -> tuple[DesignClassification, ...]:
    rows: list[DesignClassification] = []
    orbital_linear, orbital_sideband = orbital_design()
    two_tone_linear, two_tone_sideband = two_tone_design(include_difference=False)
    two_tone_diff_linear, two_tone_diff_sideband = two_tone_design(include_difference=True)
    for _, linear_order, sideband_degree, projection_nuisance in default_comparator_grid():
        rows.append(
            classify_design(
                "orbital n,2n -> 3n",
                orbital_linear,
                orbital_sideband,
                linear_order,
                sideband_degree,
                projection_nuisance,
            )
        )
        rows.append(
            classify_design(
                "two-tone sum-only",
                two_tone_linear,
                two_tone_sideband,
                linear_order,
                sideband_degree,
                projection_nuisance,
            )
        )
        rows.append(
            classify_design(
                "two-tone sum+difference",
                two_tone_diff_linear,
                two_tone_diff_sideband,
                linear_order,
                sideband_degree,
                projection_nuisance,
            )
        )
    return tuple(rows)


def _format_classification(row: DesignClassification) -> str:
    budget = row.budget
    return (
        f"{row.name} | N={budget.linear_order} M={budget.sideband_degree} "
        f"K={budget.projection_nuisance} | available=({row.linear_samples},"
        f" {row.sideband_pairs}) | required=({budget.minimum_linear_samples},"
        f" {budget.minimum_sideband_pairs}) | verdict={row.verdict}"
    )


def audit_report() -> str:
    example = compute_sample_budget(linear_order=1, sideband_degree=1, projection_nuisance=0)
    lines = [
        "Sample-budget audit",
        "",
        "Theorem budget for N=1, M=1, K=0:",
        (
            f"- linear>{example.linear_parameter_budget}, "
            f"sideband>{example.sideband_parameter_budget}, "
            f"total>{example.total_shared_budget}"
        ),
        (
            f"- canonical minimum: linear={example.minimum_linear_samples}, "
            f"sideband={example.minimum_sideband_pairs}, "
            f"total={example.minimum_total_samples}"
        ),
        "",
        "Case classifications:",
    ]
    for row in classification_rows():
        lines.append(f"- {_format_classification(row)}")
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
