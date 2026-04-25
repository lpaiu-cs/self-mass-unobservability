from __future__ import annotations

from dataclasses import dataclass

from sample_budget_audit import DesignClassification, classify_design, compute_sample_budget


@dataclass(frozen=True)
class OrbitalHarmonicDesign:
    linear_harmonics: int
    sideband_harmonic_cutoff: int
    clean_sideband_pairs: tuple[tuple[int, int, int], ...]
    classification: DesignClassification


def _check_positive(name: str, value: int) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def clean_sum_sideband_pairs(
    linear_harmonics: int,
    sideband_harmonic_cutoff: int,
) -> tuple[tuple[int, int, int], ...]:
    _check_positive("linear_harmonics", linear_harmonics)
    _check_positive("sideband_harmonic_cutoff", sideband_harmonic_cutoff)
    pairs: list[tuple[int, int, int]] = []
    for left in range(1, linear_harmonics + 1):
        for right in range(left, linear_harmonics + 1):
            output = left + right
            if linear_harmonics < output <= sideband_harmonic_cutoff:
                pairs.append((left, right, output))
    return tuple(pairs)


def classify_orbital_harmonic_design(
    linear_harmonics: int,
    sideband_harmonic_cutoff: int,
    linear_order: int,
    sideband_degree: int,
    projection_nuisance: int = 0,
) -> OrbitalHarmonicDesign:
    pairs = clean_sum_sideband_pairs(linear_harmonics, sideband_harmonic_cutoff)
    classification = classify_design(
        f"orbital H={linear_harmonics} R={sideband_harmonic_cutoff}",
        linear_harmonics,
        len(pairs),
        linear_order,
        sideband_degree,
        projection_nuisance,
    )
    return OrbitalHarmonicDesign(
        linear_harmonics=linear_harmonics,
        sideband_harmonic_cutoff=sideband_harmonic_cutoff,
        clean_sideband_pairs=pairs,
        classification=classification,
    )


def find_minimum_orbital_design(
    linear_order: int,
    sideband_degree: int,
    projection_nuisance: int = 0,
    max_linear_harmonics: int = 12,
) -> OrbitalHarmonicDesign | None:
    best: OrbitalHarmonicDesign | None = None
    for linear_harmonics in range(1, max_linear_harmonics + 1):
        for cutoff in range(linear_harmonics + 1, 2 * linear_harmonics + 1):
            design = classify_orbital_harmonic_design(
                linear_harmonics,
                cutoff,
                linear_order,
                sideband_degree,
                projection_nuisance,
            )
            if design.classification.verdict != "budget-breaking":
                continue
            if best is None:
                best = design
                continue
            current_key = (
                design.linear_harmonics + design.sideband_harmonic_cutoff,
                design.linear_harmonics,
                design.sideband_harmonic_cutoff,
            )
            best_key = (
                best.linear_harmonics + best.sideband_harmonic_cutoff,
                best.linear_harmonics,
                best.sideband_harmonic_cutoff,
            )
            if current_key < best_key:
                best = design
    return best


def exact_in_e_clean_new_line_count() -> int:
    return 0


def comparator_cases() -> tuple[tuple[str, int, int, int], ...]:
    return (
        ("trivial static sideband", 0, 0, 0),
        ("first-derivative linear + constant sideband", 1, 0, 0),
        ("first-derivative linear + linear sideband", 1, 1, 0),
        ("first-derivative linear + linear sideband + projection nuisance", 1, 1, 2),
        ("first-derivative linear + quadratic sideband", 1, 2, 0),
    )


def minimum_design_rows() -> tuple[OrbitalHarmonicDesign, ...]:
    rows: list[OrbitalHarmonicDesign] = []
    for _, linear_order, sideband_degree, projection_nuisance in comparator_cases():
        design = find_minimum_orbital_design(
            linear_order,
            sideband_degree,
            projection_nuisance,
        )
        if design is not None:
            rows.append(design)
    return tuple(rows)


def _format_design(label: str, design: OrbitalHarmonicDesign | None) -> str:
    if design is None:
        return f"{label}: no budget-breaking design found in the search window"
    budget = compute_sample_budget(
        design.classification.budget.linear_order,
        design.classification.budget.sideband_degree,
        design.classification.budget.projection_nuisance,
    )
    return (
        f"{label}: H={design.linear_harmonics}, R={design.sideband_harmonic_cutoff}, "
        f"clean_pairs={len(design.clean_sideband_pairs)}, "
        f"required=({budget.minimum_linear_samples}, {budget.minimum_sideband_pairs})"
    )


def audit_report() -> str:
    lines = [
        "Richer orbital harmonic budget audit",
        "",
        "Current small-e case:",
        f"- {_format_design('N=1 M=1 K=0 minimum', find_minimum_orbital_design(1, 1, 0))}",
        "",
        "Comparator minimum designs:",
    ]
    for label, linear_order, sideband_degree, projection_nuisance in comparator_cases():
        lines.append(
            "- "
            + _format_design(
                label,
                find_minimum_orbital_design(
                    linear_order,
                    sideband_degree,
                    projection_nuisance,
                ),
            )
        )
    lines.extend(
        [
            "",
            "Exact-in-e clean new-line count:",
            f"- {exact_in_e_clean_new_line_count()}",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
