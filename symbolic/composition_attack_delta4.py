from __future__ import annotations

from dataclasses import dataclass

from enumerate_contractions_delta4 import (
    BLOCK_TYPES,
    BlockType,
    ContractionClass,
    classify_and_label,
    enumerate_contraction_classes,
)


DELTA_MAX = 4

R2_BLOCK_TYPES = (
    BlockType("X", 4, 2, parity=0, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
)
R0A_BLOCK_TYPES = (
    BlockType("S", 5, 0, parity=0),
)
R0B_BLOCK_TYPES = (
    BlockType("DtS", 3, 0, parity=0),
    BlockType("GradS", 3, 1, parity=1),
    BlockType("Dt2S", 4, 0, parity=0),
)

FAMILY_BLOCKS = {
    "R2": R2_BLOCK_TYPES,
    "R0a": R0A_BLOCK_TYPES,
    "R0b": R0B_BLOCK_TYPES,
}

FAMILY_THRESHOLDS = {
    "R2": "w_X >= 4 unless an explicit EX = 0-type rule removes the mixed quadratic witness",
    "R0a": "w_S >= 5, or explicit exclusion of bare S",
    "R0b": "w_D >= 3, or explicit rule removing the mixed derivative witnesses",
}

FAMILY_NAMES = {"X", "S", "DtS", "GradS", "Dt2S"}


@dataclass(frozen=True)
class CompositionCaseSummary:
    combo_id: str
    families: tuple[str, ...]
    imposed_thresholds: tuple[str, ...]
    total_classes: int
    surviving_labels: tuple[str, ...]
    new_surviving_labels: tuple[str, ...]
    sufficient: bool
    smallest_surviving_cross_family: str | None


@dataclass(frozen=True)
class CompositionSummary:
    baseline_survivors: tuple[str, ...]
    cases: tuple[CompositionCaseSummary, ...]


def classify_composition_contraction(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[str, str, str]:
    if not any(name in FAMILY_NAMES for name in signature):
        return classify_and_label(signature, representative)

    if signature == ("DtS",):
        return "DtS", "Proven reducible", "total derivative"
    if signature == ("Dt2S",):
        return "Dt2S", "Proven reducible", "total derivative"
    if signature == ("GradS", "a"):
        return "aGradS", "Proven reducible", "lower-order EOM"

    raise ValueError(
        "Unhandled threshold-budget composition signature "
        f"{signature} with representative {representative}"
    )


def baseline_survivor_labels(max_weight: int = DELTA_MAX) -> tuple[str, ...]:
    classes = enumerate_contraction_classes(max_weight=max_weight)
    return tuple(item.label for item in classes if item.classification == "Surviving candidate")


def enumerate_composition_classes(
    families: tuple[str, ...],
    max_weight: int = DELTA_MAX,
) -> tuple[ContractionClass, ...]:
    family_blocks = tuple(block for family in families for block in FAMILY_BLOCKS[family])
    return enumerate_contraction_classes(
        max_weight=max_weight,
        block_types=BLOCK_TYPES + family_blocks,
        classifier=classify_composition_contraction,
        require_even_parity=True,
    )


def composition_case_summary(
    families: tuple[str, ...],
    max_weight: int = DELTA_MAX,
) -> CompositionCaseSummary:
    classes = enumerate_composition_classes(families, max_weight=max_weight)
    survivors = tuple(item for item in classes if item.classification == "Surviving candidate")
    baseline = set(baseline_survivor_labels(max_weight=max_weight))
    new_survivors = tuple(item for item in survivors if item.label not in baseline)
    smallest = min(new_survivors, key=lambda item: (item.weight, item.label), default=None)
    return CompositionCaseSummary(
        combo_id="+".join(families),
        families=families,
        imposed_thresholds=tuple(FAMILY_THRESHOLDS[family] for family in families),
        total_classes=len(classes),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=tuple(item.label for item in new_survivors),
        sufficient=not new_survivors,
        smallest_surviving_cross_family=None if smallest is None else smallest.label,
    )


def composition_summary(max_weight: int = DELTA_MAX) -> CompositionSummary:
    cases = (
        composition_case_summary(("R2", "R0a"), max_weight=max_weight),
        composition_case_summary(("R2", "R0b"), max_weight=max_weight),
        composition_case_summary(("R0a", "R0b"), max_weight=max_weight),
        composition_case_summary(("R2", "R0a", "R0b"), max_weight=max_weight),
    )
    return CompositionSummary(
        baseline_survivors=baseline_survivor_labels(max_weight=max_weight),
        cases=cases,
    )


def composition_report(max_weight: int = DELTA_MAX) -> str:
    summary = composition_summary(max_weight=max_weight)
    lines = [
        "combo_id\tfamilies\tresult\tsmallest_surviving_cross_family\timposed_thresholds",
    ]
    for case in summary.cases:
        lines.append(
            "\t".join(
                (
                    case.combo_id,
                    ",".join(case.families),
                    "sufficient" if case.sufficient else "not sufficient",
                    case.smallest_surviving_cross_family or "none",
                    " | ".join(case.imposed_thresholds),
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(composition_report())


if __name__ == "__main__":
    main()
