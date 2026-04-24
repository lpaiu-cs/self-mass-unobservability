from __future__ import annotations

from dataclasses import dataclass

from enumerate_contractions_delta4 import BlockType, ContractionClass, enumerate_contraction_classes
from eb_sector_delta4 import EB_BLOCK_TYPES, classify_eb_contraction, eb_summary


MAX_WEIGHT = 4

SHIFT_SCALAR_BLOCK_TYPES = (
    BlockType("DtS", 2, 0, parity=0),
    BlockType("GradS", 2, 1, parity=1),
    BlockType("Dt2S", 3, 0, parity=0),
)

SHIFT_SCALAR_FAMILY_NAMES = {block.name for block in SHIFT_SCALAR_BLOCK_TYPES}
EBDS_BLOCK_TYPES = EB_BLOCK_TYPES + SHIFT_SCALAR_BLOCK_TYPES


@dataclass(frozen=True)
class ShiftScalarSummary:
    total_classes: int
    surviving_labels: tuple[str, ...]
    new_surviving_labels: tuple[str, ...]
    first_new_weight: int | None
    first_new_labels: tuple[str, ...]
    canonical_new_survivor: str | None


def classify_shift_scalar_contraction(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[str, str, str]:
    if not any(name in SHIFT_SCALAR_FAMILY_NAMES for name in signature):
        return classify_eb_contraction(signature, representative)

    if signature == ("DtS",):
        return "DtS", "Proven reducible", "total derivative"
    if signature == ("Dt2S",):
        return "Dt2S", "Proven reducible", "total derivative"
    if signature == ("GradS", "a"):
        return "aGradS", "Proven reducible", "lower-order EOM"
    if signature == ("B", "B", "DtS"):
        return "DtS_B2", "Surviving candidate", "shift-symmetric normal form"
    if signature == ("DtS", "DtS"):
        return "dotS2", "Surviving candidate", "shift-symmetric normal form"
    if signature == ("DtS", "E", "E"):
        return "DtS_E2", "Surviving candidate", "shift-symmetric normal form"
    if signature == ("DtS", "a", "a"):
        return "a2DtS", "Proven reducible", "lower-order EOM"
    if signature == ("E", "GradS", "a"):
        return "aEGradS", "Proven reducible", "lower-order EOM"
    if signature == ("GradE", "GradS"):
        return "divEGradS", "Surviving candidate", "survives current rules"
    if signature == ("GradS", "GradS"):
        return "gradS2", "Surviving candidate", "shift-symmetric normal form"
    raise ValueError(
        "Unhandled E/B+derivative-only-scalar signature "
        f"{signature} with representative {representative}"
    )


def enumerate_shift_scalar_sector_classes(
    max_weight: int = MAX_WEIGHT,
) -> tuple[ContractionClass, ...]:
    return enumerate_contraction_classes(
        max_weight=max_weight,
        block_types=EBDS_BLOCK_TYPES,
        classifier=classify_shift_scalar_contraction,
        require_even_parity=True,
    )


def shift_scalar_surviving_classes(
    max_weight: int = MAX_WEIGHT,
) -> tuple[ContractionClass, ...]:
    return tuple(
        item
        for item in enumerate_shift_scalar_sector_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    )


def shift_scalar_summary(max_weight: int = MAX_WEIGHT) -> ShiftScalarSummary:
    baseline_labels = set(eb_summary(max_weight=max_weight).surviving_labels)
    survivors = shift_scalar_surviving_classes(max_weight=max_weight)
    new_survivors = tuple(item.label for item in survivors if item.label not in baseline_labels)
    first_weight = min(
        (item.weight for item in survivors if item.label not in baseline_labels),
        default=None,
    )
    first_weight_labels = tuple(
        item.label
        for item in survivors
        if item.label not in baseline_labels and item.weight == first_weight
    )
    canonical = "dotS2" if "dotS2" in first_weight_labels else (
        None if not first_weight_labels else first_weight_labels[0]
    )
    return ShiftScalarSummary(
        total_classes=len(enumerate_shift_scalar_sector_classes(max_weight=max_weight)),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=new_survivors,
        first_new_weight=first_weight,
        first_new_labels=first_weight_labels,
        canonical_new_survivor=canonical,
    )


def shift_scalar_sector_report(max_weight: int = MAX_WEIGHT) -> str:
    classes = enumerate_shift_scalar_sector_classes(max_weight=max_weight)
    survivors = shift_scalar_surviving_classes(max_weight=max_weight)
    summary = shift_scalar_summary(max_weight=max_weight)

    lines = [
        "Delta<=4 E/B+derivative-only-scalar audit",
        "",
        f"Total parity-even scalar classes: {len(classes)}",
        f"Surviving classes under the current allowed rules: {len(survivors)}",
        "",
        "Surviving Delta<=4 scalar classes:",
    ]
    current_weight: int | None = None
    for item in survivors:
        if item.weight != current_weight:
            current_weight = item.weight
            lines.append(f"- weight {current_weight}:")
        lines.append(f"  {item.label}")

    lines.extend(
        [
            "",
            "Operational verdict:",
            f"- The first new survivors beyond the corrected E/B basis appear at weight {summary.first_new_weight}: {', '.join(summary.first_new_labels)}.",
            f"- A canonical derivative-only scalar obstruction is {summary.canonical_new_survivor}.",
            "- The derivative-only scalar extension still yields a finite corrected survivor list at Delta<=4.",
            "- Therefore shift symmetry removes the bare-scalar obstruction but does not rescue minimal-sector uniqueness by itself.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(shift_scalar_sector_report())


if __name__ == "__main__":
    main()
