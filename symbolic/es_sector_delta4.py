from __future__ import annotations

from dataclasses import dataclass

from enumerate_contractions_delta4 import BlockType, ContractionClass, enumerate_contraction_classes
from eb_sector_delta4 import EB_BLOCK_TYPES, classify_eb_contraction, eb_summary


MAX_WEIGHT = 4

SCALAR_BLOCK_TYPES = (
    BlockType("S", 1, 0, parity=0),
    BlockType("DtS", 2, 0, parity=0),
    BlockType("GradS", 2, 1, parity=1),
    BlockType("Dt2S", 3, 0, parity=0),
)

ES_BLOCK_TYPES = EB_BLOCK_TYPES + SCALAR_BLOCK_TYPES
SCALAR_FAMILY_NAMES = {block.name for block in SCALAR_BLOCK_TYPES}


@dataclass(frozen=True)
class ESSummary:
    total_classes: int
    surviving_labels: tuple[str, ...]
    new_surviving_labels: tuple[str, ...]
    smallest_new_survivor: str | None


def classify_es_contraction(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[str, str, str]:
    if not any(name in SCALAR_FAMILY_NAMES for name in signature):
        return classify_eb_contraction(signature, representative)

    if signature == ("S",):
        return "S", "Surviving candidate", "normal form"
    if signature == ("DtS",):
        return "DtS", "Proven reducible", "total derivative"
    if signature == ("S", "S"):
        return "S2", "Surviving candidate", "normal form"
    if signature == ("B", "B", "S"):
        return "SB2", "Surviving candidate", "normal form"
    if signature == ("Dt2S",):
        return "Dt2S", "Proven reducible", "total derivative"
    if signature == ("DtS", "S"):
        return "S_DtS", "Proven reducible", "total derivative"
    if signature == ("E", "E", "S"):
        return "SE2", "Surviving candidate", "normal form"
    if signature == ("GradS", "a"):
        return "aGradS", "Proven reducible", "lower-order EOM"
    if signature == ("S", "S", "S"):
        return "S3", "Surviving candidate", "normal form"
    if signature == ("S", "a", "a"):
        return "a2S", "Proven reducible", "lower-order EOM"
    if signature == ("B", "B", "DtS"):
        return "DtS_B2", "Surviving candidate", "normal form"
    if signature == ("B", "B", "E", "S"):
        return "SEB2", "Surviving candidate", "normal form"
    if signature == ("B", "B", "S", "S"):
        return "S2B2", "Surviving candidate", "normal form"
    if signature == ("B", "DtB", "S"):
        return "S_BDtB", "Proven reducible", "total derivative to -DtS_B2/2"
    if signature == ("Dt2S", "S"):
        return "S_Dt2S", "Proven reducible", "total derivative to -dotS2"
    if signature == ("DtE", "E", "S"):
        return "SE_DtE", "Proven reducible", "total derivative to -DtS_E2/2"
    if signature == ("DtS", "DtS"):
        return "dotS2", "Surviving candidate", "normal form"
    if signature == ("DtS", "E", "E"):
        return "DtS_E2", "Surviving candidate", "normal form"
    if signature == ("DtS", "S", "S"):
        return "S2DtS", "Proven reducible", "total derivative"
    if signature == ("DtS", "a", "a"):
        return "a2DtS", "Proven reducible", "lower-order EOM"
    if signature == ("E", "E", "E", "S"):
        return "SE3", "Surviving candidate", "normal form"
    if signature == ("E", "E", "S", "S"):
        return "S2E2", "Surviving candidate", "normal form"
    if signature == ("E", "GradS", "a"):
        return "aEGradS", "Proven reducible", "lower-order EOM"
    if signature == ("E", "S", "a", "a"):
        return "aSEa", "Proven reducible", "lower-order EOM"
    if signature == ("GradE", "GradS"):
        return "divEGradS", "Surviving candidate", "survives current rules"
    if signature == ("GradE", "S", "a"):
        return "SaDivE", "Proven reducible", "lower-order EOM"
    if signature == ("GradS", "GradS"):
        return "gradS2", "Surviving candidate", "normal form"
    if signature == ("GradS", "S", "a"):
        return "SaGradS", "Proven reducible", "lower-order EOM"
    if signature == ("S", "S", "S", "S"):
        return "S4", "Surviving candidate", "normal form"
    if signature == ("S", "S", "a", "a"):
        return "a2S2", "Proven reducible", "lower-order EOM"
    raise ValueError(
        f"Unhandled E/B+scalar signature {signature} with representative {representative}"
    )


def enumerate_es_sector_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return enumerate_contraction_classes(
        max_weight=max_weight,
        block_types=ES_BLOCK_TYPES,
        classifier=classify_es_contraction,
        require_even_parity=True,
    )


def es_surviving_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return tuple(
        item
        for item in enumerate_es_sector_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    )


def es_summary(max_weight: int = MAX_WEIGHT) -> ESSummary:
    baseline_labels = set(eb_summary(max_weight=max_weight).surviving_labels)
    survivors = es_surviving_classes(max_weight=max_weight)
    new_survivors = tuple(item.label for item in survivors if item.label not in baseline_labels)
    smallest_new = min(
        (item for item in survivors if item.label not in baseline_labels),
        key=lambda item: (item.weight, item.label),
        default=None,
    )
    return ESSummary(
        total_classes=len(enumerate_es_sector_classes(max_weight=max_weight)),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=new_survivors,
        smallest_new_survivor=None if smallest_new is None else smallest_new.label,
    )


def es_sector_report(max_weight: int = MAX_WEIGHT) -> str:
    classes = enumerate_es_sector_classes(max_weight=max_weight)
    survivors = es_surviving_classes(max_weight=max_weight)
    summary = es_summary(max_weight=max_weight)

    lines = [
        "Delta<=4 E/B+scalar-sector audit",
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
            f"- Smallest new survivor beyond the corrected E/B basis: {summary.smallest_new_survivor}.",
            "- The scalar-family extension still yields a finite corrected candidate basis at Delta<=4.",
            "- Therefore the scalar-like family obstructs adequacy claims for the current E/B exact-current-set theorem, not fixed-order finiteness by itself.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(es_sector_report())


if __name__ == "__main__":
    main()
