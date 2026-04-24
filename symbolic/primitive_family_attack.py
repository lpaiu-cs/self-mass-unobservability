from __future__ import annotations

from dataclasses import dataclass

from enumerate_contractions_delta4 import (
    BLOCK_TYPES,
    BlockType,
    ContractionClass,
    classify_and_label,
    connected_components,
    enumerate_contraction_classes,
)


MAGNETIC_BLOCK = BlockType(
    "B",
    1,
    2,
    parity=1,
    sym_groups=((0, 1),),
    tracefree_pairs=((0, 1),),
)


@dataclass(frozen=True)
class PrimitiveAttackSummary:
    extension_name: str
    total_classes: int
    new_class_count: int
    new_survivor_labels: tuple[str, ...]
    smallest_new_survivor: str | None


def classify_with_magnetic_family(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[str, str, str]:
    if "B" not in signature:
        return classify_and_label(signature, representative)

    components = connected_components(representative)
    if signature == ("B", "B"):
        return "B2", "Surviving candidate", "new primitive family"
    if signature == ("B", "B", "E"):
        return "EB2", "Surviving candidate", "new primitive family"
    if signature == ("B", "B", "DtE"):
        return "B2DtE", "Surviving candidate", "new primitive family"
    if signature == ("B", "B", "a", "a"):
        if len(components) == 2:
            return "a2B2", "Proven reducible", "lower-order EOM"
        return "aB2a", "Proven reducible", "lower-order EOM"
    if signature == ("B", "B", "B", "B"):
        if len(components) == 2:
            return "B2^2", "Surviving candidate", "new primitive family"
        return "B4", "Proven reducible", "Cayley-Hamilton"
    if signature == ("B", "B", "E", "E"):
        bb_edges = sum(1 for left, right in representative if {left[0], right[0]} == {0, 1})
        ee_edges = sum(1 for left, right in representative if {left[0], right[0]} == {2, 3})
        if len(components) == 2 and bb_edges == 2 and ee_edges == 2:
            return "E2B2", "Surviving candidate", "new primitive family"
        if len(components) == 2 and bb_edges == 0 and ee_edges == 0:
            return "EB_sq", "Surviving candidate", "new primitive family"
        if len(components) == 1 and bb_edges == 1 and ee_edges == 1:
            return "TrE2B2", "Surviving candidate", "new primitive family"
        if len(components) == 1 and bb_edges == 0 and ee_edges == 0:
            return "EBEB", "Surviving candidate", "new primitive family"
    raise ValueError(f"Unhandled magnetic-family signature {signature} with representative {representative}")


def enumerate_with_magnetic_family(max_weight: int = 4) -> tuple[ContractionClass, ...]:
    return enumerate_contraction_classes(
        max_weight=max_weight,
        block_types=BLOCK_TYPES + (MAGNETIC_BLOCK,),
        classifier=classify_with_magnetic_family,
        require_even_parity=True,
    )


def primitive_attack_summary(max_weight: int = 4) -> PrimitiveAttackSummary:
    base_labels = {item.label for item in enumerate_contraction_classes(max_weight=max_weight)}
    extended = enumerate_with_magnetic_family(max_weight=max_weight)
    new_classes = tuple(item for item in extended if item.label not in base_labels)
    new_survivors = tuple(
        item for item in new_classes if item.classification == "Surviving candidate"
    )
    smallest = min(new_survivors, key=lambda item: (item.weight, item.label), default=None)
    return PrimitiveAttackSummary(
        extension_name="magnetic tidal STF family B_ij",
        total_classes=len(extended),
        new_class_count=len(new_classes),
        new_survivor_labels=tuple(item.label for item in new_survivors),
        smallest_new_survivor=None if smallest is None else smallest.label,
    )


def primitive_attack_report(max_weight: int = 4) -> str:
    summary = primitive_attack_summary(max_weight=max_weight)
    extended = enumerate_with_magnetic_family(max_weight=max_weight)
    base_labels = {item.label for item in enumerate_contraction_classes(max_weight=max_weight)}
    new_classes = [item for item in extended if item.label not in base_labels]

    lines = [
        "Primitive-family adequacy attack",
        "",
        f"Extension family: {summary.extension_name}",
        f"Total Delta<=4 classes after extension: {summary.total_classes}",
        f"New classes relative to the exact current primitive set: {summary.new_class_count}",
        "",
        "New Delta<=4 classes:",
    ]
    for item in new_classes:
        lines.append(
            f"- {item.label} [weight {item.weight}; {item.classification}; {item.reduction_channel}]"
        )

    lines.extend(
        [
            "",
            "Adequacy result:",
            "- A new surviving scalar appears as soon as the magnetic tidal family is admitted.",
            f"- Smallest explicit new survivor: {summary.smallest_new_survivor}.",
            "- Therefore the exact-current-set theorem does not yet upgrade to a physically justified minimal sector theorem.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(primitive_attack_report())


if __name__ == "__main__":
    main()
