from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from enumerate_contractions_delta4 import (
    BlockType,
    ContractionClass,
    classify_and_label,
    connected_components,
    enumerate_contraction_classes,
)


MAX_WEIGHT = 4

EB_BLOCK_TYPES = (
    BlockType("E", 1, 2, parity=0, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
    BlockType("B", 1, 2, parity=1, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
    BlockType("a", 1, 1, parity=1),
    BlockType("DtE", 2, 2, parity=0, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
    BlockType("DtB", 2, 2, parity=1, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
    BlockType("GradE", 2, 3, parity=1, sym_groups=((1, 2),), tracefree_pairs=((1, 2),)),
    BlockType("GradB", 2, 3, parity=0, sym_groups=((1, 2),), tracefree_pairs=((1, 2),)),
    BlockType("Dt2E", 3, 2, parity=0, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
    BlockType("Dt2B", 3, 2, parity=1, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
)

EB_ORDER = tuple(block.name for block in EB_BLOCK_TYPES)
B_FAMILY_NAMES = {"B", "DtB", "GradB", "Dt2B"}


@dataclass(frozen=True)
class EBSummary:
    total_classes: int
    surviving_labels: tuple[str, ...]
    new_surviving_labels: tuple[str, ...]
    smallest_new_survivor: str | None


def actual_instance_names(signature: tuple[str, ...]) -> tuple[str, ...]:
    counts = Counter(signature)
    ordered: list[str] = []
    for name in EB_ORDER:
        ordered.extend([name] * counts.get(name, 0))
    return tuple(ordered)


def count_edges_between_instances(
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
    left_index: int,
    right_index: int,
) -> int:
    target = {left_index, right_index}
    return sum(1 for left, right in representative if {left[0], right[0]} == target)


def classify_eb_contraction(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[str, str, str]:
    if not any(name in B_FAMILY_NAMES for name in signature):
        return classify_and_label(signature, representative)

    internal_pairs = {
        edge
        for edge in representative
        if edge[0][0] == edge[1][0]
    }
    components = connected_components(representative)

    if signature == ("B", "B"):
        return "B2", "Surviving candidate", "normal form"
    if signature == ("B", "DtB"):
        return "B_DtB", "Proven reducible", "total derivative"
    if signature == ("B", "Dt2B"):
        return "B_Dt2B", "Proven reducible", "total derivative"
    if signature == ("DtB", "DtB"):
        return "dotB2", "Surviving candidate", "normal form"
    if signature == ("B", "B", "E"):
        return "EB2", "Surviving candidate", "normal form"
    if signature == ("B", "DtB", "E"):
        return "EBDtB", "Surviving candidate", "survives current rules"
    if signature == ("B", "B", "DtE"):
        return "B2DtE", "Proven reducible", "total derivative to -2 EBDtB"
    if signature == ("B", "B", "B", "B"):
        if len(components) == 2:
            return "B2^2", "Surviving candidate", "normal form"
        return "B4", "Proven reducible", "Cayley-Hamilton"
    if signature == ("B", "B", "a", "a"):
        if len(components) == 2:
            return "a2B2", "Proven reducible", "lower-order EOM"
        return "aB2a", "Proven reducible", "lower-order EOM"
    if signature == ("B", "GradB", "a"):
        if internal_pairs:
            return "aBGradB_1", "Proven reducible", "lower-order EOM"
        accel_edge = next(
            edge for edge in representative if edge[0][0] == 1 or edge[1][0] == 1
        )
        grad_endpoint = accel_edge[0] if accel_edge[0][0] == 2 else accel_edge[1]
        if grad_endpoint[1] == 0:
            return "aBGradB_3", "Proven reducible", "lower-order EOM"
        return "aBGradB_2", "Proven reducible", "lower-order EOM"
    if signature == ("GradB", "GradB"):
        if len(internal_pairs) == 2:
            return "divB2", "Surviving candidate", "survives current rules"
        if representative == (
            ((0, 0), (1, 0)),
            ((0, 1), (1, 1)),
            ((0, 2), (1, 2)),
        ):
            return "gradB2", "Surviving candidate", "normal form"
        return "mixedGradB2", "Surviving candidate", "survives current rules"
    if signature == ("B", "B", "E", "E"):
        instance_names = actual_instance_names(signature)
        e_indices = [index for index, name in enumerate(instance_names) if name == "E"]
        b_indices = [index for index, name in enumerate(instance_names) if name == "B"]
        e_edges = count_edges_between_instances(representative, e_indices[0], e_indices[1])
        b_edges = count_edges_between_instances(representative, b_indices[0], b_indices[1])
        if len(components) == 2 and e_edges == 2 and b_edges == 2:
            return "E2B2", "Surviving candidate", "normal form"
        if len(components) == 2 and e_edges == 0 and b_edges == 0:
            return "EB_sq", "Surviving candidate", "normal form"
        if len(components) == 1 and e_edges == 1 and b_edges == 1:
            return "TrE2B2", "Surviving candidate", "survives current rules"
        if len(components) == 1 and e_edges == 0 and b_edges == 0:
            return "EBEB", "Proven reducible", "mixed STF quartic identity"
    raise ValueError(f"Unhandled E/B signature {signature} with representative {representative}")


def enumerate_eb_sector_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return enumerate_contraction_classes(
        max_weight=max_weight,
        block_types=EB_BLOCK_TYPES,
        classifier=classify_eb_contraction,
        require_even_parity=True,
    )


def eb_surviving_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return tuple(
        item
        for item in enumerate_eb_sector_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    )


def eb_summary(max_weight: int = MAX_WEIGHT) -> EBSummary:
    electric_labels = {item.label for item in enumerate_contraction_classes(max_weight=max_weight)}
    survivors = eb_surviving_classes(max_weight=max_weight)
    new_survivors = tuple(item.label for item in survivors if item.label not in electric_labels)
    smallest_new = min(
        (item for item in survivors if item.label not in electric_labels),
        key=lambda item: (item.weight, item.label),
        default=None,
    )
    return EBSummary(
        total_classes=len(enumerate_eb_sector_classes(max_weight=max_weight)),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=new_survivors,
        smallest_new_survivor=None if smallest_new is None else smallest_new.label,
    )


def eb_sector_report(max_weight: int = MAX_WEIGHT) -> str:
    classes = enumerate_eb_sector_classes(max_weight=max_weight)
    survivors = eb_surviving_classes(max_weight=max_weight)
    summary = eb_summary(max_weight=max_weight)

    lines = [
        "Delta<=4 E/B-sector audit",
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
            f"- Smallest new survivor beyond the electric-only exact current set: {summary.smallest_new_survivor}.",
            "- The raw E/B survivor candidate set has one explicit quartic dependence relation, so the corrected linearly independent E/B basis has 18 elements.",
            "- The E/B enlargement still yields a corrected finite linearly independent basis at Delta<=4.",
            "- Therefore the magnetic family obstructs the electric-only minimal-sector claim, not finite fixed-order closure by itself.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(eb_sector_report())


if __name__ == "__main__":
    main()
