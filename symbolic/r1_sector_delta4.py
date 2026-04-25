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


MAX_WEIGHT = 4

R1_BLOCK_TYPES = (
    BlockType("V", 1, 1, parity=0),
    BlockType("DtV", 2, 1, parity=0),
    BlockType("GradV", 2, 2, parity=1),
    BlockType("Dt2V", 3, 1, parity=0),
)

R1_FAMILY_NAMES = {block.name for block in R1_BLOCK_TYPES}


@dataclass(frozen=True)
class R1Summary:
    total_classes: int
    surviving_labels: tuple[str, ...]
    new_surviving_labels: tuple[str, ...]
    first_self_witness: str | None
    first_mixed_witness: str | None
    smallest_new_witness: str | None


def classify_r1_contraction(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[str, str, str]:
    if not any(name in R1_FAMILY_NAMES for name in signature):
        return classify_and_label(signature, representative)

    internal_pairs = {
        edge
        for edge in representative
        if edge[0][0] == edge[1][0]
    }
    components = connected_components(representative)

    if signature == ("V", "V"):
        return "V2", "Surviving candidate", "normal form"
    if signature == ("DtV", "V"):
        return "V_DtV", "Proven reducible", "total derivative"
    if signature == ("Dt2V", "V"):
        return "V_Dt2V", "Proven reducible", "total derivative to -dotV2"
    if signature == ("DtV", "DtV"):
        return "dotV2", "Surviving candidate", "normal form"
    if signature == ("E", "V", "V"):
        return "EVV", "Surviving candidate", "normal form"
    if signature == ("DtE", "V", "V"):
        return "DtEVV", "Proven reducible", "total derivative to -2 EVDtV"
    if signature == ("DtV", "E", "V"):
        return "EVDtV", "Surviving candidate", "survives current rules"
    if signature == ("V", "V", "V", "V"):
        return "V2^2", "Surviving candidate", "normal form"
    if signature == ("E", "E", "V", "V"):
        if len(components) == 2:
            return "E2V2", "Surviving candidate", "normal form"
        return "EVEV", "Surviving candidate", "survives current rules"
    if signature == ("V", "V", "a", "a"):
        if len(components) == 2:
            return "a2V2", "Proven reducible", "lower-order EOM"
        return "aVaV", "Proven reducible", "lower-order EOM"
    if signature == ("GradV", "V", "a"):
        if internal_pairs:
            return "aVGradV_1", "Proven reducible", "lower-order EOM"
        accel_edge = next(
            edge for edge in representative if edge[0][0] == 2 or edge[1][0] == 2
        )
        grad_endpoint = accel_edge[0] if accel_edge[0][0] == 0 else accel_edge[1]
        if grad_endpoint[1] == 0:
            return "aVGradV_2", "Proven reducible", "lower-order EOM"
        return "aVGradV_3", "Proven reducible", "lower-order EOM"
    if signature == ("GradV", "GradV"):
        if len(internal_pairs) == 2:
            return "divV2", "Surviving candidate", "survives current rules"
        if representative == (
            ((0, 0), (1, 0)),
            ((0, 1), (1, 1)),
        ):
            return "gradV2", "Surviving candidate", "normal form"
        return "mixedGradV2", "Surviving candidate", "survives current rules"
    raise ValueError(f"Unhandled R1 signature {signature} with representative {representative}")


def enumerate_r1_sector_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return enumerate_contraction_classes(
        max_weight=max_weight,
        block_types=BLOCK_TYPES + R1_BLOCK_TYPES,
        classifier=classify_r1_contraction,
        require_even_parity=True,
    )


def r1_surviving_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return tuple(
        item
        for item in enumerate_r1_sector_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    )


def _is_self_witness(item: ContractionClass) -> bool:
    return all(name in R1_FAMILY_NAMES for name in item.signature)


def r1_summary(max_weight: int = MAX_WEIGHT) -> R1Summary:
    baseline_labels = {
        item.label
        for item in enumerate_contraction_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    }
    survivors = r1_surviving_classes(max_weight=max_weight)
    new_survivors = tuple(item for item in survivors if item.label not in baseline_labels)
    self_candidates = tuple(item for item in new_survivors if _is_self_witness(item))
    mixed_candidates = tuple(item for item in new_survivors if not _is_self_witness(item))
    first_self = min(self_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_mixed = min(mixed_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_new = min(new_survivors, key=lambda item: (item.weight, item.label), default=None)
    return R1Summary(
        total_classes=len(enumerate_r1_sector_classes(max_weight=max_weight)),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=tuple(item.label for item in new_survivors),
        first_self_witness=None if first_self is None else first_self.label,
        first_mixed_witness=None if first_mixed is None else first_mixed.label,
        smallest_new_witness=None if first_new is None else first_new.label,
    )


def r1_sector_report(max_weight: int = MAX_WEIGHT) -> str:
    classes = enumerate_r1_sector_classes(max_weight=max_weight)
    survivors = r1_surviving_classes(max_weight=max_weight)
    summary = r1_summary(max_weight=max_weight)

    lines = [
        "Delta<=4 rank-1 vector-family audit",
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
            "- Derivative-generated vectors such as GradS or divergence blocks are not counted as primitive R1 admission in this audit.",
            "- The current audited-set thresholds for R2, R0a, and R0b already reduce the Delta<=4 baseline back to the electric sector, so the first mixed vector witness is read against that electric baseline.",
            f"- First self witness: {summary.first_self_witness}.",
            f"- First mixed witness: {summary.first_mixed_witness}.",
            f"- Smallest new witness beyond the audited-set baseline: {summary.smallest_new_witness}.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(r1_sector_report())


if __name__ == "__main__":
    main()
