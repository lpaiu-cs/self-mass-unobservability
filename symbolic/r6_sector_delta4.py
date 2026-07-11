from __future__ import annotations

from dataclasses import dataclass

from enumerate_contractions_delta4 import ContractionClass, enumerate_contraction_classes
from high_rank_family_enumerator import GeneratedHighRankClass, generated_family_classes


MAX_WEIGHT = 4
R6_FAMILY_NAMES = {"Z", "DtZ", "GradZ"}
R6_GRADIENT_LABELS = {
    (1, 1, 0, 0, 0, 5): "divZ2",
    (0, 0, 1, 0, 0, 6): "gradZ2",
    (0, 0, 0, 1, 1, 5): "mixedGradZ2",
}
R6_DIRECT_E2Z2_KEY = (2, 0, 0, 0, 0, 6)


@dataclass(frozen=True)
class R6Summary:
    total_classes: int
    total_new_family_classes: int
    surviving_labels: tuple[str, ...]
    new_surviving_labels: tuple[str, ...]
    first_self_witness: str | None
    first_mixed_witness: str | None
    first_mixed_layer_labels: tuple[str, ...]
    smallest_new_witness: str | None
    even_rank_additional_first_mixed_labels: tuple[str, ...]


def _normalized_r6_label(
    item: GeneratedHighRankClass,
    signature_index_maps: dict[tuple[str, ...], dict[tuple, str]],
) -> str:
    signature = item.signature
    if signature == ("GradZ", "GradZ"):
        return R6_GRADIENT_LABELS[item.key[2]]
    if signature == ("E", "E", "Z", "Z"):
        if item.key[2] == R6_DIRECT_E2Z2_KEY:
            return "E2Z2"
        return signature_index_maps[signature][item.key]
    if signature == ("Z", "Z", "Z", "Z"):
        return signature_index_maps[signature][item.key]
    if signature == ("E", "Z", "Z", "Z"):
        return signature_index_maps[signature][item.key]
    return item.label


def _signature_index_maps(items: tuple[GeneratedHighRankClass, ...]) -> dict[tuple[str, ...], dict[tuple, str]]:
    maps: dict[tuple[str, ...], dict[tuple, str]] = {}
    for signature in {
        ("E", "E", "Z", "Z"),
        ("E", "Z", "Z", "Z"),
        ("Z", "Z", "Z", "Z"),
    }:
        same_signature = [item for item in items if item.signature == signature]
        if not same_signature:
            continue
        key_to_label: dict[tuple, str] = {}
        next_index = 1
        for item in sorted(same_signature, key=lambda candidate: candidate.key):
            if signature == ("E", "E", "Z", "Z") and item.key[2] == R6_DIRECT_E2Z2_KEY:
                key_to_label[item.key] = "E2Z2"
                continue
            if signature == ("E", "E", "Z", "Z"):
                key_to_label[item.key] = f"E2Z2_mixed_{next_index}"
            elif signature == ("E", "Z", "Z", "Z"):
                key_to_label[item.key] = f"EZZZ_{next_index}"
            else:
                key_to_label[item.key] = f"Z4_{next_index}"
            next_index += 1
        maps[signature] = key_to_label
    return maps


def _generated_to_contraction(
    item: GeneratedHighRankClass,
    signature_index_maps: dict[tuple[str, ...], dict[tuple, str]],
) -> ContractionClass:
    return ContractionClass(
        signature=item.signature,
        label=_normalized_r6_label(item, signature_index_maps),
        classification="Surviving candidate",
        status="Conjectural",
        weight=item.weight,
        reduction_channel="exhaustive normal-form candidate generation",
        representative=item.representative,
    )


def enumerate_r6_family_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    if max_weight != MAX_WEIGHT:
        raise ValueError("The current rank-6 family audit is only fixed at Delta <= 4.")
    generated = generated_family_classes(6)
    signature_index_maps = _signature_index_maps(generated)
    classes = tuple(_generated_to_contraction(item, signature_index_maps) for item in generated)
    return tuple(sorted(classes, key=lambda item: (item.weight, item.signature, item.label)))


def enumerate_r6_sector_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    baseline = enumerate_contraction_classes(max_weight=max_weight)
    family = enumerate_r6_family_classes(max_weight=max_weight)
    return tuple(sorted(baseline + family, key=lambda item: (item.weight, item.signature, item.label)))


def r6_surviving_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return tuple(
        item
        for item in enumerate_r6_sector_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    )


def _is_self_witness(item: ContractionClass) -> bool:
    return all(name in R6_FAMILY_NAMES for name in item.signature)


def r6_summary(max_weight: int = MAX_WEIGHT) -> R6Summary:
    baseline_labels = {
        item.label
        for item in enumerate_contraction_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    }
    family_classes = enumerate_r6_family_classes(max_weight=max_weight)
    survivors = r6_surviving_classes(max_weight=max_weight)
    new_survivors = tuple(item for item in survivors if item.label not in baseline_labels)
    self_candidates = tuple(item for item in new_survivors if _is_self_witness(item))
    mixed_candidates = tuple(item for item in new_survivors if not _is_self_witness(item))
    first_self = min(self_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_mixed = min(mixed_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_new = min(new_survivors, key=lambda item: (item.weight, item.label), default=None)
    first_mixed_weight = None if first_mixed is None else first_mixed.weight
    first_mixed_layer_labels = tuple(
        item.label
        for item in mixed_candidates
        if first_mixed_weight is not None and item.weight == first_mixed_weight
    )
    even_rank_additional_first_mixed = tuple(
        label for label in first_mixed_layer_labels if label != (first_mixed.label if first_mixed is not None else "")
    )
    return R6Summary(
        total_classes=len(enumerate_r6_sector_classes(max_weight=max_weight)),
        total_new_family_classes=len(family_classes),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=tuple(item.label for item in new_survivors),
        first_self_witness=None if first_self is None else first_self.label,
        first_mixed_witness=None if first_mixed is None else first_mixed.label,
        first_mixed_layer_labels=first_mixed_layer_labels,
        smallest_new_witness=None if first_new is None else first_new.label,
        even_rank_additional_first_mixed_labels=even_rank_additional_first_mixed,
    )


def r6_sector_report(max_weight: int = MAX_WEIGHT) -> str:
    summary = r6_summary(max_weight=max_weight)
    lines = [
        "Delta<=4 rank-6 family audit",
        "",
        f"Total parity-even scalar classes in the exhaustive candidate layer: {summary.total_classes}",
        f"New rank-6 candidate-survivor labels beyond the electric baseline: {summary.total_new_family_classes}",
        "",
        "New rank-6 candidate-survivor labels:",
        "- " + ", ".join(summary.new_surviving_labels),
        "",
        "Operational verdict:",
        "- The primitive family audited here is a genuine parity-even fully symmetric trace-free rank-6 tensor family Z_ijklmn.",
        "- Trace descendants reducible to lower-rank audited families and derivative-generated rank-6 blocks such as GradU, second-gradient rank-4 descendants, higher-gradient vector descendants, or scalar descendants dressed only by derivatives are excluded from the primitive-family definition and are not double-counted here.",
        "- This script uses exhaustive high-rank candidate generation from the start rather than manual survivor bookkeeping.",
        f"- First self witness: {summary.first_self_witness}.",
        f"- First mixed witness: {summary.first_mixed_witness}.",
        f"- First mixed-layer labels: {', '.join(summary.first_mixed_layer_labels) if summary.first_mixed_layer_labels else 'none'}.",
        (
            "- Additional same-order even-rank mixed classes beyond the first mixed witness: "
            + (
                ", ".join(summary.even_rank_additional_first_mixed_labels)
                if summary.even_rank_additional_first_mixed_labels
                else "none"
            )
            + "."
        ),
        f"- Smallest new witness beyond the current enlarged audited-set baseline: {summary.smallest_new_witness}.",
    ]
    return "\n".join(lines)


def main() -> None:
    print(r6_sector_report())


if __name__ == "__main__":
    main()
