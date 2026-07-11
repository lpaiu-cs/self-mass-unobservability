from __future__ import annotations

from dataclasses import dataclass

from enumerate_contractions_delta4 import (
    BLOCK_TYPES,
    BlockType,
    ContractionClass,
    connected_components,
    enumerate_contraction_classes,
)


MAX_WEIGHT = 4

R3_BLOCK_TYPES = (
    BlockType(
        "T",
        1,
        3,
        parity=0,
        sym_groups=((0, 1), (1, 2)),
        tracefree_pairs=((0, 1), (0, 2), (1, 2)),
    ),
    BlockType(
        "DtT",
        2,
        3,
        parity=0,
        sym_groups=((0, 1), (1, 2)),
        tracefree_pairs=((0, 1), (0, 2), (1, 2)),
    ),
    BlockType(
        "GradT",
        2,
        4,
        parity=1,
        sym_groups=((1, 2), (2, 3)),
        tracefree_pairs=((1, 2), (1, 3), (2, 3)),
    ),
    BlockType(
        "Dt2T",
        3,
        3,
        parity=0,
        sym_groups=((0, 1), (1, 2)),
        tracefree_pairs=((0, 1), (0, 2), (1, 2)),
    ),
)

R3_FAMILY_NAMES = {block.name for block in R3_BLOCK_TYPES}

REP_GRADT_GRADT_GRAD = (
    ((0, 0), (1, 0)),
    ((0, 1), (1, 1)),
    ((0, 2), (1, 2)),
    ((0, 3), (1, 3)),
)
REP_GRADT_GRADT_MIXED = (
    ((0, 0), (1, 1)),
    ((0, 1), (1, 0)),
    ((0, 2), (1, 2)),
    ((0, 3), (1, 3)),
)
REP_TTTT_TETRA = (
    ((0, 0), (1, 0)),
    ((0, 1), (2, 0)),
    ((0, 2), (3, 0)),
    ((1, 1), (2, 1)),
    ((1, 2), (3, 1)),
    ((2, 2), (3, 2)),
)


@dataclass(frozen=True)
class R3Summary:
    total_classes: int
    surviving_labels: tuple[str, ...]
    new_surviving_labels: tuple[str, ...]
    first_self_witness: str | None
    first_mixed_witness: str | None
    smallest_new_witness: str | None


def _make_class(
    signature: tuple[str, ...],
    label: str,
    classification: str,
    reduction_channel: str,
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
    weight: int,
) -> ContractionClass:
    return ContractionClass(
        signature=signature,
        label=label,
        classification=classification,
        status="Proven" if classification == "Proven reducible" else "Conjectural",
        weight=weight,
        reduction_channel=reduction_channel,
        representative=representative,
    )


def _signature_weight(signature: tuple[str, ...]) -> int:
    weights = {block.name: block.weight for block in BLOCK_TYPES + R3_BLOCK_TYPES}
    return sum(weights[name] for name in signature)


def _r3_family_classes() -> tuple[ContractionClass, ...]:
    classes = (
        _make_class(
            ("T", "T"),
            "T2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
            ),
            2,
        ),
        _make_class(
            ("DtT", "T"),
            "T_DtT",
            "Proven reducible",
            "total derivative",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
            ),
            3,
        ),
        _make_class(
            ("E", "T", "T"),
            "ETT",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
            ),
            3,
        ),
        _make_class(
            ("DtT", "DtT"),
            "dotT2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
            ),
            4,
        ),
        _make_class(
            ("Dt2T", "T"),
            "T_Dt2T",
            "Proven reducible",
            "total derivative to -dotT2",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
            ),
            4,
        ),
        _make_class(
            ("DtE", "T", "T"),
            "DtETT",
            "Proven reducible",
            "total derivative to -2 ETDtT",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
            ),
            4,
        ),
        _make_class(
            ("DtT", "E", "T"),
            "ETDtT",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((0, 2), (2, 1)),
                ((1, 1), (2, 2)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "T", "T"),
            "E2T2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((2, 0), (3, 0)),
                ((2, 1), (3, 1)),
                ((2, 2), (3, 2)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "T", "T"),
            "E2T2_mixed_1",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (3, 0)),
                ((2, 1), (3, 1)),
                ((2, 2), (3, 2)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "T", "T"),
            "E2T2_mixed_2",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (2, 0)),
                ((0, 1), (2, 1)),
                ((1, 0), (3, 0)),
                ((1, 1), (3, 1)),
                ((2, 2), (3, 2)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "T", "T"),
            "E2T2_mixed_3",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (2, 0)),
                ((0, 1), (3, 0)),
                ((1, 0), (2, 1)),
                ((1, 1), (3, 1)),
                ((2, 2), (3, 2)),
            ),
            4,
        ),
        _make_class(
            ("T", "T", "a", "a"),
            "a2T2",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((2, 0), (3, 0)),
            ),
            4,
        ),
        _make_class(
            ("T", "T", "a", "a"),
            "aTaT",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (2, 0)),
                ((1, 2), (3, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradT", "T", "a"),
            "aTGradT_1",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (0, 1)),
                ((0, 2), (1, 0)),
                ((0, 3), (1, 1)),
                ((1, 2), (2, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradT", "T", "a"),
            "aTGradT_2",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (2, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradT", "T", "a"),
            "aTGradT_3",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (2, 0)),
                ((0, 1), (1, 0)),
                ((0, 2), (1, 1)),
                ((0, 3), (1, 2)),
            ),
            4,
        ),
        _make_class(
            ("GradT", "GradT"),
            "divT2",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (0, 1)),
                ((0, 2), (1, 1)),
                ((0, 3), (1, 2)),
                ((1, 0), (1, 3)),
            ),
            4,
        ),
        _make_class(
            ("GradT", "GradT"),
            "gradT2",
            "Surviving candidate",
            "normal form",
            REP_GRADT_GRADT_GRAD,
            4,
        ),
        _make_class(
            ("GradT", "GradT"),
            "mixedGradT2",
            "Surviving candidate",
            "survives current rules",
            REP_GRADT_GRADT_MIXED,
            4,
        ),
        _make_class(
            ("T", "T", "T", "T"),
            "T2^2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((2, 0), (3, 0)),
                ((2, 1), (3, 1)),
                ((2, 2), (3, 2)),
            ),
            4,
        ),
        _make_class(
            ("T", "T", "T", "T"),
            "T4_chain",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (2, 0)),
                ((1, 2), (3, 0)),
                ((2, 1), (3, 1)),
                ((2, 2), (3, 2)),
            ),
            4,
        ),
        _make_class(
            ("T", "T", "T", "T"),
            "T4_tetra",
            "Surviving candidate",
            "survives current rules",
            REP_TTTT_TETRA,
            4,
        ),
    )
    return tuple(sorted(classes, key=lambda item: (item.weight, item.signature, item.label)))


def enumerate_r3_family_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    if max_weight != MAX_WEIGHT:
        raise ValueError("The current rank-3 family audit is only fixed at Delta <= 4.")
    return _r3_family_classes()


def enumerate_r3_sector_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    baseline = enumerate_contraction_classes(max_weight=max_weight)
    family = enumerate_r3_family_classes(max_weight=max_weight)
    return tuple(sorted(baseline + family, key=lambda item: (item.weight, item.signature, item.label)))


def r3_surviving_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return tuple(
        item
        for item in enumerate_r3_sector_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    )


def _is_self_witness(item: ContractionClass) -> bool:
    return all(name in R3_FAMILY_NAMES for name in item.signature)


def r3_summary(max_weight: int = MAX_WEIGHT) -> R3Summary:
    baseline_labels = {
        item.label
        for item in enumerate_contraction_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    }
    survivors = r3_surviving_classes(max_weight=max_weight)
    new_survivors = tuple(item for item in survivors if item.label not in baseline_labels)
    self_candidates = tuple(item for item in new_survivors if _is_self_witness(item))
    mixed_candidates = tuple(item for item in new_survivors if not _is_self_witness(item))
    first_self = min(self_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_mixed = min(mixed_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_new = min(new_survivors, key=lambda item: (item.weight, item.label), default=None)
    return R3Summary(
        total_classes=len(enumerate_r3_sector_classes(max_weight=max_weight)),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=tuple(item.label for item in new_survivors),
        first_self_witness=None if first_self is None else first_self.label,
        first_mixed_witness=None if first_mixed is None else first_mixed.label,
        smallest_new_witness=None if first_new is None else first_new.label,
    )


def r3_sector_report(max_weight: int = MAX_WEIGHT) -> str:
    classes = enumerate_r3_sector_classes(max_weight=max_weight)
    survivors = r3_surviving_classes(max_weight=max_weight)
    summary = r3_summary(max_weight=max_weight)

    lines = [
        "Delta<=4 rank-3 family audit",
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
            "- The primitive family audited here is a genuine parity-even fully symmetric trace-free rank-3 tensor family T_ijk.",
            "- Derivative-generated rank-3 blocks such as GradE, GradB, GradV, or scalar-derivative descendants are excluded from the primitive-family definition and are not double-counted here.",
            "- The current audited-set thresholds for R2, R0a, R0b, and R1 already reduce the Delta<=4 baseline back to the electric sector, so the first mixed rank-3 witness is read against that electric baseline.",
            f"- First self witness: {summary.first_self_witness}.",
            f"- First mixed witness: {summary.first_mixed_witness}.",
            f"- Smallest new witness beyond the enlarged audited-set baseline: {summary.smallest_new_witness}.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(r3_sector_report())


if __name__ == "__main__":
    main()
