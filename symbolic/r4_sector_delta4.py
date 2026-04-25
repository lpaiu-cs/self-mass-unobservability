from __future__ import annotations

from dataclasses import dataclass

from enumerate_contractions_delta4 import (
    BLOCK_TYPES,
    BlockType,
    ContractionClass,
    enumerate_contraction_classes,
)


MAX_WEIGHT = 4

R4_BLOCK_TYPES = (
    BlockType(
        "Q",
        1,
        4,
        parity=0,
        sym_groups=((0, 1), (1, 2), (2, 3)),
        tracefree_pairs=((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)),
    ),
    BlockType(
        "DtQ",
        2,
        4,
        parity=0,
        sym_groups=((0, 1), (1, 2), (2, 3)),
        tracefree_pairs=((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)),
    ),
    BlockType(
        "GradQ",
        2,
        5,
        parity=1,
        sym_groups=((1, 2), (2, 3), (3, 4)),
        tracefree_pairs=((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)),
    ),
    BlockType(
        "Dt2Q",
        3,
        4,
        parity=0,
        sym_groups=((0, 1), (1, 2), (2, 3)),
        tracefree_pairs=((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)),
    ),
)

R4_FAMILY_NAMES = {block.name for block in R4_BLOCK_TYPES}

REP_GRADQ_GRADQ_DIV = (
    ((0, 0), (0, 1)),
    ((0, 2), (1, 1)),
    ((0, 3), (1, 2)),
    ((0, 4), (1, 3)),
    ((1, 0), (1, 4)),
)
REP_GRADQ_GRADQ_GRAD = (
    ((0, 0), (1, 0)),
    ((0, 1), (1, 1)),
    ((0, 2), (1, 2)),
    ((0, 3), (1, 3)),
    ((0, 4), (1, 4)),
)
REP_GRADQ_GRADQ_MIXED = (
    ((0, 0), (1, 1)),
    ((0, 1), (1, 0)),
    ((0, 2), (1, 2)),
    ((0, 3), (1, 3)),
    ((0, 4), (1, 4)),
)
REP_QQQQ_PAIR = (
    ((0, 0), (3, 0)),
    ((0, 1), (3, 1)),
    ((0, 2), (3, 2)),
    ((0, 3), (3, 3)),
    ((1, 0), (2, 0)),
    ((1, 1), (2, 1)),
    ((1, 2), (2, 2)),
    ((1, 3), (2, 3)),
)
REP_QQQQ_CHAIN = (
    ((0, 0), (2, 0)),
    ((0, 1), (3, 0)),
    ((0, 2), (3, 1)),
    ((0, 3), (3, 2)),
    ((1, 0), (2, 1)),
    ((1, 1), (2, 2)),
    ((1, 2), (2, 3)),
    ((1, 3), (3, 3)),
)
REP_QQQQ_BRIDGE = (
    ((0, 0), (2, 0)),
    ((0, 1), (2, 1)),
    ((0, 2), (3, 0)),
    ((0, 3), (3, 1)),
    ((1, 0), (2, 2)),
    ((1, 1), (2, 3)),
    ((1, 2), (3, 2)),
    ((1, 3), (3, 3)),
)
REP_QQQQ_TETRA = (
    ((0, 0), (1, 0)),
    ((0, 1), (2, 0)),
    ((0, 2), (3, 0)),
    ((0, 3), (3, 1)),
    ((1, 1), (2, 1)),
    ((1, 2), (2, 2)),
    ((1, 3), (3, 2)),
    ((2, 3), (3, 3)),
)


@dataclass(frozen=True)
class R4Summary:
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


def _r4_family_classes() -> tuple[ContractionClass, ...]:
    classes = (
        _make_class(
            ("Q", "Q"),
            "Q2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
            ),
            2,
        ),
        _make_class(
            ("DtQ", "Q"),
            "Q_DtQ",
            "Proven reducible",
            "total derivative",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
            ),
            3,
        ),
        _make_class(
            ("E", "Q", "Q"),
            "EQQ",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
                ((1, 3), (2, 3)),
            ),
            3,
        ),
        _make_class(
            ("DtQ", "DtQ"),
            "dotQ2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
            ),
            4,
        ),
        _make_class(
            ("Dt2Q", "Q"),
            "Q_Dt2Q",
            "Proven reducible",
            "total derivative to -dotQ2",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
            ),
            4,
        ),
        _make_class(
            ("DtE", "Q", "Q"),
            "DtEQQ",
            "Proven reducible",
            "total derivative to -2 EQDtQ",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
                ((1, 3), (2, 3)),
            ),
            4,
        ),
        _make_class(
            ("DtQ", "E", "Q"),
            "EQDtQ",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((0, 2), (2, 1)),
                ((0, 3), (2, 2)),
                ((1, 1), (2, 3)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "Q", "Q"),
            "E2Q2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((2, 0), (3, 0)),
                ((2, 1), (3, 1)),
                ((2, 2), (3, 2)),
                ((2, 3), (3, 3)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "Q", "Q"),
            "E2Q2_mixed_1",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (3, 0)),
                ((2, 1), (3, 1)),
                ((2, 2), (3, 2)),
                ((2, 3), (3, 3)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "Q", "Q"),
            "E2Q2_mixed_2",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (2, 0)),
                ((0, 1), (2, 1)),
                ((1, 0), (3, 0)),
                ((1, 1), (3, 1)),
                ((2, 2), (3, 2)),
                ((2, 3), (3, 3)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "Q", "Q"),
            "E2Q2_mixed_3",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (2, 0)),
                ((0, 1), (3, 0)),
                ((1, 0), (2, 1)),
                ((1, 1), (3, 1)),
                ((2, 2), (3, 2)),
                ((2, 3), (3, 3)),
            ),
            4,
        ),
        _make_class(
            ("Q", "Q", "a", "a"),
            "a2Q2",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((2, 0), (3, 0)),
            ),
            4,
        ),
        _make_class(
            ("Q", "Q", "a", "a"),
            "aQaQ",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (2, 0)),
                ((1, 3), (3, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradQ", "Q", "a"),
            "aQGradQ_1",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (0, 1)),
                ((0, 2), (1, 0)),
                ((0, 3), (1, 1)),
                ((0, 4), (1, 2)),
                ((1, 3), (2, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradQ", "Q", "a"),
            "aQGradQ_2",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (2, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradQ", "Q", "a"),
            "aQGradQ_3",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (2, 0)),
                ((0, 1), (1, 0)),
                ((0, 2), (1, 1)),
                ((0, 3), (1, 2)),
                ((0, 4), (1, 3)),
            ),
            4,
        ),
        _make_class(
            ("GradQ", "GradQ"),
            "divQ2",
            "Surviving candidate",
            "survives current rules",
            REP_GRADQ_GRADQ_DIV,
            4,
        ),
        _make_class(
            ("GradQ", "GradQ"),
            "gradQ2",
            "Surviving candidate",
            "normal form",
            REP_GRADQ_GRADQ_GRAD,
            4,
        ),
        _make_class(
            ("GradQ", "GradQ"),
            "mixedGradQ2",
            "Surviving candidate",
            "survives current rules",
            REP_GRADQ_GRADQ_MIXED,
            4,
        ),
        _make_class(
            ("Q", "Q", "Q", "Q"),
            "Q2^2",
            "Surviving candidate",
            "normal form",
            REP_QQQQ_PAIR,
            4,
        ),
        _make_class(
            ("Q", "Q", "Q", "Q"),
            "Q4_chain",
            "Surviving candidate",
            "survives current rules",
            REP_QQQQ_CHAIN,
            4,
        ),
        _make_class(
            ("Q", "Q", "Q", "Q"),
            "Q4_bridge",
            "Surviving candidate",
            "survives current rules",
            REP_QQQQ_BRIDGE,
            4,
        ),
        _make_class(
            ("Q", "Q", "Q", "Q"),
            "Q4_tetra",
            "Surviving candidate",
            "survives current rules",
            REP_QQQQ_TETRA,
            4,
        ),
    )
    return tuple(sorted(classes, key=lambda item: (item.weight, item.signature, item.label)))


def enumerate_r4_family_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    if max_weight != MAX_WEIGHT:
        raise ValueError("The current rank-4 family audit is only fixed at Delta <= 4.")
    return _r4_family_classes()


def enumerate_r4_sector_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    baseline = enumerate_contraction_classes(max_weight=max_weight)
    family = enumerate_r4_family_classes(max_weight=max_weight)
    return tuple(sorted(baseline + family, key=lambda item: (item.weight, item.signature, item.label)))


def r4_surviving_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return tuple(
        item
        for item in enumerate_r4_sector_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    )


def _is_self_witness(item: ContractionClass) -> bool:
    return all(name in R4_FAMILY_NAMES for name in item.signature)


def r4_summary(max_weight: int = MAX_WEIGHT) -> R4Summary:
    baseline_labels = {
        item.label
        for item in enumerate_contraction_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    }
    survivors = r4_surviving_classes(max_weight=max_weight)
    new_survivors = tuple(item for item in survivors if item.label not in baseline_labels)
    self_candidates = tuple(item for item in new_survivors if _is_self_witness(item))
    mixed_candidates = tuple(item for item in new_survivors if not _is_self_witness(item))
    first_self = min(self_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_mixed = min(mixed_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_new = min(new_survivors, key=lambda item: (item.weight, item.label), default=None)
    return R4Summary(
        total_classes=len(enumerate_r4_sector_classes(max_weight=max_weight)),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=tuple(item.label for item in new_survivors),
        first_self_witness=None if first_self is None else first_self.label,
        first_mixed_witness=None if first_mixed is None else first_mixed.label,
        smallest_new_witness=None if first_new is None else first_new.label,
    )


def r4_sector_report(max_weight: int = MAX_WEIGHT) -> str:
    classes = enumerate_r4_sector_classes(max_weight=max_weight)
    survivors = r4_surviving_classes(max_weight=max_weight)
    summary = r4_summary(max_weight=max_weight)

    lines = [
        "Delta<=4 rank-4 family audit",
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
            "- The primitive family audited here is a genuine parity-even fully symmetric trace-free rank-4 tensor family Q_ijkl.",
            "- Trace descendants reducible to scalar or rank-2 audited families, and derivative-generated rank-4 blocks such as GradT, GradGradE, GradGradV, or scalar-Hessian descendants, are excluded from the primitive-family definition and are not double-counted here.",
            "- The current audited-set thresholds for R2, R0a, R0b, R1, and Rodd+ already reduce the Delta<=4 baseline back to the electric sector, so the first mixed rank-4 witness is read against that electric baseline.",
            f"- First self witness: {summary.first_self_witness}.",
            f"- First mixed witness: {summary.first_mixed_witness}.",
            f"- Smallest new witness beyond the enlarged audited-set baseline: {summary.smallest_new_witness}.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(r4_sector_report())


if __name__ == "__main__":
    main()
