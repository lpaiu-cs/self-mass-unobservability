from __future__ import annotations

from dataclasses import dataclass

from enumerate_contractions_delta4 import (
    BlockType,
    ContractionClass,
    enumerate_contraction_classes,
)


MAX_WEIGHT = 4

R5_BLOCK_TYPES = (
    BlockType(
        "U",
        1,
        5,
        parity=0,
        sym_groups=((0, 1), (1, 2), (2, 3), (3, 4)),
        tracefree_pairs=(
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 3),
            (2, 4),
            (3, 4),
        ),
    ),
    BlockType(
        "DtU",
        2,
        5,
        parity=0,
        sym_groups=((0, 1), (1, 2), (2, 3), (3, 4)),
        tracefree_pairs=(
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 3),
            (2, 4),
            (3, 4),
        ),
    ),
    BlockType(
        "GradU",
        2,
        6,
        parity=1,
        sym_groups=((1, 2), (2, 3), (3, 4), (4, 5)),
        tracefree_pairs=(
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
            (2, 3),
            (2, 4),
            (2, 5),
            (3, 4),
            (3, 5),
            (4, 5),
        ),
    ),
    BlockType(
        "Dt2U",
        3,
        5,
        parity=0,
        sym_groups=((0, 1), (1, 2), (2, 3), (3, 4)),
        tracefree_pairs=(
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 3),
            (2, 4),
            (3, 4),
        ),
    ),
)

R5_FAMILY_NAMES = {block.name for block in R5_BLOCK_TYPES}

PAIR_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


@dataclass(frozen=True)
class R5Summary:
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


def _multigraph_representative(
    degrees: tuple[int, ...],
    multiplicities: tuple[int, int, int, int, int, int],
) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    slot_cursors = [0] * len(degrees)
    representative: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for (left, right), multiplicity in zip(PAIR_EDGES, multiplicities):
        for _ in range(multiplicity):
            representative.append(
                ((left, slot_cursors[left]), (right, slot_cursors[right]))
            )
            slot_cursors[left] += 1
            slot_cursors[right] += 1
    if tuple(slot_cursors) != degrees:
        raise ValueError(
            f"Representative construction drifted: got {tuple(slot_cursors)} expected {degrees}."
        )
    return tuple(representative)


REP_GRADU_GRADU_DIV = (
    ((0, 0), (0, 1)),
    ((0, 2), (1, 1)),
    ((0, 3), (1, 2)),
    ((0, 4), (1, 3)),
    ((0, 5), (1, 4)),
    ((1, 0), (1, 5)),
)
REP_GRADU_GRADU_GRAD = (
    ((0, 0), (1, 0)),
    ((0, 1), (1, 1)),
    ((0, 2), (1, 2)),
    ((0, 3), (1, 3)),
    ((0, 4), (1, 4)),
    ((0, 5), (1, 5)),
)
REP_GRADU_GRADU_MIXED = (
    ((0, 0), (1, 1)),
    ((0, 1), (1, 0)),
    ((0, 2), (1, 2)),
    ((0, 3), (1, 3)),
    ((0, 4), (1, 4)),
    ((0, 5), (1, 5)),
)

REP_E2U2 = _multigraph_representative((2, 2, 5, 5), (2, 0, 0, 0, 0, 5))
REP_E2U2_MIXED_1 = _multigraph_representative((2, 2, 5, 5), (1, 0, 1, 1, 0, 4))
REP_E2U2_MIXED_2 = _multigraph_representative((2, 2, 5, 5), (0, 0, 2, 2, 0, 3))
REP_E2U2_MIXED_3 = _multigraph_representative((2, 2, 5, 5), (0, 1, 1, 1, 1, 3))

REP_UUUU_PAIR = _multigraph_representative((5, 5, 5, 5), (0, 0, 5, 5, 0, 0))
REP_UUUU_CHAIN = _multigraph_representative((5, 5, 5, 5), (0, 1, 4, 4, 1, 0))
REP_UUUU_BRIDGE = _multigraph_representative((5, 5, 5, 5), (0, 2, 3, 3, 2, 0))
REP_UUUU_TETRA = _multigraph_representative((5, 5, 5, 5), (1, 1, 3, 3, 1, 1))
REP_UUUU_BALANCED = _multigraph_representative((5, 5, 5, 5), (1, 2, 2, 2, 2, 1))


def _r5_family_classes() -> tuple[ContractionClass, ...]:
    classes = (
        _make_class(
            ("U", "U"),
            "U2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
            ),
            2,
        ),
        _make_class(
            ("DtU", "U"),
            "U_DtU",
            "Proven reducible",
            "total derivative",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
            ),
            3,
        ),
        _make_class(
            ("E", "U", "U"),
            "EUU",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
                ((1, 3), (2, 3)),
                ((1, 4), (2, 4)),
            ),
            3,
        ),
        _make_class(
            ("DtU", "DtU"),
            "dotU2",
            "Surviving candidate",
            "normal form",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
            ),
            4,
        ),
        _make_class(
            ("Dt2U", "U"),
            "U_Dt2U",
            "Proven reducible",
            "total derivative to -dotU2",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
            ),
            4,
        ),
        _make_class(
            ("DtE", "U", "U"),
            "DtEUU",
            "Proven reducible",
            "total derivative to -2 EUDtU",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
                ((1, 3), (2, 3)),
                ((1, 4), (2, 4)),
            ),
            4,
        ),
        _make_class(
            ("DtU", "E", "U"),
            "EUDtU",
            "Surviving candidate",
            "survives current rules",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((0, 2), (2, 1)),
                ((0, 3), (2, 2)),
                ((0, 4), (2, 3)),
                ((1, 1), (2, 4)),
            ),
            4,
        ),
        _make_class(
            ("E", "E", "U", "U"),
            "E2U2",
            "Surviving candidate",
            "normal form",
            REP_E2U2,
            4,
        ),
        _make_class(
            ("E", "E", "U", "U"),
            "E2U2_mixed_1",
            "Surviving candidate",
            "survives current rules",
            REP_E2U2_MIXED_1,
            4,
        ),
        _make_class(
            ("E", "E", "U", "U"),
            "E2U2_mixed_2",
            "Surviving candidate",
            "survives current rules",
            REP_E2U2_MIXED_2,
            4,
        ),
        _make_class(
            ("E", "E", "U", "U"),
            "E2U2_mixed_3",
            "Surviving candidate",
            "survives current rules",
            REP_E2U2_MIXED_3,
            4,
        ),
        _make_class(
            ("U", "U", "a", "a"),
            "a2U2",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
                ((2, 0), (3, 0)),
            ),
            4,
        ),
        _make_class(
            ("U", "U", "a", "a"),
            "aUaU",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (2, 0)),
                ((1, 4), (3, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradU", "U", "a"),
            "aUGradU_1",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (0, 1)),
                ((0, 2), (1, 0)),
                ((0, 3), (1, 1)),
                ((0, 4), (1, 2)),
                ((0, 5), (1, 3)),
                ((1, 4), (2, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradU", "U", "a"),
            "aUGradU_2",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
                ((0, 5), (2, 0)),
            ),
            4,
        ),
        _make_class(
            ("GradU", "U", "a"),
            "aUGradU_3",
            "Proven reducible",
            "lower-order EOM",
            (
                ((0, 0), (2, 0)),
                ((0, 1), (1, 0)),
                ((0, 2), (1, 1)),
                ((0, 3), (1, 2)),
                ((0, 4), (1, 3)),
                ((0, 5), (1, 4)),
            ),
            4,
        ),
        _make_class(
            ("GradU", "GradU"),
            "divU2",
            "Surviving candidate",
            "survives current rules",
            REP_GRADU_GRADU_DIV,
            4,
        ),
        _make_class(
            ("GradU", "GradU"),
            "gradU2",
            "Surviving candidate",
            "normal form",
            REP_GRADU_GRADU_GRAD,
            4,
        ),
        _make_class(
            ("GradU", "GradU"),
            "mixedGradU2",
            "Surviving candidate",
            "survives current rules",
            REP_GRADU_GRADU_MIXED,
            4,
        ),
        _make_class(
            ("U", "U", "U", "U"),
            "U2^2",
            "Surviving candidate",
            "normal form",
            REP_UUUU_PAIR,
            4,
        ),
        _make_class(
            ("U", "U", "U", "U"),
            "U4_chain",
            "Surviving candidate",
            "survives current rules",
            REP_UUUU_CHAIN,
            4,
        ),
        _make_class(
            ("U", "U", "U", "U"),
            "U4_bridge",
            "Surviving candidate",
            "survives current rules",
            REP_UUUU_BRIDGE,
            4,
        ),
        _make_class(
            ("U", "U", "U", "U"),
            "U4_tetra",
            "Surviving candidate",
            "survives current rules",
            REP_UUUU_TETRA,
            4,
        ),
        _make_class(
            ("U", "U", "U", "U"),
            "U4_balanced",
            "Surviving candidate",
            "survives current rules",
            REP_UUUU_BALANCED,
            4,
        ),
    )
    return tuple(sorted(classes, key=lambda item: (item.weight, item.signature, item.label)))


def enumerate_r5_family_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    if max_weight != MAX_WEIGHT:
        raise ValueError("The current rank-5 family audit is only fixed at Delta <= 4.")
    return _r5_family_classes()


def enumerate_r5_sector_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    baseline = enumerate_contraction_classes(max_weight=max_weight)
    family = enumerate_r5_family_classes(max_weight=max_weight)
    return tuple(sorted(baseline + family, key=lambda item: (item.weight, item.signature, item.label)))


def r5_surviving_classes(max_weight: int = MAX_WEIGHT) -> tuple[ContractionClass, ...]:
    return tuple(
        item
        for item in enumerate_r5_sector_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    )


def _is_self_witness(item: ContractionClass) -> bool:
    return all(name in R5_FAMILY_NAMES for name in item.signature)


def r5_summary(max_weight: int = MAX_WEIGHT) -> R5Summary:
    baseline_labels = {
        item.label
        for item in enumerate_contraction_classes(max_weight=max_weight)
        if item.classification == "Surviving candidate"
    }
    survivors = r5_surviving_classes(max_weight=max_weight)
    new_survivors = tuple(item for item in survivors if item.label not in baseline_labels)
    self_candidates = tuple(item for item in new_survivors if _is_self_witness(item))
    mixed_candidates = tuple(item for item in new_survivors if not _is_self_witness(item))
    first_self = min(self_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_mixed = min(mixed_candidates, key=lambda item: (item.weight, item.label), default=None)
    first_new = min(new_survivors, key=lambda item: (item.weight, item.label), default=None)
    return R5Summary(
        total_classes=len(enumerate_r5_sector_classes(max_weight=max_weight)),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=tuple(item.label for item in new_survivors),
        first_self_witness=None if first_self is None else first_self.label,
        first_mixed_witness=None if first_mixed is None else first_mixed.label,
        smallest_new_witness=None if first_new is None else first_new.label,
    )


def r5_sector_report(max_weight: int = MAX_WEIGHT) -> str:
    classes = enumerate_r5_sector_classes(max_weight=max_weight)
    survivors = r5_surviving_classes(max_weight=max_weight)
    summary = r5_summary(max_weight=max_weight)

    lines = [
        "Delta<=4 rank-5 family audit",
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
            "- The primitive family audited here is a genuine parity-even fully symmetric trace-free rank-5 tensor family U_ijklm.",
            "- Trace descendants reducible to lower-rank audited families, and derivative-generated rank-5 blocks such as GradQ, GradT, higher-gradient vector descendants, or scalar descendants dressed only by derivatives, are excluded from the primitive-family definition and are not double-counted here.",
            "- The current audited-set thresholds for R2, R0a, R0b, R1, Rodd+, and Reven4+ already reduce the Delta<=4 baseline back to the electric sector, so the first mixed rank-5 witness is read against that electric baseline.",
            f"- First self witness: {summary.first_self_witness}.",
            f"- First mixed witness: {summary.first_mixed_witness}.",
            f"- Smallest new witness beyond the enlarged audited-set baseline: {summary.smallest_new_witness}.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(r5_sector_report())


if __name__ == "__main__":
    main()
