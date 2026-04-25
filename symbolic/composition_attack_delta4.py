from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

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
R1_BLOCK_TYPES = (
    BlockType("V", 3, 1, parity=0),
    BlockType("DtV", 4, 1, parity=0),
    BlockType("GradV", 4, 2, parity=1),
)
R3_BLOCK_TYPES = (
    BlockType(
        "T",
        3,
        3,
        parity=0,
        sym_groups=((0, 1), (1, 2)),
        tracefree_pairs=((0, 1), (0, 2), (1, 2)),
    ),
    BlockType(
        "DtT",
        4,
        3,
        parity=0,
        sym_groups=((0, 1), (1, 2)),
        tracefree_pairs=((0, 1), (0, 2), (1, 2)),
    ),
    BlockType(
        "GradT",
        4,
        4,
        parity=1,
        sym_groups=((1, 2), (2, 3)),
        tracefree_pairs=((1, 2), (1, 3), (2, 3)),
    ),
)
R4_BLOCK_TYPES = (
    BlockType(
        "Q",
        3,
        4,
        parity=0,
        sym_groups=((0, 1), (1, 2), (2, 3)),
        tracefree_pairs=((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)),
    ),
    BlockType(
        "DtQ",
        4,
        4,
        parity=0,
        sym_groups=((0, 1), (1, 2), (2, 3)),
        tracefree_pairs=((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)),
    ),
    BlockType(
        "GradQ",
        4,
        5,
        parity=1,
        sym_groups=((1, 2), (2, 3), (3, 4)),
        tracefree_pairs=((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)),
    ),
)
R5_BLOCK_TYPES = (
    BlockType(
        "U",
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
    BlockType(
        "DtU",
        4,
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
        4,
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
)

FAMILY_BLOCKS = {
    "R2": R2_BLOCK_TYPES,
    "R0a": R0A_BLOCK_TYPES,
    "R0b": R0B_BLOCK_TYPES,
    "R1": R1_BLOCK_TYPES,
    "Rodd+": R3_BLOCK_TYPES,
    "Reven4+": R4_BLOCK_TYPES,
    "Rodd5+": R5_BLOCK_TYPES,
}

FAMILY_THRESHOLDS = {
    "R2": "w_X >= 4 unless an explicit EX = 0-type rule removes the mixed quadratic witness",
    "R0a": "w_S >= 5, or explicit exclusion of bare S",
    "R0b": "w_D >= 3, or explicit rule removing the mixed derivative witnesses",
    "R1": "w_V >= 3, or explicit rule excluding or absorbing the primitive vector family",
    "Rodd+": "w_T >= 3, or explicit rule excluding or absorbing the primitive rank-3 STF family",
    "Reven4+": "w_Q >= 3, or explicit rule excluding or absorbing the primitive rank-4 STF family",
    "Rodd5+": "w_U >= 3, or explicit rule excluding or absorbing the primitive rank-5 STF family",
}

FAMILY_NAMES = {
    "X",
    "S",
    "DtS",
    "GradS",
    "Dt2S",
    "V",
    "DtV",
    "GradV",
    "T",
    "DtT",
    "GradT",
    "Q",
    "DtQ",
    "GradQ",
    "U",
    "DtU",
    "GradU",
}

PRE_R1_AUDITED_SET_CASES = (
    ("R2", "R0a"),
    ("R2", "R0b"),
    ("R0a", "R0b"),
    ("R2", "R0a", "R0b"),
)

PRE_RODD_AUDITED_SET_CASES = (
    ("R1", "R2"),
    ("R1", "R0a"),
    ("R1", "R0b"),
    ("R1", "R2", "R0a"),
    ("R1", "R2", "R0b"),
    ("R1", "R0a", "R0b"),
    ("R2", "R0a", "R0b", "R1"),
)

PRE_REVEN_AUDITED_SET_CASES = (
    ("Rodd+", "R2"),
    ("Rodd+", "R0a"),
    ("Rodd+", "R0b"),
    ("Rodd+", "R1"),
    ("Rodd+", "R2", "R0a"),
    ("Rodd+", "R2", "R0b"),
    ("Rodd+", "R2", "R1"),
    ("Rodd+", "R0a", "R0b"),
    ("Rodd+", "R0a", "R1"),
    ("Rodd+", "R0b", "R1"),
    ("Rodd+", "R2", "R0a", "R0b"),
    ("Rodd+", "R2", "R0a", "R1"),
    ("Rodd+", "R2", "R0b", "R1"),
    ("Rodd+", "R0a", "R0b", "R1"),
    ("R2", "R0a", "R0b", "R1", "Rodd+"),
)

PRE_REVEN_BASE_FAMILIES = ("R2", "R0a", "R0b", "R1", "Rodd+")

PRE_RODD5_PAIRWISE_CASES = (
    ("Reven4+", "R2"),
    ("Reven4+", "R0a"),
    ("Reven4+", "R0b"),
    ("Reven4+", "R1"),
    ("Reven4+", "Rodd+"),
)
PRE_RODD5_TRIPLE_CASES = tuple(
    ("Reven4+",) + combo for combo in combinations(PRE_REVEN_BASE_FAMILIES, 2)
)
PRE_RODD5_QUADRUPLE_CASES = tuple(
    ("Reven4+",) + combo for combo in combinations(PRE_REVEN_BASE_FAMILIES, 3)
)
PRE_RODD5_QUINTUPLE_CASES = tuple(
    ("Reven4+",) + combo for combo in combinations(PRE_REVEN_BASE_FAMILIES, 4)
)
PRE_RODD5_FULL_SET_CASES = (
    PRE_REVEN_BASE_FAMILIES + ("Reven4+",),
)
PRE_RODD5_AUDITED_SET_CASES = (
    PRE_RODD5_PAIRWISE_CASES
    + PRE_RODD5_TRIPLE_CASES
    + PRE_RODD5_QUADRUPLE_CASES
    + PRE_RODD5_QUINTUPLE_CASES
    + PRE_RODD5_FULL_SET_CASES
)

PRE_RODD5_BASE_FAMILIES = PRE_REVEN_BASE_FAMILIES + ("Reven4+",)

POST_RODD5_PAIRWISE_CASES = (
    ("Rodd5+", "R2"),
    ("Rodd5+", "R0a"),
    ("Rodd5+", "R0b"),
    ("Rodd5+", "R1"),
    ("Rodd5+", "Rodd+"),
    ("Rodd5+", "Reven4+"),
)
POST_RODD5_TRIPLE_CASES = tuple(
    ("Rodd5+",) + combo for combo in combinations(PRE_RODD5_BASE_FAMILIES, 2)
)
POST_RODD5_QUADRUPLE_CASES = tuple(
    ("Rodd5+",) + combo for combo in combinations(PRE_RODD5_BASE_FAMILIES, 3)
)
POST_RODD5_QUINTUPLE_CASES = tuple(
    ("Rodd5+",) + combo for combo in combinations(PRE_RODD5_BASE_FAMILIES, 4)
)
POST_RODD5_SIX_FAMILY_CASES = tuple(
    ("Rodd5+",) + combo for combo in combinations(PRE_RODD5_BASE_FAMILIES, 5)
)
POST_RODD5_FULL_SET_CASES = (
    PRE_RODD5_BASE_FAMILIES + ("Rodd5+",),
)
POST_RODD5_AUDITED_SET_CASES = (
    POST_RODD5_PAIRWISE_CASES
    + POST_RODD5_TRIPLE_CASES
    + POST_RODD5_QUADRUPLE_CASES
    + POST_RODD5_QUINTUPLE_CASES
    + POST_RODD5_SIX_FAMILY_CASES
    + POST_RODD5_FULL_SET_CASES
)


@dataclass(frozen=True)
class CompositionCaseSummary:
    case_scope: str
    combination_kind: str
    combo_id: str
    families: tuple[str, ...]
    imposed_thresholds: tuple[str, ...]
    total_classes: int
    surviving_labels: tuple[str, ...]
    new_surviving_labels: tuple[str, ...]
    sufficient: bool
    smallest_surviving_cross_family: str | None
    smallest_surviving_cross_family_weight: int | None


@dataclass(frozen=True)
class CompositionSummary:
    baseline_survivors: tuple[str, ...]
    pre_r1_audited_set_joint_sufficient: bool
    pre_r1_smallest_surviving_cross_family: str | None
    pre_rodd_audited_set_joint_sufficient: bool
    pre_rodd_smallest_surviving_cross_family: str | None
    pre_reven_audited_set_joint_sufficient: bool
    pre_reven_smallest_surviving_cross_family: str | None
    pre_rodd5_audited_set_joint_sufficient: bool
    pre_rodd5_smallest_surviving_cross_family: str | None
    post_rodd5_audited_set_joint_sufficient: bool
    post_rodd5_smallest_surviving_cross_family: str | None
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


def case_scope(families: tuple[str, ...]) -> str:
    if "Rodd5+" in families:
        return "post_rodd5_audited_set"
    if "Reven4+" in families:
        return "pre_rodd5_audited_set"
    if "Rodd+" in families:
        return "pre_reven_audited_set"
    if "R1" in families:
        return "pre_rodd_audited_set"
    return "pre_r1_audited_set"


def combination_kind(families: tuple[str, ...]) -> str:
    scope = case_scope(families)
    size = len(families)
    if scope == "post_rodd5_audited_set":
        if size == 2:
            return "pairwise"
        if size == 3:
            return "triple"
        if size == 4:
            return "quadruple"
        if size == 5:
            return "quintuple"
        if size == 6:
            return "six-family"
        if size == 7:
            return "full-set"
    if scope == "pre_rodd5_audited_set":
        if size == 2:
            return "pairwise"
        if size == 3:
            return "triple"
        if size == 4:
            return "quadruple"
        if size == 5:
            return "quintuple"
        if size == 6:
            return "full-set"
    if size == 2:
        return "pairwise"
    if size == 3:
        return "triple"
    if size == 4:
        return "quadruple"
    if size == 5:
        return "full-set"
    raise ValueError(f"Unsupported composition arity {size} for families {families}")


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
        case_scope=case_scope(families),
        combination_kind=combination_kind(families),
        combo_id="+".join(families),
        families=families,
        imposed_thresholds=tuple(FAMILY_THRESHOLDS[family] for family in families),
        total_classes=len(classes),
        surviving_labels=tuple(item.label for item in survivors),
        new_surviving_labels=tuple(item.label for item in new_survivors),
        sufficient=not new_survivors,
        smallest_surviving_cross_family=None if smallest is None else smallest.label,
        smallest_surviving_cross_family_weight=None if smallest is None else smallest.weight,
    )


def _smallest_surviving_case(
    cases: tuple[CompositionCaseSummary, ...],
) -> CompositionCaseSummary | None:
    return min(
        (
            case
            for case in cases
            if case.smallest_surviving_cross_family is not None
            and case.smallest_surviving_cross_family_weight is not None
        ),
        key=lambda case: (
            case.smallest_surviving_cross_family_weight,
            case.smallest_surviving_cross_family or "",
        ),
        default=None,
    )


def composition_summary(max_weight: int = DELTA_MAX) -> CompositionSummary:
    pre_r1_cases = tuple(
        composition_case_summary(families, max_weight=max_weight)
        for families in PRE_R1_AUDITED_SET_CASES
    )
    pre_rodd_cases = tuple(
        composition_case_summary(families, max_weight=max_weight)
        for families in PRE_RODD_AUDITED_SET_CASES
    )
    pre_reven_cases = tuple(
        composition_case_summary(families, max_weight=max_weight)
        for families in PRE_REVEN_AUDITED_SET_CASES
    )
    pre_rodd5_cases = tuple(
        composition_case_summary(families, max_weight=max_weight)
        for families in PRE_RODD5_AUDITED_SET_CASES
    )
    post_rodd5_cases = tuple(
        composition_case_summary(families, max_weight=max_weight)
        for families in POST_RODD5_AUDITED_SET_CASES
    )

    pre_r1_smallest_case = _smallest_surviving_case(pre_r1_cases)
    pre_rodd_smallest_case = _smallest_surviving_case(pre_rodd_cases)
    pre_reven_smallest_case = _smallest_surviving_case(pre_reven_cases)
    pre_rodd5_smallest_case = _smallest_surviving_case(pre_rodd5_cases)
    post_rodd5_smallest_case = _smallest_surviving_case(post_rodd5_cases)

    return CompositionSummary(
        baseline_survivors=baseline_survivor_labels(max_weight=max_weight),
        pre_r1_audited_set_joint_sufficient=all(case.sufficient for case in pre_r1_cases),
        pre_r1_smallest_surviving_cross_family=(
            None
            if pre_r1_smallest_case is None
            else pre_r1_smallest_case.smallest_surviving_cross_family
        ),
        pre_rodd_audited_set_joint_sufficient=all(case.sufficient for case in pre_rodd_cases),
        pre_rodd_smallest_surviving_cross_family=(
            None
            if pre_rodd_smallest_case is None
            else pre_rodd_smallest_case.smallest_surviving_cross_family
        ),
        pre_reven_audited_set_joint_sufficient=all(case.sufficient for case in pre_reven_cases),
        pre_reven_smallest_surviving_cross_family=(
            None
            if pre_reven_smallest_case is None
            else pre_reven_smallest_case.smallest_surviving_cross_family
        ),
        pre_rodd5_audited_set_joint_sufficient=all(case.sufficient for case in pre_rodd5_cases),
        pre_rodd5_smallest_surviving_cross_family=(
            None
            if pre_rodd5_smallest_case is None
            else pre_rodd5_smallest_case.smallest_surviving_cross_family
        ),
        post_rodd5_audited_set_joint_sufficient=all(case.sufficient for case in post_rodd5_cases),
        post_rodd5_smallest_surviving_cross_family=(
            None
            if post_rodd5_smallest_case is None
            else post_rodd5_smallest_case.smallest_surviving_cross_family
        ),
        cases=pre_r1_cases + pre_rodd_cases + pre_reven_cases + pre_rodd5_cases + post_rodd5_cases,
    )


def composition_report(max_weight: int = DELTA_MAX) -> str:
    summary = composition_summary(max_weight=max_weight)
    lines = [
        "key\tvalue",
        (
            "pre_r1_audited_set_joint_sufficient\t"
            f"{str(summary.pre_r1_audited_set_joint_sufficient).lower()}"
        ),
        (
            "pre_r1_smallest_surviving_cross_family\t"
            f"{summary.pre_r1_smallest_surviving_cross_family or 'none'}"
        ),
        (
            "pre_rodd_audited_set_joint_sufficient\t"
            f"{str(summary.pre_rodd_audited_set_joint_sufficient).lower()}"
        ),
        (
            "pre_rodd_smallest_surviving_cross_family\t"
            f"{summary.pre_rodd_smallest_surviving_cross_family or 'none'}"
        ),
        (
            "pre_reven_audited_set_joint_sufficient\t"
            f"{str(summary.pre_reven_audited_set_joint_sufficient).lower()}"
        ),
        (
            "pre_reven_smallest_surviving_cross_family\t"
            f"{summary.pre_reven_smallest_surviving_cross_family or 'none'}"
        ),
        (
            "pre_rodd5_audited_set_joint_sufficient\t"
            f"{str(summary.pre_rodd5_audited_set_joint_sufficient).lower()}"
        ),
        (
            "pre_rodd5_smallest_surviving_cross_family\t"
            f"{summary.pre_rodd5_smallest_surviving_cross_family or 'none'}"
        ),
        (
            "post_rodd5_audited_set_joint_sufficient\t"
            f"{str(summary.post_rodd5_audited_set_joint_sufficient).lower()}"
        ),
        (
            "post_rodd5_smallest_surviving_cross_family\t"
            f"{summary.post_rodd5_smallest_surviving_cross_family or 'none'}"
        ),
        "",
        (
            "case_scope\tcombination_kind\tcombo_id\tfamilies\tresult\t"
            "smallest_surviving_cross_family\timposed_thresholds"
        ),
    ]
    for case in summary.cases:
        lines.append(
            "\t".join(
                (
                    case.case_scope,
                    case.combination_kind,
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
