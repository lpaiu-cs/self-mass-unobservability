from __future__ import annotations

from dataclasses import dataclass

from eb_sector_delta4 import eb_summary
from eb_survivor_rank_check import eb_rank_summary
from enumerate_contractions_delta4 import enumerate_contraction_classes
from es_sector_delta4 import es_summary
from es_survivor_rank_check import es_rank_summary
from family_envelope_census import family_envelope_summary
from irrep_family_census import irrep_family_summary
from r1_sector_delta4 import r1_summary
from r1_survivor_rank_check import r1_rank_summary
from r3_sector_delta4 import r3_summary
from r3_survivor_rank_check import r3_rank_summary
from r4_sector_delta4 import r4_summary
from r4_survivor_rank_check import r4_rank_summary
from r5_sector_delta4 import r5_summary
from r5_survivor_rank_check import r5_rank_summary
from r6_sector_delta4 import r6_summary
from r6_survivor_rank_check import r6_rank_summary
from shift_scalar_sector_delta4 import shift_scalar_summary
from survivor_rank_check import rank_summary


DELTA_MAX = 4
REDUCTION_LAYERS = (
    "irreducible family-envelope closure",
    "trace-descendant absorption",
    "total derivatives",
    "lower-order EOM",
    "explicit algebraic identities",
    "sector-specific linear-dependence quotient when present",
)


@dataclass(frozen=True)
class OperatorCountEntry:
    catalog_id: str
    primitive_family_classes: tuple[str, ...]
    candidate_operator_count: int
    reduced_operator_count: int
    reduced_count_kind: str
    finite_before_reduction: bool
    finite_after_reduction: bool
    reduction_layers_applied: tuple[str, ...]
    note: str


@dataclass(frozen=True)
class FixedFamilyOperatorCountSummary:
    delta_max: int
    theorem_domain_classes: tuple[str, ...]
    baseline_sector: str
    envelope_closed: bool
    candidate_operator_space_finite: bool
    reduced_operator_space_finite: bool
    reduction_layers_applied: tuple[str, ...]
    entries: tuple[OperatorCountEntry, ...]


def fixed_family_operator_entries() -> tuple[OperatorCountEntry, ...]:
    electric = rank_summary()
    eb = eb_summary()
    eb_rank = eb_rank_summary()
    es = es_summary()
    es_rank = es_rank_summary()
    shift = shift_scalar_summary()
    r1 = r1_summary()
    r1_rank = r1_rank_summary()
    r3 = r3_summary()
    r3_rank = r3_rank_summary()
    r4 = r4_summary()
    r4_rank = r4_rank_summary()
    r5 = r5_summary()
    r5_rank = r5_rank_summary()
    r6 = r6_summary()
    r6_rank = r6_rank_summary()

    return (
        OperatorCountEntry(
            catalog_id="electric_exact_current_set",
            primitive_family_classes=("baseline electric sector",),
            candidate_operator_count=len(enumerate_contraction_classes()),
            reduced_operator_count=electric.total_rank,
            reduced_count_kind="exact polynomial rank of the corrected seven-element basis",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Current exact-current-set electric witness of finite scalar normal-form closure.",
        ),
        OperatorCountEntry(
            catalog_id="rank2_stf_special_case",
            primitive_family_classes=("STF2 special case",),
            candidate_operator_count=eb.total_classes,
            reduced_operator_count=eb_rank.corrected_rank,
            reduced_count_kind="corrected exact rank after removing the explicit mixed quartic dependence",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Representative corrected E/B sector for the mixed-aware rank-2 STF special case.",
        ),
        OperatorCountEntry(
            catalog_id="scalar_bare_source_extension",
            primitive_family_classes=("Scalar", "STF2 special case"),
            candidate_operator_count=es.total_classes,
            reduced_operator_count=es_rank.rank,
            reduced_count_kind="exact rank of the corrected E/B+scalar survivor set",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Representative scalar-family sector with bare source admitted.",
        ),
        OperatorCountEntry(
            catalog_id="scalar_derivative_only_extension",
            primitive_family_classes=("Scalar derivative-only", "STF2 special case"),
            candidate_operator_count=shift.total_classes,
            reduced_operator_count=len(shift.surviving_labels),
            reduced_count_kind="explicit surviving-label count under the derivative-only reduction rules",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Derivative-only scalar representative sector; finiteness does not rely on a further independence theorem.",
        ),
        OperatorCountEntry(
            catalog_id="vector_representative_sector",
            primitive_family_classes=("Vector",),
            candidate_operator_count=r1.total_classes,
            reduced_operator_count=r1_rank.rank,
            reduced_count_kind="exact rank of the corrected vector-extended survivor set",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Representative genuine vector-family sector.",
        ),
        OperatorCountEntry(
            catalog_id="stf_rank3_representative_sector",
            primitive_family_classes=("STFge3",),
            candidate_operator_count=r3.total_classes,
            reduced_operator_count=r3_rank.rank,
            reduced_count_kind="exact rank after verified null-relation quotient",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Representative odd-rank STF sector at rank 3.",
        ),
        OperatorCountEntry(
            catalog_id="stf_rank4_representative_sector",
            primitive_family_classes=("STFge3",),
            candidate_operator_count=r4.total_classes,
            reduced_operator_count=r4_rank.rank,
            reduced_count_kind="sample-verified corrected quotient rank after the exhaustive rank-4 patch",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Representative even-rank STF exception sector at rank 4.",
        ),
        OperatorCountEntry(
            catalog_id="stf_rank5_representative_sector",
            primitive_family_classes=("STFge3",),
            candidate_operator_count=r5.total_classes,
            reduced_operator_count=r5_rank.rank,
            reduced_count_kind="sample-verified corrected quotient rank after null-relation extraction",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Representative odd-rank STF sector at rank 5.",
        ),
        OperatorCountEntry(
            catalog_id="stf_rank6_representative_sector",
            primitive_family_classes=("STFge3",),
            candidate_operator_count=r6.total_classes,
            reduced_operator_count=r6_rank.rank,
            reduced_count_kind="sample-stable corrected quotient rank after signature-wise null-relation extraction",
            finite_before_reduction=True,
            finite_after_reduction=True,
            reduction_layers_applied=REDUCTION_LAYERS,
            note="Representative even-rank STF sector at rank 6.",
        ),
    )


def fixed_family_operator_count_summary(
    delta_max: int = DELTA_MAX,
) -> FixedFamilyOperatorCountSummary:
    if delta_max != 4:
        raise ValueError("The current fixed-family operator count is only fixed at Delta <= 4.")
    irrep = irrep_family_summary(delta_max=delta_max)
    envelope = family_envelope_summary(delta_max=delta_max)
    entries = fixed_family_operator_entries()
    return FixedFamilyOperatorCountSummary(
        delta_max=delta_max,
        theorem_domain_classes=("Scalar", "Vector", "STF2", "STFge3"),
        baseline_sector="electric parity-even free-fall scalar sector",
        envelope_closed=irrep.envelope_closes_on_audited_classes and envelope.envelope_closed,
        candidate_operator_space_finite=all(entry.finite_before_reduction for entry in entries),
        reduced_operator_space_finite=all(entry.finite_after_reduction for entry in entries),
        reduction_layers_applied=REDUCTION_LAYERS,
        entries=entries,
    )


def fixed_family_operator_count_report(delta_max: int = DELTA_MAX) -> str:
    summary = fixed_family_operator_count_summary(delta_max=delta_max)
    lines = [
        "key\tvalue",
        f"delta_max\t{summary.delta_max}",
        f"theorem_domain_classes\t{','.join(summary.theorem_domain_classes)}",
        f"baseline_sector\t{summary.baseline_sector}",
        f"envelope_closed\t{str(summary.envelope_closed).lower()}",
        f"candidate_operator_space_finite\t{str(summary.candidate_operator_space_finite).lower()}",
        f"reduced_operator_space_finite\t{str(summary.reduced_operator_space_finite).lower()}",
        f"reduction_layers_applied\t{' | '.join(summary.reduction_layers_applied)}",
        "",
        (
            "catalog_id\tprimitive_family_classes\tcandidate_operator_count\t"
            "reduced_operator_count\treduced_count_kind\tfinite_before_reduction\t"
            "finite_after_reduction\tnote"
        ),
    ]
    for entry in summary.entries:
        lines.append(
            "\t".join(
                (
                    entry.catalog_id,
                    ",".join(entry.primitive_family_classes),
                    str(entry.candidate_operator_count),
                    str(entry.reduced_operator_count),
                    entry.reduced_count_kind,
                    str(entry.finite_before_reduction).lower(),
                    str(entry.finite_after_reduction).lower(),
                    entry.note,
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(fixed_family_operator_count_report())


if __name__ == "__main__":
    main()
