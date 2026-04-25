from __future__ import annotations

from dataclasses import dataclass

from eb_sector_delta4 import eb_summary
from es_sector_delta4 import es_summary
from r1_sector_delta4 import r1_summary
from r3_sector_delta4 import r3_summary
from r4_sector_delta4 import r4_summary
from shift_scalar_sector_delta4 import shift_scalar_summary


@dataclass(frozen=True)
class FamilyWitness:
    family_class: str
    audited_profile: str
    smallest_surviving_operator: str
    audited_instance: str
    weight: int
    obstructed_theorem_layer: str
    finite_family_collapse_obstructed: bool
    harmless_without_extra_assumptions: bool


@dataclass(frozen=True)
class FamilyWitnessSummary:
    entries: tuple[FamilyWitness, ...]
    any_harmless: bool


def family_witness_entries() -> tuple[FamilyWitness, ...]:
    eb = eb_summary()
    es = es_summary()
    shift = shift_scalar_summary()
    r1 = r1_summary()
    r3 = r3_summary()
    r4 = r4_summary()

    entries = (
        FamilyWitness(
            family_class="rank2_stf",
            audited_profile="unsuppressed weight-1 STF rank-2 family admitted to parity-even scalar sector",
            smallest_surviving_operator="X2",
            audited_instance=eb.smallest_new_survivor or "B2",
            weight=2,
            obstructed_theorem_layer="physically justified minimal-sector theorem",
            finite_family_collapse_obstructed=False,
            harmless_without_extra_assumptions=False,
        ),
        FamilyWitness(
            family_class="rank0_scalar_unsuppressed",
            audited_profile="unsuppressed bare scalar family",
            smallest_surviving_operator=es.smallest_new_survivor or "S",
            audited_instance=es.smallest_new_survivor or "S",
            weight=1,
            obstructed_theorem_layer="promotion of corrected E/B exact-current-set theorem to a physically justified minimal-sector theorem",
            finite_family_collapse_obstructed=False,
            harmless_without_extra_assumptions=False,
        ),
        FamilyWitness(
            family_class="rank0_scalar_derivative_only",
            audited_profile="shift-symmetric or derivative-only scalar family",
            smallest_surviving_operator=shift.canonical_new_survivor or "dotS2",
            audited_instance=", ".join(shift.first_new_labels),
            weight=shift.first_new_weight or 4,
            obstructed_theorem_layer="rescue of minimal-sector uniqueness after removing bare S",
            finite_family_collapse_obstructed=False,
            harmless_without_extra_assumptions=False,
        ),
        FamilyWitness(
            family_class="rank1_vector",
            audited_profile="genuine local parity-even vector family excluding derivative-generated vector blocks",
            smallest_surviving_operator=r1.first_self_witness or "V2",
            audited_instance=r1.smallest_new_witness or "V2",
            weight=2,
            obstructed_theorem_layer="promotion of the current audited-set result to MVP-envelope sufficiency",
            finite_family_collapse_obstructed=False,
            harmless_without_extra_assumptions=False,
        ),
        FamilyWitness(
            family_class="rank3_tensor_stf",
            audited_profile="genuine local parity-even fully symmetric trace-free rank-3 family excluding derivative-generated rank-3 blocks",
            smallest_surviving_operator=r3.first_self_witness or "T2",
            audited_instance=r3.smallest_new_witness or "T2",
            weight=2,
            obstructed_theorem_layer="promotion of the enlarged audited-set result to MVP-envelope sufficiency",
            finite_family_collapse_obstructed=False,
            harmless_without_extra_assumptions=False,
        ),
        FamilyWitness(
            family_class="rank4_tensor_stf",
            audited_profile="genuine local parity-even fully symmetric trace-free rank-4 family excluding trace descendants and derivative-generated rank-4 blocks",
            smallest_surviving_operator=r4.first_self_witness or "Q2",
            audited_instance=r4.smallest_new_witness or "Q2",
            weight=2,
            obstructed_theorem_layer="promotion of the enlarged audited-set result to MVP-envelope sufficiency; enlarged audited-set composition must be re-closed next",
            finite_family_collapse_obstructed=False,
            harmless_without_extra_assumptions=False,
        ),
    )
    return entries


def family_witness_summary() -> FamilyWitnessSummary:
    entries = family_witness_entries()
    return FamilyWitnessSummary(
        entries=entries,
        any_harmless=any(entry.harmless_without_extra_assumptions for entry in entries),
    )


def family_witness_report() -> str:
    summary = family_witness_summary()
    lines = [
        "Family-admission witness map",
        "",
        "Audited family classes:",
    ]
    for entry in summary.entries:
        lines.extend(
            [
                f"- family_class: {entry.family_class}",
                f"  profile: {entry.audited_profile}",
                f"  smallest_witness: {entry.smallest_surviving_operator} (audited as {entry.audited_instance})",
                f"  weight: {entry.weight}",
                f"  obstructs: {entry.obstructed_theorem_layer}",
                f"  obstructs_finite_family_collapse: {entry.finite_family_collapse_obstructed}",
                f"  harmless_without_extra_assumptions: {entry.harmless_without_extra_assumptions}",
            ]
        )
    lines.extend(
        [
            "",
            f"Any audited family class harmless without extra assumptions: {summary.any_harmless}",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(family_witness_report())


if __name__ == "__main__":
    main()
