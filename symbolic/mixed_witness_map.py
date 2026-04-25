from __future__ import annotations

from dataclasses import dataclass


DELTA_MAX = 4


@dataclass(frozen=True)
class MixedWitnessEntry:
    family_class: str
    first_self_witness: str
    self_weight: int
    first_mixed_witness: str
    mixed_weight: int
    w_min: int
    current_threshold_statement: str
    theorem_layer_obstructed: str
    sharpness_status: str
    confidence_level: str


@dataclass(frozen=True)
class MixedWitnessSummary:
    delta_max: int
    entries: tuple[MixedWitnessEntry, ...]
    mixed_dominant_family_class: str | None


def mixed_witness_entries(delta_max: int = DELTA_MAX) -> tuple[MixedWitnessEntry, ...]:
    if delta_max != 4:
        raise ValueError("The current mixed-witness audit is only fixed at Delta_max = 4.")

    return (
        MixedWitnessEntry(
            family_class="rank2_stf",
            first_self_witness="X2",
            self_weight=2,
            first_mixed_witness="EX",
            mixed_weight=2,
            w_min=2,
            current_threshold_statement="w_X >= 4, unless an explicit EX = 0-type rule removes the mixed quadratic witness",
            theorem_layer_obstructed=(
                "promotion of electric-only exact-current-set theorem to a physically justified minimal-sector theorem"
            ),
            sharpness_status="mixed-aware; current unsuppressed case is tied",
            confidence_level="Proven",
        ),
        MixedWitnessEntry(
            family_class="rank0_scalar_unsuppressed",
            first_self_witness="S",
            self_weight=1,
            first_mixed_witness="SE2",
            mixed_weight=3,
            w_min=1,
            current_threshold_statement="w_S >= 5, or explicit rule forbidding bare S",
            theorem_layer_obstructed=(
                "promotion of corrected E/B exact-current-set theorem to a physically justified minimal-sector theorem"
            ),
            sharpness_status="self-only",
            confidence_level="Proven",
        ),
        MixedWitnessEntry(
            family_class="rank0_scalar_derivative_only",
            first_self_witness="dotS2",
            self_weight=4,
            first_mixed_witness="DtS_E2",
            mixed_weight=4,
            w_min=4,
            current_threshold_statement="w_D >= 3, or explicit rule removing the mixed derivative witnesses",
            theorem_layer_obstructed="rescue of minimal-sector uniqueness after removing bare S",
            sharpness_status="tied-sharp",
            confidence_level="Proven",
        ),
        MixedWitnessEntry(
            family_class="rank1_vector",
            first_self_witness="V2",
            self_weight=2,
            first_mixed_witness="EVV",
            mixed_weight=3,
            w_min=2,
            current_threshold_statement="w_V >= 3, or explicit rule excluding or absorbing the primitive vector family",
            theorem_layer_obstructed="promotion of the current audited-set result to MVP-envelope sufficiency",
            sharpness_status="self-only",
            confidence_level="Proven",
        ),
        MixedWitnessEntry(
            family_class="rank3_tensor_stf",
            first_self_witness="T2",
            self_weight=2,
            first_mixed_witness="ETT",
            mixed_weight=3,
            w_min=2,
            current_threshold_statement="w_T >= 3, or explicit rule excluding or absorbing the primitive rank-3 family",
            theorem_layer_obstructed="promotion of the enlarged audited-set result to MVP-envelope sufficiency",
            sharpness_status="self-only",
            confidence_level="Proven",
        ),
        MixedWitnessEntry(
            family_class="rank4_tensor_stf",
            first_self_witness="Q2",
            self_weight=2,
            first_mixed_witness="EQQ",
            mixed_weight=3,
            w_min=2,
            current_threshold_statement="w_Q >= 3, or explicit rule excluding or absorbing the primitive rank-4 family",
            theorem_layer_obstructed="promotion of the enlarged audited-set result to MVP-envelope sufficiency",
            sharpness_status="self-only",
            confidence_level="Proven",
        ),
    )


def mixed_witness_summary(delta_max: int = DELTA_MAX) -> MixedWitnessSummary:
    entries = mixed_witness_entries(delta_max)
    mixed_dominant = next(
        (entry.family_class for entry in entries if entry.mixed_weight < entry.self_weight),
        None,
    )
    return MixedWitnessSummary(
        delta_max=delta_max,
        entries=entries,
        mixed_dominant_family_class=mixed_dominant,
    )


def mixed_witness_report(delta_max: int = DELTA_MAX) -> str:
    summary = mixed_witness_summary(delta_max)
    lines = [
        "family_class\tfirst_self_witness\tself_weight\tfirst_mixed_witness\tmixed_weight\tW_min\tcurrent_threshold_statement\ttheorem_layer_obstructed",
    ]
    for entry in summary.entries:
        lines.append(
            "\t".join(
                (
                    entry.family_class,
                    entry.first_self_witness,
                    str(entry.self_weight),
                    entry.first_mixed_witness,
                    str(entry.mixed_weight),
                    str(entry.w_min),
                    entry.current_threshold_statement,
                    entry.theorem_layer_obstructed,
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(mixed_witness_report())


if __name__ == "__main__":
    main()
