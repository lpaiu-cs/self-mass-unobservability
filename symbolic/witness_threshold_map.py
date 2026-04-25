from __future__ import annotations

from dataclasses import dataclass


DELTA_MAX = 4


@dataclass(frozen=True)
class WitnessThresholdEntry:
    family_class: str
    witness: str
    witness_weight: int
    self_only_lower_bound_formula: str
    self_only_lower_bound_at_delta4: str
    necessary_threshold_formula: str
    necessary_threshold_at_delta4: str
    threshold_type: str
    theorem_layer_at_risk: str
    sufficient_for_uniqueness: bool


@dataclass(frozen=True)
class WitnessThresholdSummary:
    delta_max: int
    entries: tuple[WitnessThresholdEntry, ...]


def witness_threshold_entries(delta_max: int = DELTA_MAX) -> tuple[WitnessThresholdEntry, ...]:
    if delta_max != 4:
        raise ValueError("The current audited threshold map is only fixed at Delta_max = 4.")

    return (
        WitnessThresholdEntry(
            family_class="rank2_stf",
            witness="X2",
            witness_weight=2,
            self_only_lower_bound_formula="2 w_X > Delta_max",
            self_only_lower_bound_at_delta4="w_X >= 3",
            necessary_threshold_formula="min(2 w_X, w_X + 1) > Delta_max",
            necessary_threshold_at_delta4="w_X >= 4, unless an explicit EX = 0-type rule removes the mixed quadratic witness",
            threshold_type="mixed-aware",
            theorem_layer_at_risk="minimal-sector uniqueness",
            sufficient_for_uniqueness=False,
        ),
        WitnessThresholdEntry(
            family_class="rank0_scalar_unsuppressed",
            witness="S",
            witness_weight=1,
            self_only_lower_bound_formula="w_S > Delta_max",
            self_only_lower_bound_at_delta4="w_S >= 5",
            necessary_threshold_formula="w_S > Delta_max",
            necessary_threshold_at_delta4="w_S >= 5, or explicit rule forbidding bare S",
            threshold_type="self-only",
            theorem_layer_at_risk="minimal-sector uniqueness",
            sufficient_for_uniqueness=False,
        ),
        WitnessThresholdEntry(
            family_class="rank0_scalar_derivative_only",
            witness="dotS2",
            witness_weight=4,
            self_only_lower_bound_formula="2 w_D > Delta_max",
            self_only_lower_bound_at_delta4="w_D >= 3",
            necessary_threshold_formula="min(2 w_D, w_D + 2) > Delta_max",
            necessary_threshold_at_delta4="w_D >= 3, or explicit rule removing the mixed derivative witnesses",
            threshold_type="tied-sharp",
            theorem_layer_at_risk="minimal-sector uniqueness after removing bare S",
            sufficient_for_uniqueness=False,
        ),
        WitnessThresholdEntry(
            family_class="rank1_vector",
            witness="V2",
            witness_weight=2,
            self_only_lower_bound_formula="2 w_V > Delta_max",
            self_only_lower_bound_at_delta4="w_V >= 3",
            necessary_threshold_formula="2 w_V > Delta_max",
            necessary_threshold_at_delta4="w_V >= 3, or explicit rule excluding or absorbing the primitive vector family",
            threshold_type="self-only",
            theorem_layer_at_risk="MVP-envelope sufficiency",
            sufficient_for_uniqueness=False,
        ),
        WitnessThresholdEntry(
            family_class="rank3_tensor_stf",
            witness="T2",
            witness_weight=2,
            self_only_lower_bound_formula="2 w_T > Delta_max",
            self_only_lower_bound_at_delta4="w_T >= 3",
            necessary_threshold_formula="2 w_T > Delta_max",
            necessary_threshold_at_delta4="w_T >= 3, or explicit rule excluding or absorbing the primitive rank-3 family",
            threshold_type="self-only",
            theorem_layer_at_risk="MVP-envelope sufficiency after the R1-closed audited set",
            sufficient_for_uniqueness=False,
        ),
    )


def witness_threshold_summary(delta_max: int = DELTA_MAX) -> WitnessThresholdSummary:
    return WitnessThresholdSummary(delta_max=delta_max, entries=witness_threshold_entries(delta_max))


def witness_threshold_report(delta_max: int = DELTA_MAX) -> str:
    summary = witness_threshold_summary(delta_max)
    lines = [
        "Witness-threshold map",
        "",
        f"Delta_max: {summary.delta_max}",
        "",
        "Audited class thresholds:",
    ]
    for entry in summary.entries:
        lines.extend(
            [
                f"- family_class: {entry.family_class}",
                f"  witness: {entry.witness}",
                f"  witness_weight: {entry.witness_weight}",
                f"  self_only_lower_bound_formula: {entry.self_only_lower_bound_formula}",
                f"  self_only_lower_bound_at_delta4: {entry.self_only_lower_bound_at_delta4}",
                f"  necessary_threshold_formula: {entry.necessary_threshold_formula}",
                f"  necessary_threshold_at_delta4: {entry.necessary_threshold_at_delta4}",
                f"  threshold_type: {entry.threshold_type}",
                f"  theorem_layer_at_risk: {entry.theorem_layer_at_risk}",
                f"  sufficient_for_uniqueness: {entry.sufficient_for_uniqueness}",
            ]
        )
    return "\n".join(lines)


def main() -> None:
    print(witness_threshold_report())


if __name__ == "__main__":
    main()
