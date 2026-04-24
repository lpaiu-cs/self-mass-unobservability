from __future__ import annotations

from dataclasses import dataclass


DELTA_MAX = 4


@dataclass(frozen=True)
class WitnessThresholdEntry:
    family_class: str
    witness: str
    witness_weight: int
    necessary_threshold_formula: str
    necessary_threshold_at_delta4: str
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
            necessary_threshold_formula="2 w_X > Delta_max",
            necessary_threshold_at_delta4="w_X >= 3, or explicit family exclusion/background restriction",
            theorem_layer_at_risk="minimal-sector uniqueness",
            sufficient_for_uniqueness=False,
        ),
        WitnessThresholdEntry(
            family_class="rank0_scalar_unsuppressed",
            witness="S",
            witness_weight=1,
            necessary_threshold_formula="w_S > Delta_max",
            necessary_threshold_at_delta4="w_S >= 5, or explicit rule forbidding bare S",
            theorem_layer_at_risk="minimal-sector uniqueness",
            sufficient_for_uniqueness=False,
        ),
        WitnessThresholdEntry(
            family_class="rank0_scalar_derivative_only",
            witness="dotS2",
            witness_weight=4,
            necessary_threshold_formula="raise derivative-family insertions so the current weight-4 witnesses move above Delta_max",
            necessary_threshold_at_delta4="increase derivative-family block weight above 2, for example to 3, or explicitly remove the mixed derivative witnesses",
            theorem_layer_at_risk="minimal-sector uniqueness after removing bare S",
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
                f"  necessary_threshold_formula: {entry.necessary_threshold_formula}",
                f"  necessary_threshold_at_delta4: {entry.necessary_threshold_at_delta4}",
                f"  theorem_layer_at_risk: {entry.theorem_layer_at_risk}",
                f"  sufficient_for_uniqueness: {entry.sufficient_for_uniqueness}",
            ]
        )
    return "\n".join(lines)


def main() -> None:
    print(witness_threshold_report())


if __name__ == "__main__":
    main()
