from __future__ import annotations

from dataclasses import dataclass

from mixed_witness_map import mixed_witness_summary


DELTA_MAX = 4


@dataclass(frozen=True)
class ThresholdFormulaEntry:
    family_class: str
    self_formula: str
    mixed_formula: str
    w_min_formula: str
    current_case_value: str
    sharpness_status: str
    threshold_type: str


@dataclass(frozen=True)
class ThresholdFormulaSummary:
    delta_max: int
    entries: tuple[ThresholdFormulaEntry, ...]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def threshold_formula_entries(delta_max: int = DELTA_MAX) -> tuple[ThresholdFormulaEntry, ...]:
    if delta_max != 4:
        raise ValueError("The current threshold-formula audit is only fixed at Delta_max = 4.")

    mixed = mixed_witness_summary(delta_max)
    rank2, bare_scalar, derivative_only = mixed.entries

    _require(rank2.self_weight == 2 and rank2.mixed_weight == 2 and rank2.w_min == 2, "Rank-2 current-case values drifted.")
    _require(bare_scalar.self_weight == 1 and bare_scalar.mixed_weight == 3 and bare_scalar.w_min == 1, "Bare-scalar current-case values drifted.")
    _require(derivative_only.self_weight == 4 and derivative_only.mixed_weight == 4 and derivative_only.w_min == 4, "Derivative-only current-case values drifted.")

    return (
        ThresholdFormulaEntry(
            family_class="rank2_stf",
            self_formula="2*w_X",
            mixed_formula="w_X + 1",
            w_min_formula="min(2*w_X, w_X + 1)",
            current_case_value="w_X=1 -> (W_self, W_mix, W_min) = (2, 2, 2)",
            sharpness_status="mixed-aware sharp; current unsuppressed case is tied",
            threshold_type="mixed-aware",
        ),
        ThresholdFormulaEntry(
            family_class="rank0_scalar_unsuppressed",
            self_formula="w_S",
            mixed_formula="w_S + 2",
            w_min_formula="w_S",
            current_case_value="w_S=1 -> (W_self, W_mix, W_min) = (1, 3, 1)",
            sharpness_status="self-only sharp",
            threshold_type="self-only",
        ),
        ThresholdFormulaEntry(
            family_class="rank0_scalar_derivative_only",
            self_formula="2*w_D",
            mixed_formula="w_D + 2",
            w_min_formula="min(2*w_D, w_D + 2)",
            current_case_value="w_D=2 -> (W_self, W_mix, W_min) = (4, 4, 4)",
            sharpness_status="tied-sharp",
            threshold_type="tied-sharp",
        ),
    )


def threshold_formula_summary(delta_max: int = DELTA_MAX) -> ThresholdFormulaSummary:
    return ThresholdFormulaSummary(delta_max=delta_max, entries=threshold_formula_entries(delta_max))


def threshold_formula_report(delta_max: int = DELTA_MAX) -> str:
    summary = threshold_formula_summary(delta_max)
    lines = [
        "family_class\tself_formula\tmixed_formula\tW_min_formula\tcurrent_case_value\tsharpness_status",
    ]
    for entry in summary.entries:
        lines.append(
            "\t".join(
                (
                    entry.family_class,
                    entry.self_formula,
                    entry.mixed_formula,
                    entry.w_min_formula,
                    entry.current_case_value,
                    entry.sharpness_status,
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(threshold_formula_report())


if __name__ == "__main__":
    main()
