from __future__ import annotations

from dataclasses import asdict, dataclass
import json


DELTA_MAX = 4


@dataclass(frozen=True)
class WeightSpectrumCase:
    case_id: str
    model_class: str
    total_family_count: str
    finite_primitive_family_content: bool
    locally_finite_below_delta_max: bool
    candidate_operator_space_finite: bool
    normal_form_quotient_finite: bool
    theorem_layer_first_broken: str | None
    representative_data: str


@dataclass(frozen=True)
class WeightSpectrumSummary:
    delta_max: int
    local_weight_spectrum_finiteness_suffices: bool
    a8_stronger_than_necessary: bool
    sharp_a8_failure_mode: str
    cases: tuple[WeightSpectrumCase, ...]


def weight_spectrum_cases() -> tuple[WeightSpectrumCase, ...]:
    return (
        WeightSpectrumCase(
            case_id="finite_family_catalog_control",
            model_class="finite catalog control",
            total_family_count="finite",
            finite_primitive_family_content=True,
            locally_finite_below_delta_max=True,
            candidate_operator_space_finite=True,
            normal_form_quotient_finite=True,
            theorem_layer_first_broken=None,
            representative_data=(
                "Representative audited scalar/vector/STF catalog with finitely many admitted species overall."
            ),
        ),
        WeightSpectrumCase(
            case_id="infinite_but_locally_finite_weight_spectrum",
            model_class="infinite but locally finite weight spectrum",
            total_family_count="infinite",
            finite_primitive_family_content=False,
            locally_finite_below_delta_max=True,
            candidate_operator_space_finite=True,
            normal_form_quotient_finite=True,
            theorem_layer_first_broken=None,
            representative_data=(
                "Infinite species tower {Y^(n)} with intrinsic weights w_n = n+1, so only finitely many satisfy w_n <= 4."
            ),
        ),
        WeightSpectrumCase(
            case_id="infinite_low_weight_stf_tower",
            model_class="infinite low-weight STF tower",
            total_family_count="infinite",
            finite_primitive_family_content=False,
            locally_finite_below_delta_max=False,
            candidate_operator_space_finite=False,
            normal_form_quotient_finite=False,
            theorem_layer_first_broken="Lemma 53 fixed-order candidate operator finiteness",
            representative_data=(
                "Infinite genuine STF species tower {T^(n)} with identical intrinsic weight 1, giving infinitely many quadratic self witnesses (T^(n))^2 at weight 2."
            ),
        ),
    )


def weight_spectrum_summary(delta_max: int = DELTA_MAX) -> WeightSpectrumSummary:
    if delta_max != 4:
        raise ValueError("The weight-spectrum demo is only fixed for Delta <= 4.")
    return WeightSpectrumSummary(
        delta_max=delta_max,
        local_weight_spectrum_finiteness_suffices=True,
        a8_stronger_than_necessary=True,
        sharp_a8_failure_mode="infinite_low_weight_stf_tower",
        cases=weight_spectrum_cases(),
    )


def weight_spectrum_report(delta_max: int = DELTA_MAX) -> str:
    summary = weight_spectrum_summary(delta_max=delta_max)
    return json.dumps(
        {
            "delta_max": summary.delta_max,
            "local_weight_spectrum_finiteness_suffices": summary.local_weight_spectrum_finiteness_suffices,
            "a8_stronger_than_necessary": summary.a8_stronger_than_necessary,
            "sharp_a8_failure_mode": summary.sharp_a8_failure_mode,
            "cases": [asdict(case) for case in summary.cases],
        },
        indent=2,
        sort_keys=True,
    )


def main() -> None:
    print(weight_spectrum_report())


if __name__ == "__main__":
    main()
