from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math


DELTA_MAX = 4


def _flat_activation(y: float) -> float:
    if y <= 0.0:
        return 0.0
    return math.exp(-1.0 / (y * y))


def _threshold_sqrt_activation(y: float, y_c: float = 0.0) -> float:
    if y <= y_c:
        return 0.0
    return math.sqrt(y - y_c)


@dataclass(frozen=True)
class NonanalyticJetCase:
    case_id: str
    model_class: str
    formula: str
    locality_kept: bool
    finite_family_operator_closure_kept: bool
    finite_taylor_jet_valid: bool
    first_exact_failure_mode: str | None
    generalized_replacement_data: str
    theorem_layer_broken: str | None
    canonical_counterexample: bool
    sample_point: float
    model_value_at_sample: float
    jet_value_at_sample: float | None


@dataclass(frozen=True)
class NonanalyticJetSummary:
    delta_max: int
    locality_kept_in_all_cases: bool
    finite_family_operator_closure_kept_in_all_cases: bool
    smallest_local_nonanalytic_counterexample: str
    broken_layer: str
    cases: tuple[NonanalyticJetCase, ...]


def nonanalytic_jet_cases() -> tuple[NonanalyticJetCase, ...]:
    analytic_sample = 0.5
    analytic_value = 1.0 + 2.0 * analytic_sample + 3.0 * analytic_sample**2

    flat_sample = 0.5
    flat_value = 1.0 + _flat_activation(flat_sample)

    threshold_sample = 0.25
    threshold_value = 1.0 + _threshold_sqrt_activation(threshold_sample)

    return (
        NonanalyticJetCase(
            case_id="analytic_control_quadratic",
            model_class="analytic control",
            formula="m_A(Y) = m0 + alpha*Y + beta*Y^2",
            locality_kept=True,
            finite_family_operator_closure_kept=True,
            finite_taylor_jet_valid=True,
            first_exact_failure_mode=None,
            generalized_replacement_data="ordinary finite Taylor coefficients",
            theorem_layer_broken=None,
            canonical_counterexample=False,
            sample_point=analytic_sample,
            model_value_at_sample=analytic_value,
            jet_value_at_sample=analytic_value,
        ),
        NonanalyticJetCase(
            case_id="smooth_flat_single_coordinate",
            model_class="smooth flat nonanalytic",
            formula="m_A(Y) = m0 + alpha*exp(-1/Y^2)*Theta(Y)",
            locality_kept=True,
            finite_family_operator_closure_kept=True,
            finite_taylor_jet_valid=False,
            first_exact_failure_mode=(
                "all Taylor coefficients about Y=0 vanish, but the response is nonzero for Y>0"
            ),
            generalized_replacement_data="non-Taylor monopole germ data (m0, alpha, flat-profile label)",
            theorem_layer_broken="analytic monopole jet collapse (Lemma 55 / A5)",
            canonical_counterexample=True,
            sample_point=flat_sample,
            model_value_at_sample=flat_value,
            jet_value_at_sample=1.0,
        ),
        NonanalyticJetCase(
            case_id="threshold_sqrt_activation",
            model_class="threshold branch point",
            formula="m_A(Y) = m0 + alpha*Theta(Y-Yc)*sqrt(Y-Yc)",
            locality_kept=True,
            finite_family_operator_closure_kept=True,
            finite_taylor_jet_valid=False,
            first_exact_failure_mode=(
                "branch point at Y=Yc prevents a valid analytic Taylor jet at the activation point"
            ),
            generalized_replacement_data="threshold parameters (m0, alpha, Yc, branch exponent, branch label)",
            theorem_layer_broken="analytic monopole jet collapse (Lemma 55 / A5)",
            canonical_counterexample=False,
            sample_point=threshold_sample,
            model_value_at_sample=threshold_value,
            jet_value_at_sample=None,
        ),
    )


def nonanalytic_jet_summary(delta_max: int = DELTA_MAX) -> NonanalyticJetSummary:
    if delta_max != 4:
        raise ValueError("The nonanalytic jet demo is only fixed for Delta <= 4.")
    cases = nonanalytic_jet_cases()
    return NonanalyticJetSummary(
        delta_max=delta_max,
        locality_kept_in_all_cases=all(case.locality_kept for case in cases),
        finite_family_operator_closure_kept_in_all_cases=all(
            case.finite_family_operator_closure_kept for case in cases
        ),
        smallest_local_nonanalytic_counterexample="smooth_flat_single_coordinate",
        broken_layer="analytic monopole jet collapse (Lemma 55 / A5)",
        cases=cases,
    )


def nonanalytic_jet_report(delta_max: int = DELTA_MAX) -> str:
    summary = nonanalytic_jet_summary(delta_max=delta_max)
    return json.dumps(
        {
            "delta_max": summary.delta_max,
            "locality_kept_in_all_cases": summary.locality_kept_in_all_cases,
            "finite_family_operator_closure_kept_in_all_cases": summary.finite_family_operator_closure_kept_in_all_cases,
            "smallest_local_nonanalytic_counterexample": summary.smallest_local_nonanalytic_counterexample,
            "broken_layer": summary.broken_layer,
            "cases": [asdict(case) for case in summary.cases],
        },
        indent=2,
        sort_keys=True,
    )


def main() -> None:
    print(nonanalytic_jet_report())


if __name__ == "__main__":
    main()
