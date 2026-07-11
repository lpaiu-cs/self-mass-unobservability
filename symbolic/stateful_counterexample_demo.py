from __future__ import annotations

from dataclasses import asdict, dataclass
import json


DELTA_MAX = 4
M0 = 1.0
ALPHA = 0.5
BETA = 0.4
LAMBDA = 0.75


def _monopole_no_state(y_value: float) -> float:
    return M0 + ALPHA * y_value + 0.25 * y_value * y_value


def _monopole_stateful(y_value: float, chi_value: float) -> float:
    return M0 + ALPHA * y_value + LAMBDA * chi_value


@dataclass(frozen=True)
class StatefulCounterexampleCase:
    case_id: str
    model_class: str
    formula: str
    locality_kept: bool
    analyticity_kept: bool
    finite_primitive_family_envelope_kept: bool
    y_only_finite_jet_valid: bool
    finite_state_augmented_bookkeeping_valid: bool
    first_exact_failure_mode: str | None
    bookkeeping_replacement_data: str
    theorem_layer_first_broken: str | None
    canonical_counterexample: bool
    current_coordinate_value: float
    same_point_response_difference: float


@dataclass(frozen=True)
class StatefulCounterexampleSummary:
    delta_max: int
    locality_kept_in_all_cases: bool
    analyticity_kept_in_all_cases: bool
    finite_primitive_family_envelope_kept_in_all_cases: bool
    smallest_local_finite_state_counterexample: str
    broken_layer: str
    finite_state_augmented_collapse_survives: bool
    cases: tuple[StatefulCounterexampleCase, ...]


def stateful_counterexample_cases() -> tuple[StatefulCounterexampleCase, ...]:
    y_now = 1.0
    no_state_a = _monopole_no_state(y_now)
    no_state_b = _monopole_no_state(y_now)

    chi_slaved = BETA * y_now
    slaved_a = _monopole_stateful(y_now, chi_slaved)
    slaved_b = _monopole_stateful(y_now, chi_slaved)

    chi_dynamic_a = BETA * y_now
    chi_dynamic_b = 0.0
    dynamic_a = _monopole_stateful(y_now, chi_dynamic_a)
    dynamic_b = _monopole_stateful(y_now, chi_dynamic_b)

    return (
        StatefulCounterexampleCase(
            case_id="local_analytic_no_state_control",
            model_class="local analytic no-state control",
            formula="m_A(Y) = m0 + alpha*Y + beta*Y^2",
            locality_kept=True,
            analyticity_kept=True,
            finite_primitive_family_envelope_kept=True,
            y_only_finite_jet_valid=True,
            finite_state_augmented_bookkeeping_valid=True,
            first_exact_failure_mode=None,
            bookkeeping_replacement_data="ordinary finite Taylor sensitivities in Y^I",
            theorem_layer_first_broken=None,
            canonical_counterexample=False,
            current_coordinate_value=y_now,
            same_point_response_difference=abs(no_state_a - no_state_b),
        ),
        StatefulCounterexampleCase(
            case_id="adiabatic_slaved_local_state",
            model_class="adiabatic/slaved local-state control",
            formula="m_A(Y, chi) = m0 + alpha*Y + lambda*chi with chi = beta*Y",
            locality_kept=True,
            analyticity_kept=True,
            finite_primitive_family_envelope_kept=True,
            y_only_finite_jet_valid=True,
            finite_state_augmented_bookkeeping_valid=True,
            first_exact_failure_mode=None,
            bookkeeping_replacement_data=(
                "Y-only sensitivities after eliminating chi, or equivalently finite augmented coordinates with a slaved state"
            ),
            theorem_layer_first_broken=None,
            canonical_counterexample=False,
            current_coordinate_value=y_now,
            same_point_response_difference=abs(slaved_a - slaved_b),
        ),
        StatefulCounterexampleCase(
            case_id="dynamical_one_state_chi",
            model_class="genuinely dynamical one-state local model",
            formula=(
                "m_A(Y, chi) = m0 + alpha*Y + lambda*chi, "
                "dot chi = -(1/T_h)*chi + (beta/T_h)*Y"
            ),
            locality_kept=True,
            analyticity_kept=True,
            finite_primitive_family_envelope_kept=True,
            y_only_finite_jet_valid=False,
            finite_state_augmented_bookkeeping_valid=True,
            first_exact_failure_mode=(
                "the same instantaneous Y with different orbital-timescale state values chi gives different monopole responses"
            ),
            bookkeeping_replacement_data=(
                "finite augmented state-space data (Y^I, chi^a) plus local state-evolution parameters, kept separate from Wilson coefficients"
            ),
            theorem_layer_first_broken=(
                "no-state reduction and Lemma 55 in Y-only form fail; the original Lemma 56 sensitivity split also fails in Y-only form"
            ),
            canonical_counterexample=True,
            current_coordinate_value=y_now,
            same_point_response_difference=abs(dynamic_a - dynamic_b),
        ),
    )


def stateful_counterexample_summary(delta_max: int = DELTA_MAX) -> StatefulCounterexampleSummary:
    if delta_max != 4:
        raise ValueError("The stateful counterexample demo is only fixed for Delta <= 4.")
    cases = stateful_counterexample_cases()
    return StatefulCounterexampleSummary(
        delta_max=delta_max,
        locality_kept_in_all_cases=all(case.locality_kept for case in cases),
        analyticity_kept_in_all_cases=all(case.analyticity_kept for case in cases),
        finite_primitive_family_envelope_kept_in_all_cases=all(
            case.finite_primitive_family_envelope_kept for case in cases
        ),
        smallest_local_finite_state_counterexample="dynamical_one_state_chi",
        broken_layer="no-state reduction to m_A(Y) alone (A4; Lemma 55 in Y-only form fails)",
        finite_state_augmented_collapse_survives=True,
        cases=cases,
    )


def stateful_counterexample_report(delta_max: int = DELTA_MAX) -> str:
    summary = stateful_counterexample_summary(delta_max=delta_max)
    return json.dumps(
        {
            "delta_max": summary.delta_max,
            "locality_kept_in_all_cases": summary.locality_kept_in_all_cases,
            "analyticity_kept_in_all_cases": summary.analyticity_kept_in_all_cases,
            "finite_primitive_family_envelope_kept_in_all_cases": summary.finite_primitive_family_envelope_kept_in_all_cases,
            "smallest_local_finite_state_counterexample": summary.smallest_local_finite_state_counterexample,
            "broken_layer": summary.broken_layer,
            "finite_state_augmented_collapse_survives": summary.finite_state_augmented_collapse_survives,
            "cases": [asdict(case) for case in summary.cases],
        },
        indent=2,
        sort_keys=True,
    )


def main() -> None:
    print(stateful_counterexample_report())


if __name__ == "__main__":
    main()
