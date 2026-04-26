from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math


DELTA_MAX = 4
WINDOW = 1.0
GRID_SIZE = 4000
M0 = 1.0
ALPHA = 0.5
BETA = 0.25
LAMBDA_EXP = 1.0
LAMBDA_GAMMA = 1.0
T_H = 1.0
GAMMA = 0.5


def _history_constant(tau: float) -> float:
    if -WINDOW <= tau <= 0.0:
        return 1.0
    return 0.0


def _history_oscillatory(tau: float) -> float:
    if -WINDOW <= tau <= 0.0:
        return 1.0 + math.sin(2.0 * math.pi * tau)
    return 0.0


def _midpoint_integral(kernel, history, *, upper: float = WINDOW, steps: int = GRID_SIZE) -> float:
    width = upper / steps
    total = 0.0
    for index in range(steps):
        sigma = (index + 0.5) * width
        total += kernel(sigma) * history(-sigma)
    return total * width


def _kernel_exponential(sigma: float) -> float:
    return math.exp(-sigma / T_H) / T_H


def _kernel_power_law(sigma: float) -> float:
    scaled = sigma / T_H
    return scaled ** (-GAMMA) / (math.gamma(1.0 - GAMMA) * T_H)


def _local_analytic_response(y_now: float) -> float:
    return M0 + ALPHA * y_now + BETA * y_now * y_now


def _memory_response(history, kernel, coupling: float) -> float:
    y_now = history(0.0)
    return M0 + ALPHA * y_now + coupling * _midpoint_integral(kernel, history)


@dataclass(frozen=True)
class HereditaryKernelCase:
    case_id: str
    model_class: str
    formula: str
    locality_kept: bool
    analyticity_kept: bool
    finite_primitive_family_envelope_kept: bool
    finite_local_jet_valid: bool
    finite_state_markovianizable: bool
    first_exact_failure_mode: str | None
    generalized_replacement_data: str
    theorem_layer_first_broken: str | None
    canonical_counterexample: bool
    current_coordinate_value: float
    same_point_response_difference: float


@dataclass(frozen=True)
class HereditaryKernelSummary:
    delta_max: int
    finite_primitive_family_envelope_kept_in_all_cases: bool
    analyticity_kept_in_all_cases: bool
    smallest_genuinely_hereditary_counterexample: str
    broken_layer: str
    cases: tuple[HereditaryKernelCase, ...]


def hereditary_kernel_cases() -> tuple[HereditaryKernelCase, ...]:
    y_now = _history_constant(0.0)

    local_a = _local_analytic_response(y_now)
    local_b = _local_analytic_response(_history_oscillatory(0.0))

    exp_a = _memory_response(_history_constant, _kernel_exponential, LAMBDA_EXP)
    exp_b = _memory_response(_history_oscillatory, _kernel_exponential, LAMBDA_EXP)

    pow_a = _memory_response(_history_constant, _kernel_power_law, LAMBDA_GAMMA)
    pow_b = _memory_response(_history_oscillatory, _kernel_power_law, LAMBDA_GAMMA)

    return (
        HereditaryKernelCase(
            case_id="local_analytic_control",
            model_class="local analytic control",
            formula="m_A(Y) = m0 + alpha*Y + beta*Y^2",
            locality_kept=True,
            analyticity_kept=True,
            finite_primitive_family_envelope_kept=True,
            finite_local_jet_valid=True,
            finite_state_markovianizable=False,
            first_exact_failure_mode=None,
            generalized_replacement_data="ordinary finite Taylor coefficients",
            theorem_layer_first_broken=None,
            canonical_counterexample=False,
            current_coordinate_value=y_now,
            same_point_response_difference=abs(local_a - local_b),
        ),
        HereditaryKernelCase(
            case_id="exponential_memory_control",
            model_class="exponential memory control",
            formula=(
                "m_A[Y](tau) = m0 + alpha*Y(tau) + lambda_exp*int exp(-(tau-tau')/T_h)"
                " Y(tau') d tau' / T_h"
            ),
            locality_kept=False,
            analyticity_kept=True,
            finite_primitive_family_envelope_kept=True,
            finite_local_jet_valid=False,
            finite_state_markovianizable=True,
            first_exact_failure_mode=(
                "two histories with the same Y(tau0) give different responses, so no local jet in Y(tau) alone exists"
            ),
            generalized_replacement_data="finite auxiliary-state data (m0, alpha, lambda_exp, T_h, chi_A)",
            theorem_layer_first_broken=(
                "local monopole reduction to m_A(Y(tau)) alone fails, but the kernel collapses into a finite-state A4-type extension"
            ),
            canonical_counterexample=False,
            current_coordinate_value=y_now,
            same_point_response_difference=abs(exp_a - exp_b),
        ),
        HereditaryKernelCase(
            case_id="power_law_single_coordinate",
            model_class="genuinely hereditary power-law kernel",
            formula=(
                "m_A[Y](tau) = m0 + alpha*Y(tau) + lambda_gamma*int K_gamma(tau-tau')"
                " Y(tau') d tau', K_gamma(s) ~ Theta(s) s^(-gamma)"
            ),
            locality_kept=False,
            analyticity_kept=True,
            finite_primitive_family_envelope_kept=True,
            finite_local_jet_valid=False,
            finite_state_markovianizable=False,
            first_exact_failure_mode=(
                "the causal power-law kernel distinguishes histories with the same instantaneous Y(tau0) and has no finite-state realization"
            ),
            generalized_replacement_data=(
                "kernel branch data (m0, alpha, lambda_gamma, gamma, T_h) or a continuous spectral density"
            ),
            theorem_layer_first_broken=(
                "local monopole reduction to a function of instantaneous normal-form coordinates (A3; Lemmas 55 and 56 become inapplicable)"
            ),
            canonical_counterexample=True,
            current_coordinate_value=y_now,
            same_point_response_difference=abs(pow_a - pow_b),
        ),
    )


def hereditary_kernel_summary(delta_max: int = DELTA_MAX) -> HereditaryKernelSummary:
    if delta_max != 4:
        raise ValueError("The hereditary kernel demo is only fixed for Delta <= 4.")
    cases = hereditary_kernel_cases()
    return HereditaryKernelSummary(
        delta_max=delta_max,
        finite_primitive_family_envelope_kept_in_all_cases=all(
            case.finite_primitive_family_envelope_kept for case in cases
        ),
        analyticity_kept_in_all_cases=all(case.analyticity_kept for case in cases),
        smallest_genuinely_hereditary_counterexample="power_law_single_coordinate",
        broken_layer=(
            "local monopole reduction to a function of instantaneous normal-form coordinates (A3; Lemmas 55 and 56 become inapplicable)"
        ),
        cases=cases,
    )


def hereditary_kernel_report(delta_max: int = DELTA_MAX) -> str:
    summary = hereditary_kernel_summary(delta_max=delta_max)
    return json.dumps(
        {
            "delta_max": summary.delta_max,
            "finite_primitive_family_envelope_kept_in_all_cases": summary.finite_primitive_family_envelope_kept_in_all_cases,
            "analyticity_kept_in_all_cases": summary.analyticity_kept_in_all_cases,
            "smallest_genuinely_hereditary_counterexample": summary.smallest_genuinely_hereditary_counterexample,
            "broken_layer": summary.broken_layer,
            "cases": [asdict(case) for case in summary.cases],
        },
        indent=2,
        sort_keys=True,
    )


def main() -> None:
    print(hereditary_kernel_report())


if __name__ == "__main__":
    main()
