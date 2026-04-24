from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

SYMBOLIC_ROOT = Path(__file__).resolve().parents[1]
if str(SYMBOLIC_ROOT) not in sys.path:
    sys.path.insert(0, str(SYMBOLIC_ROOT))

from chi_relaxation_response import (  # noqa: E402
    relaxation_residual,
    small_x_series,
    steady_state_response,
    transfer_function,
)
from chi_two_frequency_response import (  # noqa: E402
    linear_sideband_residuals,
    relaxation_residual as two_frequency_residual,
    two_frequency_response,
)


def test_monochromatic_response_solves_relaxation_equation() -> None:
    response = steady_state_response()
    residual = relaxation_residual(response)
    assert sp.simplify(residual) == 0


def test_transfer_function_components_match_real_solution() -> None:
    response = steady_state_response()
    omega, tau, alpha, f0 = response.omega, response.tau, response.alpha, response.f0
    x = omega * tau
    h = transfer_function(alpha=alpha, omega=omega, tau=tau)

    assert sp.simplify(sp.re(h) * f0 - response.in_phase) == 0
    assert sp.simplify(-sp.im(h) * f0 - response.quadrature) == 0
    assert sp.simplify(response.quadrature / response.in_phase - x) == 0


def test_small_frequency_limit() -> None:
    series = small_x_series(order=5)
    x, alpha, f0 = sp.symbols("x alpha F0", positive=True, real=True)

    assert sp.expand(series["in_phase"] - alpha * f0 * (1 - x**2 + x**4)) == 0
    assert sp.expand(series["quadrature"] - alpha * f0 * (x - x**3)) == 0


def test_static_limits_remove_quadrature() -> None:
    response = steady_state_response()
    assert response.quadrature.subs(response.tau, 0) == 0
    assert response.quadrature.subs(response.omega, 0) == 0


def test_large_frequency_asymptotics() -> None:
    x = sp.symbols("x", positive=True, real=True)
    alpha, f0 = sp.symbols("alpha F0", positive=True, real=True)
    in_phase = alpha * f0 / (1 + x**2)
    quadrature = alpha * f0 * x / (1 + x**2)

    assert sp.limit(x**2 * in_phase, x, sp.oo) == alpha * f0
    assert sp.limit(x * quadrature, x, sp.oo) == alpha * f0


def test_two_frequency_response_has_only_input_frequencies() -> None:
    response = two_frequency_response()
    assert set(response.components) == {response.omega1, response.omega2}


def test_two_frequency_response_solves_relaxation_equation() -> None:
    response = two_frequency_response()
    residual = two_frequency_residual(response)
    assert sp.simplify(residual) == 0


def test_linear_two_frequency_model_has_no_sidebands() -> None:
    residuals = linear_sideband_residuals()
    assert residuals["sum"] == 0
    assert residuals["difference"] == 0


def main() -> None:
    test_monochromatic_response_solves_relaxation_equation()
    test_transfer_function_components_match_real_solution()
    test_small_frequency_limit()
    test_static_limits_remove_quadrature()
    test_large_frequency_asymptotics()
    test_two_frequency_response_has_only_input_frequencies()
    test_two_frequency_response_solves_relaxation_equation()
    test_linear_two_frequency_model_has_no_sidebands()
    print("dynamic chi symbolic checks passed")


if __name__ == "__main__":
    main()
