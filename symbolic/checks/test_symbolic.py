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
from chi_basis_audit import (  # noqa: E402
    derivative_series_residual,
    polynomial_obstruction,
    rows_as_dicts,
    single_frequency_derivative_fit,
    static_only_residual,
    symbols as basis_symbols,
    validate_rows,
)
from chi_two_frequency_response import (  # noqa: E402
    linear_sideband_residuals,
    relaxation_residual as two_frequency_residual,
    two_frequency_response,
)
from frequency_sweep_distinguishability import (  # noqa: E402
    payload as frequency_sweep_payload,
    interpolation_residual,
    linear_interpolation_residual_example,
    low_frequency_error_bound,
    polynomial_sample_limit,
    real_coefficient_sample_limit,
    real_odd_channel_residual,
    real_odd_symbols,
    rows_as_dicts as frequency_sweep_rows,
    symbols as sweep_symbols,
    taylor_residual,
    taylor_residual_closed_form,
    validate_rows as validate_frequency_sweep_rows,
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


def test_basis_audit_rows_are_classified() -> None:
    rows = rows_as_dicts()
    validate_rows(rows)
    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["observable_status"] for row in rows)


def test_static_basis_leaves_quadrature_except_collapse_limits() -> None:
    sym = basis_symbols()
    residual = static_only_residual()["unmatched_quadrature"]

    assert residual != 0
    assert sp.simplify(residual.subs(sym["tau_chi"], 0)) == 0
    assert sp.simplify(residual.subs(sym["Omega"], 0)) == 0
    assert sp.simplify(residual.subs(sym["alpha"], 0)) == 0
    assert sp.simplify(residual.subs(sym["c_chi"], 0)) == 0


def test_single_frequency_dotF_basis_fits_both_quadratures() -> None:
    sym = basis_symbols()
    omega = sym["Omega"]
    f0 = sym["F0"]
    fit = single_frequency_derivative_fit()
    residual = static_only_residual()

    fitted_in_phase = sp.simplify(fit["a0"] * f0)
    fitted_quadrature = sp.simplify(-fit["a1"] * omega * f0)

    assert sp.simplify(fitted_in_phase - residual["matched_static_amplitude"]) == 0
    assert sp.simplify(fitted_quadrature - residual["unmatched_quadrature"]) == 0


def test_derivative_series_residual_starts_after_truncation_order() -> None:
    sym = basis_symbols()
    omega = sym["Omega"]
    residual = derivative_series_residual(3).subs(
        {
            sym["c_Y"]: 0,
            sym["c_chi"]: 1,
            sym["alpha"]: 1,
            sym["tau_chi"]: 1,
        }
    )

    assert sp.series(residual, omega, 0, 4).removeO() == 0


def test_finite_polynomial_obstruction_records_pole_mismatch() -> None:
    obstruction = polynomial_obstruction(3)
    assert "d3" in obstruction["highest_power_coefficient"]
    assert obstruction["constant_coefficient_after_recursion"] == "-alpha"


def test_frequency_sweep_sample_count_boundary() -> None:
    boundary = polynomial_sample_limit(4)
    assert boundary["exact_interpolation_limit"] == 5
    assert boundary["first_exact_obstruction_count"] == 6


def test_real_coefficient_boundary_is_sharper() -> None:
    boundary = real_coefficient_sample_limit(4)
    assert boundary["positive_frequency_exact_interpolation_limit"] == 2
    assert boundary["first_positive_frequency_obstruction_count"] == 3


def test_frequency_sweep_residual_has_expected_collapse_limits() -> None:
    sym = sweep_symbols(2)
    residual = interpolation_residual(2)

    assert residual != 0
    assert sp.simplify(residual.subs(sym.beta, 0)) == 0
    assert sp.simplify(residual.subs(sym.tau, 0)) == 0
    assert sp.simplify(residual.subs(sym.z_star, sym.z_nodes[0])) == 0


def test_real_odd_channel_residual_has_expected_collapse_limits() -> None:
    sym = real_odd_symbols(3)
    residual = real_odd_channel_residual(3)

    assert residual != 0
    assert sp.simplify(residual.subs(sym.beta, 0)) == 0
    assert sp.simplify(residual.subs(sym.tau, 0)) == 0
    assert sp.simplify(residual.subs(sym.u_star, sym.u_nodes[0])) == 0


def test_frequency_sweep_residual_matches_direct_interpolation() -> None:
    example = linear_interpolation_residual_example()
    assert sp.simplify(example["direct_residual"] - example["formula_residual"]) == 0
    assert example["direct_residual"] != 0


def test_taylor_residual_closed_form() -> None:
    assert sp.simplify(taylor_residual(5) - taylor_residual_closed_form(5)) == 0


def test_low_frequency_error_bound_power() -> None:
    rho, beta_abs = sp.symbols("rho abs_beta", nonnegative=True, real=True)
    assert low_frequency_error_bound(3) == beta_abs * rho**4


def test_multifrequency_audit_rows_are_classified() -> None:
    rows = frequency_sweep_rows()
    validate_frequency_sweep_rows(rows)
    assert rows
    assert all(row["N_frequencies"] for row in rows)
    assert all(row["comparator_basis"] for row in rows)
    assert all(row["verdict"] for row in rows)
    assert all(row["surviving_target"] for row in rows)


def test_multifrequency_audit_records_minimality_boundary() -> None:
    rows = {row["term"]: row for row in frequency_sweep_rows()}

    assert rows["single_frequency_dotF"]["verdict"] == "degenerate"
    assert (
        rows["real_shared_coefficients_first_obstruction"]["N_frequencies"]
        == "floor((N+1)/2)+1"
    )
    assert (
        rows["complex_shared_coefficients_first_obstruction"]["N_frequencies"]
        == "N+2"
    )


def test_multifrequency_payload_has_required_table_columns() -> None:
    data = frequency_sweep_payload()
    schema = set(data["schema"])

    assert "N_frequencies" in schema
    assert "comparator_basis" in schema
    assert "verdict" in schema
    assert "surviving_target" in schema


def main() -> None:
    test_monochromatic_response_solves_relaxation_equation()
    test_transfer_function_components_match_real_solution()
    test_small_frequency_limit()
    test_static_limits_remove_quadrature()
    test_large_frequency_asymptotics()
    test_two_frequency_response_has_only_input_frequencies()
    test_two_frequency_response_solves_relaxation_equation()
    test_linear_two_frequency_model_has_no_sidebands()
    test_basis_audit_rows_are_classified()
    test_static_basis_leaves_quadrature_except_collapse_limits()
    test_single_frequency_dotF_basis_fits_both_quadratures()
    test_derivative_series_residual_starts_after_truncation_order()
    test_finite_polynomial_obstruction_records_pole_mismatch()
    test_frequency_sweep_sample_count_boundary()
    test_real_coefficient_boundary_is_sharper()
    test_frequency_sweep_residual_has_expected_collapse_limits()
    test_real_odd_channel_residual_has_expected_collapse_limits()
    test_frequency_sweep_residual_matches_direct_interpolation()
    test_taylor_residual_closed_form()
    test_low_frequency_error_bound_power()
    test_multifrequency_audit_rows_are_classified()
    test_multifrequency_audit_records_minimality_boundary()
    test_multifrequency_payload_has_required_table_columns()
    print("dynamic chi symbolic checks passed")


if __name__ == "__main__":
    main()
