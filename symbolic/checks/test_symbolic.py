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
from forcing_observable_dictionary import (  # noqa: E402
    observable_harmonic_components,
    orbital_two_harmonic_dotf_obstruction,
    payload as forcing_dictionary_payload,
    power_law_orbital_harmonics_to_e2,
    power_law_series_residual_to_e2,
    rows_as_dicts as forcing_dictionary_rows,
    symbols as forcing_symbols,
    validate_rows as validate_forcing_dictionary_rows,
)
from nonlinear_sideband_test import (  # noqa: E402
    linear_projection_sideband_output,
    nonlinear_drive_chi_sidebands,
    nonlinear_drive_rhs_sidebands,
    nonlinear_readout_difference_sideband,
    nonlinear_readout_sum_sideband,
    orbital_3n_nonlinear_drive_amplitude,
    orbital_n_2n_input_amplitudes,
    payload as nonlinear_sideband_payload,
    rows_as_dicts as nonlinear_sideband_rows,
    symbols as nonlinear_sideband_symbols,
    validate_rows as validate_nonlinear_sideband_rows,
)
from projection_channel_audit import (  # noqa: E402
    arbitrary_projection_solution,
    observed_range_transfer,
    payload as projection_channel_payload,
    projection_pole_condition,
    range_deprojection_residual,
    range_equation_residual,
    range_projection,
    range_relaxation_pole_residue,
    relaxation_transfer as projection_relaxation_transfer,
    rows_as_dicts as projection_channel_rows,
    symbols as projection_symbols,
    validate_rows as validate_projection_channel_rows,
)
from triple_shared_tau_bridge import (  # noqa: E402
    complex_two_frequency_verdict,
    payload as triple_bridge_payload,
    real_two_frequency_verdict,
    rows_as_dicts as triple_bridge_rows,
    transfer_at_carriers,
    two_carrier_real_odd_residual_order1,
    validate_rows as validate_triple_bridge_rows,
)
from triple_gr_carrier_inventory import (  # noqa: E402
    carrier_inventory as triple_gr_carrier_inventory,
    complex_k_frequency_verdict as complex_k_carrier_verdict,
    payload as triple_gr_inventory_payload,
    real_k_frequency_verdict as real_k_carrier_verdict,
    rows_as_dicts as triple_gr_inventory_rows,
    three_carrier_conditions,
    validate_rows as validate_triple_gr_inventory_rows,
)
from triple_projection_nuisance_gate import (  # noqa: E402
    arbitrary_projection_reconstruction_residual as gate_arbitrary_projection_reconstruction_residual,
    calibrated_deprojection_residual as gate_calibrated_deprojection_residual,
    payload as triple_projection_gate_payload,
    projection_model_forms,
    range_deprojection_residual as gate_range_deprojection_residual,
    rows_as_dicts as triple_projection_gate_rows,
    triple_projection_mapping,
    validate_rows as validate_triple_projection_gate_rows,
)
from triple_projection_manifold_gate import (  # noqa: E402
    payload as triple_projection_manifold_payload,
    phase_lock_residual,
    rank_summary as triple_projection_rank_summary,
    rows_as_dicts as triple_projection_manifold_rows,
    validate_rows as validate_triple_projection_manifold_rows,
)
from named_timing_model_projection_audit import (  # noqa: E402
    payload as named_timing_model_payload,
    rows_as_dicts as named_timing_model_rows,
    source_release_ledger as named_timing_source_release_ledger,
    validate_rows as validate_named_timing_model_rows,
    verdict_summary as named_timing_verdict_summary,
)
from nutimo_runtime_worthiness_pilot import (  # noqa: E402
    jacobian_rank_gate_verdict,
    minimal_external_artifacts as nutimo_pilot_minimal_external_artifacts,
    payload as nutimo_runtime_pilot_payload,
    pilot_stage_order as nutimo_pilot_stage_order,
    rows_as_dicts as nutimo_runtime_pilot_rows,
    validate_rows as validate_nutimo_runtime_pilot_rows,
)
from nutimo_external_handoff_packet import (  # noqa: E402
    expected_return_files as nutimo_handoff_expected_return_files,
    payload as nutimo_handoff_payload,
    rows_as_dicts as nutimo_handoff_rows,
    target_carriers as nutimo_handoff_target_carriers,
    validate_rows as validate_nutimo_handoff_rows,
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


def test_power_law_orbital_harmonics_to_e2() -> None:
    sym = forcing_symbols()
    harmonics = power_law_orbital_harmonics_to_e2()
    p, e, n, t = sym["p"], sym["e"], sym["n"], sym["t"]

    assert sp.simplify(harmonics["constant"] - (1 + p * (p - 1) * e**2 / 4)) == 0
    assert sp.simplify(harmonics["first_harmonic"] - p * e * sp.cos(n * t)) == 0
    assert (
        sp.simplify(
            harmonics["second_harmonic"]
            - p * (p + 3) * e**2 * sp.cos(2 * n * t) / 4
        )
        == 0
    )


def test_power_law_orbital_harmonics_match_series_derivation() -> None:
    assert power_law_series_residual_to_e2() == 0


def test_observable_projection_recovers_calibrated_transfer_pair() -> None:
    sym = forcing_symbols()
    components = observable_harmonic_components()

    assert (
        sp.simplify(
            components["cos_coefficient"] / (sym["Lambda"] * sym["F_k"])
            - components["calibrated_cos"]
        )
        == 0
    )
    assert (
        sp.simplify(
            components["sin_coefficient"] / (sym["Lambda"] * sym["F_k"])
            - components["calibrated_sin"]
        )
        == 0
    )


def test_orbital_two_harmonic_dotf_obstruction() -> None:
    sym = forcing_symbols()
    obstruction = orbital_two_harmonic_dotf_obstruction()
    expected = (
        3
        * sym["beta"]
        * sym["tau_chi"] ** 3
        * sym["n"] ** 2
        / (
            (1 + sym["tau_chi"] ** 2 * sym["n"] ** 2)
            * (1 + 4 * sym["tau_chi"] ** 2 * sym["n"] ** 2)
        )
    )

    assert sp.simplify(obstruction - expected) == 0
    assert sp.simplify(obstruction.subs(sym["beta"], 0)) == 0
    assert sp.simplify(obstruction.subs(sym["tau_chi"], 0)) == 0
    assert sp.simplify(obstruction.subs(sym["n"], 0)) == 0


def test_forcing_dictionary_rows_are_classified() -> None:
    rows = forcing_dictionary_rows()
    validate_forcing_dictionary_rows(rows)
    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["forcing_class"] for row in rows)
    assert all(row["observable_projection"] for row in rows)


def test_forcing_dictionary_payload_has_required_columns() -> None:
    data = forcing_dictionary_payload()
    schema = set(data["schema"])

    assert "forcing_class" in schema
    assert "frequency_source" in schema
    assert "observable_projection" in schema
    assert "missing_assumption" in schema


def test_nonlinear_drive_generates_expected_rhs_sidebands() -> None:
    sym = nonlinear_sideband_symbols()
    sidebands = nonlinear_drive_rhs_sidebands()

    assert sidebands["2Omega1"] == sym["beta_F2"] * sym["F1"] ** 2 / 2
    assert sidebands["2Omega2"] == sym["beta_F2"] * sym["F2"] ** 2 / 2
    assert sidebands["Omega1_plus_Omega2"] == sym["beta_F2"] * sym["F1"] * sym["F2"]
    assert sidebands["Omega1_minus_Omega2"] == sym["beta_F2"] * sym["F1"] * sym["F2"]
    assert sp.simplify(sidebands["Omega1_plus_Omega2"].subs(sym["beta_F2"], 0)) == 0


def test_nonlinear_drive_sideband_response_has_relaxation_phase() -> None:
    sym = nonlinear_sideband_symbols()
    sidebands = nonlinear_drive_chi_sidebands()
    cos_coeff, sin_coeff = sidebands["Omega1_plus_Omega2"]
    nu = sym["Omega1"] + sym["Omega2"]

    assert sp.simplify(sin_coeff / cos_coeff - nu * sym["tau_chi"]) == 0
    assert sp.simplify(cos_coeff.subs(sym["F1"], 0)) == 0
    assert sp.simplify(sin_coeff.subs(sym["F2"], 0)) == 0


def test_nonlinear_readout_sum_sideband_collapses_when_lambdas_zero() -> None:
    sym = nonlinear_sideband_symbols()
    sideband = nonlinear_readout_sum_sideband()
    substitutions = {
        sym["lambda_Fchi"]: 0,
        sym["lambda_chi2"]: 0,
    }

    assert sideband["cos"] != 0
    assert sideband["sin"] != 0
    assert sp.simplify(sideband["cos"].subs(substitutions)) == 0
    assert sp.simplify(sideband["sin"].subs(substitutions)) == 0


def test_nonlinear_readout_difference_sideband_collapses_when_drive_missing() -> None:
    sym = nonlinear_sideband_symbols()
    sideband = nonlinear_readout_difference_sideband()

    assert sideband["cos"] != 0
    assert sideband["sin"] != 0
    assert sp.simplify(sideband["cos"].subs(sym["F1"], 0)) == 0
    assert sp.simplify(sideband["sin"].subs(sym["F2"], 0)) == 0


def test_orbital_n_2n_mixing_creates_3n_sideband() -> None:
    sym = nonlinear_sideband_symbols()
    amplitudes = orbital_n_2n_input_amplitudes()
    sideband = orbital_3n_nonlinear_drive_amplitude()
    expected = sym["beta_F2"] * sym["p"] ** 2 * (sym["p"] + 3) * sym["e"] ** 3 / 4

    assert amplitudes["n"] == sym["p"] * sym["e"]
    assert amplitudes["2n"] == sym["p"] * (sym["p"] + 3) * sym["e"] ** 2 / 4
    assert sp.simplify(sideband - expected) == 0
    assert sp.simplify(sideband.subs(sym["e"], 0)) == 0
    assert sp.simplify(sideband.subs(sym["p"], 0)) == 0
    assert sp.simplify(sideband.subs(sym["p"], -3)) == 0


def test_linear_projection_cannot_create_absent_sideband() -> None:
    sym = nonlinear_sideband_symbols()

    assert linear_projection_sideband_output() == 0
    assert linear_projection_sideband_output(sym["F1"]) == sym["F1"] * sym["Lambda_side"]


def test_nonlinear_sideband_rows_are_classified() -> None:
    rows = nonlinear_sideband_rows()
    validate_nonlinear_sideband_rows(rows)
    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["generated_frequencies"] for row in rows)
    assert all(row["collapse_condition"] for row in rows)


def test_nonlinear_sideband_payload_has_required_columns() -> None:
    data = nonlinear_sideband_payload()
    schema = set(data["schema"])

    assert "generated_frequencies" in schema
    assert "projection_note" in schema
    assert data["linear_projection_of_absent_sideband"] == "0"


def test_range_projection_deprojects_to_relaxation_transfer() -> None:
    assert range_deprojection_residual() == 0


def test_range_projection_satisfies_linear_channel_equation() -> None:
    assert range_equation_residual() == 0


def test_range_relaxation_pole_survives_unless_couplings_collapse() -> None:
    sym = projection_symbols()
    residue = range_relaxation_pole_residue()

    assert residue != 0
    assert sp.simplify(residue.subs(sym["Gamma"], 0)) == 0
    assert sp.simplify(residue.subs(sym["beta"], 0)) == 0


def test_observed_range_transfer_contains_projection_factor() -> None:
    observed = observed_range_transfer()
    expected = sp.factor(range_projection() * projection_relaxation_transfer())

    assert sp.simplify(observed - expected) == 0


def test_arbitrary_projection_solution_fits_pointwise_observable() -> None:
    sym = projection_symbols()
    lambda_solution = arbitrary_projection_solution()
    reconstructed = sp.simplify(
        lambda_solution * projection_relaxation_transfer() * sym["F_hat"]
    )

    assert sp.simplify(reconstructed - sym["O_hat"]) == 0


def test_projection_channel_rows_are_classified() -> None:
    rows = projection_channel_rows()
    validate_projection_channel_rows(rows)
    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["projection_channel"] for row in rows)
    assert all(row["verdict"] for row in rows)


def test_projection_channel_payload_has_required_columns() -> None:
    data = projection_channel_payload()
    schema = set(data["schema"])

    assert "projection_channel" in schema
    assert "lambda_form" in schema
    assert "collapse_condition" in schema
    assert data["range_deprojection_residual"] == "0"
    assert data["range_equation_residual"] == "0"


def test_projection_pole_condition_is_range_denominator() -> None:
    sym = projection_symbols()
    assert projection_pole_condition() == sym["kappa"] ** 2 + sym["z"] ** 2


def test_triple_bridge_rows_are_classified() -> None:
    rows = triple_bridge_rows()
    validate_triple_bridge_rows(rows)

    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["N_frequencies"] for row in rows)
    assert all(row["verdict"] for row in rows)


def test_triple_bridge_two_carriers_break_low_real_derivative_orders() -> None:
    assert real_two_frequency_verdict(1)["verdict"] == "distinguishable"
    assert real_two_frequency_verdict(2)["verdict"] == "distinguishable"
    assert real_two_frequency_verdict(3)["verdict"] == "degenerate at two carriers"


def test_triple_bridge_two_carriers_are_not_enough_for_complex_degree1() -> None:
    assert complex_two_frequency_verdict(0)["verdict"] == "distinguishable"
    assert complex_two_frequency_verdict(1)["verdict"] == "degenerate at two carriers"


def test_triple_bridge_transfer_uses_distinct_inner_outer_samples() -> None:
    transfer = transfer_at_carriers()

    assert transfer["G_in"] != transfer["G_out"]
    assert transfer["difference"] != 0


def test_triple_bridge_real_odd_residual_has_expected_collapse_limits() -> None:
    residual = two_carrier_real_odd_residual_order1()
    sym = real_odd_symbols(1)
    u0 = sym.u_nodes[0]

    assert residual != 0
    assert sp.simplify(residual.subs(sym.beta, 0)) == 0
    assert sp.simplify(residual.subs(sym.tau, 0)) == 0
    assert sp.simplify(residual.subs(sym.u_star, u0)) == 0


def test_triple_bridge_payload_records_carrier_inventory() -> None:
    data = triple_bridge_payload()

    assert "forcing_model" in data
    assert "transfer_at_carriers" in data
    assert "real_degree_boundaries" in data
    assert data["real_degree_boundaries"]["1"]["verdict"] == "distinguishable"


def test_triple_gr_inventory_rows_are_classified() -> None:
    rows = triple_gr_inventory_rows()
    validate_triple_gr_inventory_rows(rows)

    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["carrier_set"] for row in rows)
    assert all(row["verdict"] for row in rows)


def test_triple_gr_inventory_has_three_base_carriers() -> None:
    inventory = triple_gr_carrier_inventory()

    assert inventory["inner_monopole"] == "Omega_in"
    assert inventory["outer_monopole"] == "Omega_out"
    assert inventory["outer_dipole_combination"] == "abs(Omega_in - Omega_out)"


def test_three_carrier_real_degree_boundary() -> None:
    assert real_k_carrier_verdict(4, 3)["verdict"] == "distinguishable"
    assert real_k_carrier_verdict(5, 3)["verdict"] == "degenerate at 3 carriers"


def test_three_carrier_complex_degree_boundary() -> None:
    assert complex_k_carrier_verdict(1, 3)["verdict"] == "distinguishable"
    assert complex_k_carrier_verdict(2, 3)["verdict"] == "degenerate at 3 carriers"


def test_three_carrier_conditions_include_resonance_exclusions() -> None:
    conditions = three_carrier_conditions()

    assert "Omega_in != Omega_out" in conditions
    assert "Omega_in != 2 Omega_out" in conditions
    assert "Omega_out != 2 Omega_in" in conditions


def test_triple_gr_inventory_payload_records_boundaries() -> None:
    data = triple_gr_inventory_payload()

    assert "real_degree_boundaries_for_three_carriers" in data
    assert "complex_degree_boundaries_for_three_carriers" in data
    assert (
        data["real_degree_boundaries_for_three_carriers"]["4"]["verdict"]
        == "distinguishable"
    )
    assert (
        data["complex_degree_boundaries_for_three_carriers"]["1"]["verdict"]
        == "distinguishable"
    )


def test_triple_projection_gate_rows_are_classified() -> None:
    rows = triple_projection_gate_rows()
    validate_triple_projection_gate_rows(rows)

    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["bridge_verdict"] for row in rows)
    assert all(row["runtime_relevance"] for row in rows)


def test_triple_projection_gate_preserves_calibrated_projection() -> None:
    rows = {row["term"]: row for row in triple_projection_gate_rows()}

    assert rows["calibrated_common_real_scale"]["bridge_verdict"] == "distinguishable"
    assert rows["common_real_scale_unknown"]["runtime_relevance"] == "runtime-motivated"
    assert gate_calibrated_deprojection_residual() == 0
    assert gate_range_deprojection_residual() == 0


def test_triple_projection_gate_arbitrary_projection_collapses() -> None:
    rows = {row["term"]: row for row in triple_projection_gate_rows()}

    arbitrary = rows["arbitrary_per_carrier_complex_projection"]
    assert arbitrary["bridge_verdict"] == "collapse"
    assert arbitrary["runtime_relevance"] == "not-runtime-motivated"
    assert gate_arbitrary_projection_reconstruction_residual() == 0


def test_triple_projection_gate_keeps_finite_shared_nuisance_conditional() -> None:
    rows = {row["term"]: row for row in triple_projection_gate_rows()}

    assert rows["finite_real_shared_geometry"]["bridge_verdict"] == "conditional"
    assert rows["range_like_shared_projection_unknown"]["bridge_verdict"] == "conditional"
    assert rows["finite_complex_shared_polynomial"]["bridge_verdict"] == "conditional"


def test_triple_projection_gate_payload_records_projection_boundary() -> None:
    data = triple_projection_gate_payload()

    assert "projection_mapping" in data
    assert "projection_model_forms" in data
    assert "nuisance_dimension_ledger" in data
    assert data["arbitrary_projection_reconstruction_residual"] == "0"
    assert (
        triple_projection_mapping()["observed_carrier"]
        == "O_k=Lambda_k(theta) G(z_k) F_k"
    )
    assert (
        projection_model_forms()["arbitrary_per_carrier_complex"]
        == "Lambda_k independent for each carrier"
    )


def test_triple_projection_manifold_rows_are_classified() -> None:
    rows = triple_projection_manifold_rows()
    validate_triple_projection_manifold_rows(rows)

    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["generic_jacobian_rank"] for row in rows)
    assert all(row["runtime_worthiness"] for row in rows)


def test_triple_projection_manifold_rank_boundary() -> None:
    ranks = triple_projection_rank_summary()

    assert ranks["observable_vector_real_dimension"] == 6
    assert ranks["phase_locked_rank"] == 5
    assert ranks["linked_amplitude_rank"] == 5
    assert ranks["arbitrary_complex_rank"] == 6


def test_triple_projection_manifold_phase_lock_constraint() -> None:
    assert phase_lock_residual() == 0


def test_triple_projection_manifold_runtime_worthiness_gate() -> None:
    rows = {row["term"]: row for row in triple_projection_manifold_rows()}

    assert rows["phase_locked_outer_dipole_projection"]["bridge_verdict"] == "conditional"
    assert (
        rows["phase_locked_outer_dipole_projection"]["runtime_worthiness"]
        == "runtime-motivated"
    )
    assert (
        rows["effective_per_carrier_complex_projection"]["bridge_verdict"]
        == "collapse"
    )
    assert (
        rows["effective_per_carrier_complex_projection"]["runtime_worthiness"]
        == "not-runtime-motivated"
    )


def test_triple_projection_manifold_payload_records_gate_rule() -> None:
    data = triple_projection_manifold_payload()

    assert data["phase_lock_residual"] == "0"
    assert "runtime_worthiness_rule" in data
    assert "rank_summary" in data
    assert data["rank_summary"]["arbitrary_complex_rank"] == 6


def test_named_timing_model_rows_are_classified() -> None:
    rows = named_timing_model_rows()
    validate_named_timing_model_rows(rows)

    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["bridge_verdict"] for row in rows)
    assert all(row["runtime_worthiness"] for row in rows)


def test_named_timing_model_core_is_conditional_runtime_motivated() -> None:
    rows = {row["term"]: row for row in named_timing_model_rows()}
    core = rows["named_model_projection_class_verdict"]

    assert core["bridge_verdict"] == "conditional"
    assert core["runtime_worthiness"] == "runtime-motivated"
    assert "finite shared" in core["projection_nuisance_class"]


def test_named_timing_model_harmonic_specialcase_collapses() -> None:
    rows = {row["term"]: row for row in named_timing_model_rows()}
    specialcase = rows["rn_pl_specialcase_harmonic_nuisance"]

    assert specialcase["bridge_verdict"] == "collapse"
    assert specialcase["runtime_worthiness"] == "not-runtime-motivated"
    assert "per-harmonic" in specialcase["projection_nuisance_class"]


def test_named_timing_model_payload_records_sources_and_verdict() -> None:
    data = named_timing_model_payload()
    summary = named_timing_verdict_summary(data["rows"])
    releases = named_timing_source_release_ledger()

    assert data["source_releases"]["2025_zenodo"] == releases["2025_zenodo"]
    assert "source_releases" in data
    assert summary["standard_core_runtime_worthiness"] == "runtime-motivated"
    assert summary["specialcase_bridge_verdict"] == "collapse"


def test_nutimo_runtime_pilot_rows_are_classified() -> None:
    rows = nutimo_runtime_pilot_rows()
    validate_nutimo_runtime_pilot_rows(rows)

    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["bridge_verdict_if_pass"] for row in rows)
    assert all(row["bridge_verdict_if_fail"] for row in rows)


def test_nutimo_runtime_pilot_enforces_configuration_closure() -> None:
    rows = {row["gate"]: row for row in nutimo_runtime_pilot_rows()}
    config = rows["configuration_closure_standard_core"]

    assert config["bridge_verdict_if_pass"] == "conditional"
    assert config["bridge_verdict_if_fail"] == "collapse"
    assert "harmonic" in config["fail_condition"]


def test_nutimo_runtime_pilot_rank_gate_boundary() -> None:
    alive = jacobian_rank_gate_verdict(
        projection_rank=5,
        dynamic_chi_column_in_span=False,
    )
    rank_collapse = jacobian_rank_gate_verdict(
        projection_rank=6,
        dynamic_chi_column_in_span=False,
    )
    span_collapse = jacobian_rank_gate_verdict(
        projection_rank=5,
        dynamic_chi_column_in_span=True,
    )

    assert alive["bridge_verdict"] == "conditional"
    assert alive["runtime_worthiness"] == "runtime-motivated"
    assert rank_collapse["bridge_verdict"] == "collapse"
    assert span_collapse["bridge_verdict"] == "collapse"


def test_nutimo_runtime_pilot_payload_records_required_artifacts() -> None:
    data = nutimo_runtime_pilot_payload()
    artifacts = nutimo_pilot_minimal_external_artifacts()

    assert data["stage_order"] == nutimo_pilot_stage_order()
    assert "finite_jacobian" in data["minimal_external_artifacts"]
    assert "dynamic_chi_column" in data["minimal_external_artifacts"]
    assert data["minimal_external_artifacts"]["finite_jacobian"] == artifacts["finite_jacobian"]
    assert "rank_gate_examples" in data


def test_nutimo_handoff_rows_are_classified() -> None:
    rows = nutimo_handoff_rows()
    validate_nutimo_handoff_rows(rows)

    assert rows
    assert all(row["claim_status"] for row in rows)
    assert all(row["requested_artifact"] for row in rows)
    assert all(row["validation_rule"] for row in rows)


def test_nutimo_handoff_requires_configuration_and_jacobian() -> None:
    rows = {row["packet_item"]: row for row in nutimo_handoff_rows()}

    assert rows["configuration_manifest"]["requested_artifact"] == "configuration_manifest.json"
    assert rows["finite_parameter_jacobian"]["requested_artifact"] == "finite_jacobian.npy_or_tsv"
    assert "RN_PL" in rows["configuration_manifest"]["validation_rule"]


def test_nutimo_handoff_payload_records_return_files_and_stops() -> None:
    data = nutimo_handoff_payload()
    expected = nutimo_handoff_expected_return_files()

    assert data["target_carriers"] == nutimo_handoff_target_carriers()
    assert data["expected_return_files"] == expected
    assert "carrier projection rank is 6/6" in data["hard_stop_rules"]
    assert "configuration_manifest.json" in data["expected_return_files"]


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
    test_power_law_orbital_harmonics_to_e2()
    test_power_law_orbital_harmonics_match_series_derivation()
    test_observable_projection_recovers_calibrated_transfer_pair()
    test_orbital_two_harmonic_dotf_obstruction()
    test_forcing_dictionary_rows_are_classified()
    test_forcing_dictionary_payload_has_required_columns()
    test_nonlinear_drive_generates_expected_rhs_sidebands()
    test_nonlinear_drive_sideband_response_has_relaxation_phase()
    test_nonlinear_readout_sum_sideband_collapses_when_lambdas_zero()
    test_nonlinear_readout_difference_sideband_collapses_when_drive_missing()
    test_orbital_n_2n_mixing_creates_3n_sideband()
    test_linear_projection_cannot_create_absent_sideband()
    test_nonlinear_sideband_rows_are_classified()
    test_nonlinear_sideband_payload_has_required_columns()
    test_range_projection_deprojects_to_relaxation_transfer()
    test_range_projection_satisfies_linear_channel_equation()
    test_range_relaxation_pole_survives_unless_couplings_collapse()
    test_observed_range_transfer_contains_projection_factor()
    test_arbitrary_projection_solution_fits_pointwise_observable()
    test_projection_channel_rows_are_classified()
    test_projection_channel_payload_has_required_columns()
    test_projection_pole_condition_is_range_denominator()
    test_triple_bridge_rows_are_classified()
    test_triple_bridge_two_carriers_break_low_real_derivative_orders()
    test_triple_bridge_two_carriers_are_not_enough_for_complex_degree1()
    test_triple_bridge_transfer_uses_distinct_inner_outer_samples()
    test_triple_bridge_real_odd_residual_has_expected_collapse_limits()
    test_triple_bridge_payload_records_carrier_inventory()
    test_triple_gr_inventory_rows_are_classified()
    test_triple_gr_inventory_has_three_base_carriers()
    test_three_carrier_real_degree_boundary()
    test_three_carrier_complex_degree_boundary()
    test_three_carrier_conditions_include_resonance_exclusions()
    test_triple_gr_inventory_payload_records_boundaries()
    test_triple_projection_gate_rows_are_classified()
    test_triple_projection_gate_preserves_calibrated_projection()
    test_triple_projection_gate_arbitrary_projection_collapses()
    test_triple_projection_gate_keeps_finite_shared_nuisance_conditional()
    test_triple_projection_gate_payload_records_projection_boundary()
    test_triple_projection_manifold_rows_are_classified()
    test_triple_projection_manifold_rank_boundary()
    test_triple_projection_manifold_phase_lock_constraint()
    test_triple_projection_manifold_runtime_worthiness_gate()
    test_triple_projection_manifold_payload_records_gate_rule()
    test_named_timing_model_rows_are_classified()
    test_named_timing_model_core_is_conditional_runtime_motivated()
    test_named_timing_model_harmonic_specialcase_collapses()
    test_named_timing_model_payload_records_sources_and_verdict()
    test_nutimo_runtime_pilot_rows_are_classified()
    test_nutimo_runtime_pilot_enforces_configuration_closure()
    test_nutimo_runtime_pilot_rank_gate_boundary()
    test_nutimo_runtime_pilot_payload_records_required_artifacts()
    test_nutimo_handoff_rows_are_classified()
    test_nutimo_handoff_requires_configuration_and_jacobian()
    test_nutimo_handoff_payload_records_return_files_and_stops()
    print("dynamic chi symbolic checks passed")


if __name__ == "__main__":
    main()
