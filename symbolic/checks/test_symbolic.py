from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

SYMBOLIC_ROOT = Path(__file__).resolve().parents[1]
if str(SYMBOLIC_ROOT) not in sys.path:
    sys.path.insert(0, str(SYMBOLIC_ROOT))

from enumerate_contractions_delta4 import (
    enumerate_contraction_classes,
    gradient_sector_classes,
    mixed_time_derivative_classes,
)
from enumerate_basis import enumerate_minimal_scalar_monomials
from eb_sector_delta4 import eb_summary
from eb_survivor_rank_check import eb_rank_summary
from es_sector_delta4 import es_summary
from es_survivor_rank_check import es_rank_summary
from primitive_family_attack import primitive_attack_summary
from shift_scalar_sector_delta4 import shift_scalar_summary
from nonlinear_comparator_audit import (
    audit_rows,
    dynamic_sideband_transfer,
    linear_dynamic_transfer,
    orbital_three_n_ratio,
    shared_tau_phase_law,
    sideband_pole_interpolation_obstruction,
    static_single_line_mimic_coefficient,
    symbols as nonlinear_symbols,
    two_tone_sideband_ratio,
)
from shared_tau_ratio_audit import (
    bivariate_monomial_count,
    g_linear as shared_g_linear,
    interpolation_budget,
    no_sideband_limit,
    ratio_audit_rows,
    shared_tau_diagnostic_residual,
    shared_tau_ratio,
    static_ratio_template,
    symbols as shared_tau_symbols,
    tau_zero_limit,
)
from sample_budget_audit import (
    classification_rows,
    classify_design,
    compute_sample_budget,
    orbital_design,
    two_tone_design,
)
from orbital_harmonic_budget_audit import (
    clean_sum_sideband_pairs,
    exact_in_e_clean_new_line_count,
    find_minimum_orbital_design,
    classify_orbital_harmonic_design,
)
from second_order_mode_response import (
    cos_sin_components as second_order_cos_sin_components,
    first_order_limit as second_order_first_order_limit,
    low_frequency_series as second_order_low_frequency_series,
    phase_tangent as second_order_phase_tangent,
    resonance_frequency_squared,
    symbols as second_order_symbols,
    transfer_function as second_order_transfer_function,
)
from second_order_projection_audit import (
    acceleration_diagnostic_residual,
    deproject_acceleration,
    deproject_range,
    projection_audit_rows,
    projection_nuisance_budget,
    range_diagnostic_residual,
)
from resonant_comparator_audit import (
    dynamic_peak_condition,
    exact_polynomial_collapse_residual,
    line_shape_budget,
    phase_wrap_pole,
    resonance_audit_rows,
)
from exact_in_e_resonant_forcing_audit import (
    classify_resonant_forcing_design,
    default_design_rows as exact_in_e_design_rows,
    harmonic_bracket,
    harmonics,
    minimum_cutoff_for_bracket,
    minimum_cutoff_for_budget_and_bracket,
)
from amplitude_weighted_resonant_design import (
    classify_amplitude_design,
    exact_in_e_coefficient_formula,
    forcing_cosine_coefficient,
    minimum_amplitude_cutoff_design,
    weighted_harmonics,
)
from physical_detectability_map import (
    classify_physical_design,
    default_physical_design_rows,
    harmonic_noise_sigma,
    minimum_physical_amplitude_scale,
)
from nonlinear_second_order_mode import (
    audit_rows as nonlinear_second_order_rows,
    denominator as nonlinear_second_order_denominator,
    drive_sideband_interpolation_obstruction,
    internal_linear_transfer as nonlinear_second_order_internal_transfer,
    nonlinear_drive_sideband_transfer,
    nonlinear_readout_sideband_transfer,
    nonlinear_sideband_budget,
    observable_linear_transfer as nonlinear_second_order_observable_transfer,
    shared_denominator_residual,
    symbols as nonlinear_second_order_symbols,
)
from nonlinear_second_order_detectability import (
    classify_nonlinear_detectability,
    default_relative_rows as nonlinear_detectability_relative_rows,
    default_snr_rows as nonlinear_detectability_snr_rows,
    generated_forcing_coefficient,
    generated_line_rows,
    minimum_generated_amplitude_scale,
    nonlinear_drive_transfer_magnitude,
)
from component_separability_audit import (
    JointModelParameters,
    classify_separability_design,
    default_design_rows as component_separability_rows,
    jacobian_matrix,
    joint_observable,
    minimum_separable_cutoff,
    usable_harmonics as separability_usable_harmonics,
)
from nonlinear_robustness_map import (
    classify_robustness_row,
    representative_rows,
    robustness_grid,
    summarize_rows,
)
from fractions import Fraction
from normal_form_reduce import (
    operator_symbols,
    reduce_algebraic_identities,
    reduce_lower_order_eom,
    reduce_to_normal_form,
    reduce_total_derivatives,
)
from survivor_rank_check import rank_summary
from sensitivity_expand import make_quadratic_jet
from worldline_expand import build_worldline_model


def test_symmetric_quadratic_jet() -> None:
    y1, y2 = sp.symbols("Y1 Y2")
    jet = make_quadratic_jet("s", (y1, y2), mass_scale=sp.Symbol("m0"))
    assert jet.quadratic[0, 1] == jet.quadratic[1, 0]


def test_worldline_force_structure() -> None:
    model = build_worldline_model()
    x = model["x"]
    y1, y2 = model["invariants"]
    e2 = model["tidal_scalar"]
    m0 = model["m0"]
    lambda_e = model["lambda_E"]
    jet = model["jet"]

    s1 = jet.linear[0, 0]
    s2 = jet.linear[1, 0]
    s11 = jet.quadratic[0, 0]
    s12 = jet.quadratic[0, 1]
    s22 = jet.quadratic[1, 1]

    expected_force = sp.expand(
        m0 * (s1 + s11 * y1 + s12 * y2) * sp.diff(y1, x)
        + m0 * (s2 + s12 * y1 + s22 * y2) * sp.diff(y2, x)
        + sp.Rational(1, 2) * lambda_e * sp.diff(e2, x)
    )
    assert sp.simplify(model["force"] - expected_force) == 0


def test_basis_enumeration_order_four() -> None:
    monomials = enumerate_minimal_scalar_monomials(max_weight=4)
    labels = [monomial.label for monomial in monomials]
    assert labels == ["1", "E2", "E3", "E2^2", "dotE2", "gradE2"]
    assert all(monomial.weight <= 4 for monomial in monomials)


def test_total_derivative_reduction() -> None:
    ops = operator_symbols()
    assert reduce_total_derivatives(ops["E_DtE"]) == 0
    assert reduce_total_derivatives(ops["Dt2_E2"]) == 0
    assert reduce_total_derivatives(ops["E_Dt2E"]) == -ops["dotE2"]
    assert reduce_total_derivatives(ops["TrE2DtE"]) == 0


def test_lower_order_eom_reduction() -> None:
    ops = operator_symbols()
    expr = (
        ops["a2"]
        + ops["aEa"]
        + ops["aDivE"]
        + ops["aDtEa"]
        + ops["a2E2"]
        + ops["aE2a"]
        + ops["a4"]
        + ops["aEGradE_1"]
        + ops["aEGradE_2"]
        + ops["aEGradE_3"]
    )
    assert reduce_lower_order_eom(expr) == 0


def test_algebraic_reduction() -> None:
    ops = operator_symbols()
    assert reduce_algebraic_identities(ops["E4"]) == sp.Rational(1, 2) * ops["E2"] ** 2


def test_combined_normal_form_reduction() -> None:
    ops = operator_symbols()
    expr = (
        ops["E4"]
        + ops["E_Dt2E"]
        + ops["TrE2DtE"]
        + ops["a2"]
        + ops["gradE2"]
        + ops["divE2"]
    )
    reduced = reduce_to_normal_form(expr)
    expected = (
        sp.Rational(1, 2) * ops["E2"] ** 2
        - ops["dotE2"]
        + ops["gradE2"]
        + ops["divE2"]
    )
    assert sp.expand(reduced - expected) == 0


def test_contraction_enumeration_counts() -> None:
    classes = enumerate_contraction_classes()
    assert len(classes) == 21


def test_gradient_sector_audit() -> None:
    labels = [item.label for item in gradient_sector_classes()]
    assert labels == ["divE2", "gradE2", "mixedGradE2"]


def test_mixed_time_derivative_audit() -> None:
    items = mixed_time_derivative_classes()
    assert len(items) == 1
    assert items[0].label == "TrE2DtE"
    assert items[0].classification == "Proven reducible"


def test_a_e_grade_labels_are_unique() -> None:
    labels = sorted(
        item.label
        for item in enumerate_contraction_classes()
        if item.signature == ("E", "GradE", "a")
    )
    assert labels == ["aEGradE_1", "aEGradE_2", "aEGradE_3"]


def test_survivor_rank_independence() -> None:
    summary = rank_summary()
    assert summary.total_rank == 7
    assert summary.e_sector_rank == 3
    assert summary.dt_sector_rank == 1
    assert summary.gradient_sector_rank == 3


def test_magnetic_family_attack_finds_new_survivor() -> None:
    summary = primitive_attack_summary()
    assert summary.smallest_new_survivor == "B2"
    assert "B2" in summary.new_survivor_labels
    assert "EB2" in summary.new_survivor_labels


def test_eb_sector_survivor_list() -> None:
    summary = eb_summary()
    assert summary.total_classes == 42
    assert summary.smallest_new_survivor == "B2"
    assert summary.surviving_labels == (
        "B2",
        "E2",
        "EB2",
        "E3",
        "B2^2",
        "E2B2",
        "EB_sq",
        "TrE2B2",
        "EBDtB",
        "dotB2",
        "dotE2",
        "E2^2",
        "divB2",
        "gradB2",
        "mixedGradB2",
        "divE2",
        "gradE2",
        "mixedGradE2",
    )


def test_eb_survivor_rank_correction() -> None:
    summary = eb_rank_summary()
    assert summary.raw_rank == 18
    assert summary.raw_count == 19
    assert summary.corrected_rank == 18
    assert summary.corrected_count == 18
    relation = str(summary.null_relation)
    assert "EBEB" in relation
    assert "TrE2B2" in relation


def test_es_sector_survivor_list() -> None:
    summary = es_summary()
    assert summary.total_classes == 72
    assert summary.smallest_new_survivor == "S"
    assert summary.surviving_labels == (
        "S",
        "B2",
        "E2",
        "S2",
        "EB2",
        "SB2",
        "E3",
        "SE2",
        "S3",
        "B2^2",
        "DtS_B2",
        "E2B2",
        "EB_sq",
        "TrE2B2",
        "SEB2",
        "S2B2",
        "EBDtB",
        "dotB2",
        "dotE2",
        "dotS2",
        "DtS_E2",
        "E2^2",
        "SE3",
        "S2E2",
        "divB2",
        "gradB2",
        "mixedGradB2",
        "divE2",
        "gradE2",
        "mixedGradE2",
        "divEGradS",
        "gradS2",
        "S4",
    )


def test_es_survivor_rank_independence() -> None:
    summary = es_rank_summary()
    assert summary.rank == 33
    assert summary.count == 33
    assert summary.nullity == 0
    assert "divEGradS" in summary.labels
    assert "S" in summary.labels


def test_shift_scalar_sector_survivors() -> None:
    summary = shift_scalar_summary()
    assert summary.total_classes == 52
    assert summary.first_new_weight == 4
    assert summary.first_new_labels == (
        "DtS_B2",
        "dotS2",
        "DtS_E2",
        "divEGradS",
        "gradS2",
    )
    assert summary.canonical_new_survivor == "dotS2"
    assert summary.surviving_labels == (
        "B2",
        "E2",
        "EB2",
        "E3",
        "B2^2",
        "DtS_B2",
        "E2B2",
        "EB_sq",
        "TrE2B2",
        "EBDtB",
        "dotB2",
        "dotE2",
        "dotS2",
        "DtS_E2",
        "E2^2",
        "divB2",
        "gradB2",
        "mixedGradB2",
        "divE2",
        "gradE2",
        "mixedGradE2",
        "divEGradS",
        "gradS2",
    )


def test_dynamic_sideband_phase_law_returns_tau() -> None:
    syms = nonlinear_symbols()
    assert sp.simplify(shared_tau_phase_law() - syms["tau_chi"]) == 0


def test_single_sideband_static_mimic_is_exact() -> None:
    syms = nonlinear_symbols()
    nu = syms["nu"]
    assert sp.simplify(
        static_single_line_mimic_coefficient(nu) - dynamic_sideband_transfer(sp.I * nu)
    ) == 0


def test_sideband_pole_interpolation_residual_collapse_limits() -> None:
    syms = nonlinear_symbols()
    obstruction = sideband_pole_interpolation_obstruction(order=1)
    assert obstruction.order == 1
    assert len(obstruction.nodes) == 2
    assert sp.simplify(obstruction.residual.subs(syms["C_dyn"], 0)) == 0
    assert sp.simplify(obstruction.residual.subs(syms["tau_chi"], 0)) == 0
    assert sp.simplify(obstruction.residual.subs(obstruction.target, obstruction.nodes[0])) == 0


def test_two_tone_ratio_has_shared_sideband_pole() -> None:
    syms = nonlinear_symbols()
    omega1 = syms["Omega1"]
    omega2 = syms["Omega2"]
    tau_chi = syms["tau_chi"]
    expected = syms["C_dyn"] / (
        (1 + sp.I * (omega1 + omega2) * tau_chi)
        * linear_dynamic_transfer(omega1)
        * linear_dynamic_transfer(omega2)
    )
    assert sp.simplify(two_tone_sideband_ratio() - expected) == 0


def test_orbital_three_n_ratio_has_shared_sideband_pole() -> None:
    syms = nonlinear_symbols()
    n = syms["n"]
    tau_chi = syms["tau_chi"]
    expected = syms["C_dyn"] / (
        (1 + 3 * sp.I * n * tau_chi)
        * linear_dynamic_transfer(n)
        * linear_dynamic_transfer(2 * n)
    )
    assert sp.simplify(orbital_three_n_ratio() - expected) == 0


def test_nonlinear_comparator_audit_rows_classify_boundaries() -> None:
    rows = audit_rows()
    assert len(rows) == 4
    assert rows[0].status == "Proven"
    assert "not unique" in rows[0].verdict
    assert any("shared sideband pole" in row.observable for row in rows)


def test_shared_tau_ratio_diagnostic_residual_vanishes() -> None:
    assert shared_tau_diagnostic_residual() == 0


def test_shared_tau_ratio_collapse_limits() -> None:
    syms = shared_tau_symbols()
    expected_tau_zero = syms["C_side"] / (syms["beta"] + syms["c_Y"]) ** 2
    assert sp.simplify(tau_zero_limit() - expected_tau_zero) == 0
    assert no_sideband_limit() == 0


def test_shared_tau_ratio_matches_definition() -> None:
    syms = shared_tau_symbols()
    u = syms["u"]
    v = syms["v"]
    tau_chi = syms["tau_chi"]
    expected = syms["C_side"] / (
        (1 + tau_chi * (u + v)) * shared_g_linear(u) * shared_g_linear(v)
    )
    assert sp.cancel(shared_tau_ratio() - expected) == 0


def test_static_ratio_template_and_budget_counts() -> None:
    syms = shared_tau_symbols()
    ratio = static_ratio_template(linear_order=1, sideband_degree=2)
    assert syms["u"] in ratio.free_symbols
    assert syms["v"] in ratio.free_symbols
    assert bivariate_monomial_count(2) == 6
    assert interpolation_budget(1, 2) == {
        "linear_complex_samples": 2,
        "linear_obstruction_sample": 3,
        "sideband_pair_samples": 6,
        "sideband_obstruction_pair": 7,
    }


def test_shared_tau_ratio_audit_rows_record_no_go_boundary() -> None:
    rows = ratio_audit_rows()
    assert len(rows) == 3
    assert rows[0].status == "Proven"
    assert "interpolation budget" in rows[1].verdict
    assert "per-sideband nuisance" in rows[2].collapse_boundary


def test_sample_budget_minimum_for_realistic_first_derivative_case() -> None:
    budget = compute_sample_budget(linear_order=1, sideband_degree=1, projection_nuisance=0)
    assert budget.linear_parameter_budget == 2
    assert budget.sideband_parameter_budget == 3
    assert budget.total_shared_budget == 5
    assert budget.minimum_linear_samples == 3
    assert budget.minimum_sideband_pairs == 4
    assert budget.minimum_total_samples == 7


def test_sample_budget_projection_nuisance_increases_sideband_requirement() -> None:
    budget = compute_sample_budget(linear_order=1, sideband_degree=1, projection_nuisance=2)
    assert budget.minimum_linear_samples == 3
    assert budget.minimum_sideband_pairs == 5
    assert budget.minimum_total_samples == 8


def test_orbital_design_is_underbudget_for_trivial_sideband_comparator() -> None:
    linear_samples, sideband_pairs = orbital_design()
    classification = classify_design(
        "orbital",
        linear_samples,
        sideband_pairs,
        linear_order=0,
        sideband_degree=0,
        projection_nuisance=0,
    )
    assert classification.verdict == "underbudget-no-go"
    assert classification.beats_linear_budget
    assert not classification.beats_sideband_budget


def test_two_tone_design_beats_only_trivial_comparator() -> None:
    linear_samples, sideband_pairs = two_tone_design(include_difference=False)
    trivial = classify_design(
        "two-tone",
        linear_samples,
        sideband_pairs,
        linear_order=0,
        sideband_degree=0,
        projection_nuisance=0,
    )
    realistic = classify_design(
        "two-tone",
        linear_samples,
        sideband_pairs,
        linear_order=1,
        sideband_degree=1,
        projection_nuisance=0,
    )
    assert trivial.verdict == "budget-breaking"
    assert realistic.verdict == "underbudget-no-go"
    assert not realistic.beats_linear_budget


def test_sample_budget_classification_rows_include_current_cases() -> None:
    rows = classification_rows()
    assert len(rows) == 15
    assert any(
        row.name == "two-tone sum-only"
        and row.budget.linear_order == 0
        and row.budget.sideband_degree == 0
        and row.verdict == "budget-breaking"
        for row in rows
    )
    assert all(
        row.verdict == "underbudget-no-go"
        for row in rows
        if row.name == "orbital n,2n -> 3n"
    )


def test_clean_orbital_pair_count_reproduces_current_case() -> None:
    assert clean_sum_sideband_pairs(2, 3) == ((1, 2, 3),)
    current = classify_orbital_harmonic_design(
        2,
        3,
        linear_order=1,
        sideband_degree=1,
        projection_nuisance=0,
    )
    assert current.classification.linear_samples == 2
    assert current.classification.sideband_pairs == 1
    assert current.classification.verdict == "underbudget-no-go"


def test_richer_orbital_minimum_for_first_derivative_linear_sideband() -> None:
    design = find_minimum_orbital_design(
        linear_order=1,
        sideband_degree=1,
        projection_nuisance=0,
    )
    assert design is not None
    assert design.linear_harmonics == 3
    assert design.sideband_harmonic_cutoff == 6
    assert len(design.clean_sideband_pairs) == 4
    assert design.classification.verdict == "budget-breaking"


def test_richer_orbital_projection_nuisance_minimum() -> None:
    design = find_minimum_orbital_design(
        linear_order=1,
        sideband_degree=1,
        projection_nuisance=2,
    )
    assert design is not None
    assert design.linear_harmonics == 4
    assert design.sideband_harmonic_cutoff == 6
    assert len(design.clean_sideband_pairs) == 4


def test_exact_in_e_has_no_clean_new_line_count() -> None:
    assert exact_in_e_clean_new_line_count() == 0


def test_second_order_transfer_and_first_order_limit() -> None:
    syms = second_order_symbols()
    z = syms["z"]
    expected = syms["alpha"] * syms["omega_chi_sq"] / (
        syms["mu_chi"] * z**2 + syms["gamma_chi"] * z + syms["omega_chi_sq"]
    )
    assert sp.simplify(second_order_transfer_function() - expected) == 0
    first_order_expected = syms["alpha"] / (
        1 + syms["gamma_chi"] * z / syms["omega_chi_sq"]
    )
    assert sp.simplify(second_order_first_order_limit() - first_order_expected) == 0


def test_second_order_phase_tangent_and_resonance() -> None:
    syms = second_order_symbols()
    omega = syms["Omega"]
    expected_phase = syms["gamma_chi"] * omega / (
        syms["omega_chi_sq"] - syms["mu_chi"] * omega**2
    )
    assert sp.simplify(second_order_phase_tangent() - expected_phase) == 0
    expected_peak = (
        syms["omega_chi_sq"] / syms["mu_chi"]
        - syms["gamma_chi"] ** 2 / (2 * syms["mu_chi"] ** 2)
    )
    assert sp.simplify(resonance_frequency_squared() - expected_peak) == 0


def test_second_order_low_frequency_series_coefficients() -> None:
    syms = second_order_symbols()
    z = syms["z"]
    series = second_order_low_frequency_series(order=4)
    assert sp.simplify(series.subs(z, 0) - syms["alpha"]) == 0
    assert sp.simplify(
        series.coeff(z, 1) + syms["alpha"] * syms["gamma_chi"] / syms["omega_chi_sq"]
    ) == 0
    assert sp.simplify(
        series.coeff(z, 2)
        - syms["alpha"]
        * (
            syms["gamma_chi"] ** 2 / syms["omega_chi_sq"] ** 2
            - syms["mu_chi"] / syms["omega_chi_sq"]
        )
    ) == 0


def test_second_order_cos_sin_components_have_common_denominator() -> None:
    syms = second_order_symbols()
    omega = syms["Omega"]
    cos_component, sin_component = second_order_cos_sin_components()
    expected_denominator = (
        (syms["omega_chi_sq"] - syms["mu_chi"] * omega**2) ** 2
        + syms["gamma_chi"] ** 2 * omega**2
    )
    assert sp.simplify(sp.denom(cos_component) - expected_denominator) == 0
    assert sp.simplify(sp.denom(sin_component) - expected_denominator) == 0


def test_second_order_projection_diagnostics_deproject_to_internal_transfer() -> None:
    assert acceleration_diagnostic_residual() == 0
    assert range_diagnostic_residual() == 0
    assert sp.simplify(deproject_acceleration() - second_order_transfer_function()) == 0
    assert sp.simplify(deproject_range() - second_order_transfer_function()) == 0


def test_second_order_projection_nuisance_budgets_and_rows() -> None:
    assert projection_nuisance_budget("acceleration") == 1
    assert projection_nuisance_budget("range") == 2
    rows = projection_audit_rows()
    assert len(rows) == 2
    assert rows[0].status == "Proven"
    assert "resonance survives" in rows[1].verdict


def test_resonant_line_shape_budget_counts_projection_nuisance() -> None:
    calibrated = line_shape_budget(polynomial_order=1, projection_nuisance=0)
    range_channel = line_shape_budget(polynomial_order=1, projection_nuisance=2)
    assert calibrated.static_parameter_budget == 2
    assert calibrated.minimum_frequency_samples == 3
    assert range_channel.static_parameter_budget == 4
    assert range_channel.minimum_frequency_samples == 5
    assert range_channel.resonance_bracket_samples == 5


def test_resonance_condition_and_phase_wrap_pole() -> None:
    syms = second_order_symbols()
    assert sp.simplify(
        dynamic_peak_condition()
        - (2 * syms["mu_chi"] * syms["omega_chi_sq"] - syms["gamma_chi"] ** 2)
    ) == 0
    assert sp.simplify(phase_wrap_pole() - syms["omega_chi_sq"] / syms["mu_chi"]) == 0


def test_exact_polynomial_collapse_residual_is_not_identity_generically() -> None:
    syms = second_order_symbols()
    residual = exact_polynomial_collapse_residual()
    assert residual.coeff(syms["z"], 4) == syms["mu_chi"] * sp.Symbol("p2")
    assert residual.subs({
        syms["mu_chi"]: 0,
        syms["gamma_chi"]: 0,
        sp.Symbol("p0"): syms["alpha"],
        sp.Symbol("p1"): 0,
        sp.Symbol("p2"): 0,
    }) == 0


def test_resonant_comparator_rows_record_budget_and_bracket_conditions() -> None:
    rows = resonance_audit_rows()
    assert len(rows) == 2
    assert "N+1+K" in rows[0].verdict
    assert "bracketing" in rows[1].verdict


def test_exact_in_e_harmonic_bracket_examples() -> None:
    assert harmonics(4) == (1, 2, 3, 4)
    half = harmonic_bracket(Fraction(3, 2), 4)
    assert half.lower_count == 1
    assert half.upper_count == 3
    assert not half.exact_hit
    assert half.bracketed
    integer_hit = harmonic_bracket(Fraction(4, 1), 4)
    assert integer_hit.exact_hit
    assert not integer_hit.bracketed
    integer_bracket = harmonic_bracket(Fraction(4, 1), 5)
    assert integer_bracket.bracketed


def test_exact_in_e_minimum_cutoffs_for_budget_and_bracket() -> None:
    assert minimum_cutoff_for_bracket(Fraction(3, 2)) == 2
    assert minimum_cutoff_for_bracket(Fraction(4, 1)) == 5
    assert minimum_cutoff_for_bracket(Fraction(3, 4)) is None
    assert minimum_cutoff_for_budget_and_bracket(Fraction(3, 2), 1, 1) == 4
    assert minimum_cutoff_for_budget_and_bracket(Fraction(3, 2), 1, 2) == 5
    assert minimum_cutoff_for_budget_and_bracket(Fraction(4, 1), 1, 1) == 5


def test_exact_in_e_design_classification_rows() -> None:
    acceleration = classify_resonant_forcing_design(
        "near 3/2 acceleration",
        Fraction(3, 2),
        polynomial_order=1,
        projection_nuisance=1,
    )
    sub_fundamental = classify_resonant_forcing_design(
        "subfundamental",
        Fraction(3, 4),
        polynomial_order=1,
        projection_nuisance=1,
    )
    assert acceleration.minimum_harmonic_cutoff == 4
    assert acceleration.verdict == "budget-breaking-if-usable"
    assert sub_fundamental.minimum_harmonic_cutoff is None
    assert sub_fundamental.verdict == "no single-system harmonic bracket"
    rows = exact_in_e_design_rows()
    assert len(rows) == 12
    assert any(row.minimum_harmonic_cutoff == 6 for row in rows)


def test_amplitude_weighted_exact_in_e_coefficient_formula_and_circular_limit() -> None:
    assert "A_k(p,e)" in exact_in_e_coefficient_formula()
    circular = forcing_cosine_coefficient(
        p=2.0,
        eccentricity=0.0,
        harmonic=1,
        samples=1024,
    )
    eccentric = forcing_cosine_coefficient(
        p=2.0,
        eccentricity=0.3,
        harmonic=1,
        samples=1024,
    )
    assert abs(circular) < 1.0e-12
    assert abs(eccentric) > 0.1


def test_amplitude_weighted_design_splits_low_and_moderate_eccentricity() -> None:
    low = classify_amplitude_design(
        "low-e",
        p=2.0,
        eccentricity=0.1,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=6,
        relative_cutoff=1.0e-3,
        polynomial_order=1,
        projection_nuisance=1,
        samples=2048,
    )
    moderate = classify_amplitude_design(
        "moderate-e",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=6,
        relative_cutoff=1.0e-3,
        polynomial_order=1,
        projection_nuisance=1,
        samples=2048,
    )
    assert low.verdict == "amplitude-underbudget-no-go"
    assert moderate.verdict == "amplitude-budget-breaking"
    assert moderate.lower_usable_count > 0
    assert moderate.upper_usable_count > 0


def test_amplitude_weighted_minimum_cutoff_example() -> None:
    design = minimum_amplitude_cutoff_design(
        "minimum",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        relative_cutoff=1.0e-3,
        polynomial_order=1,
        projection_nuisance=1,
        max_harmonic_cutoff=8,
        samples=2048,
    )
    assert design is not None
    assert design.harmonic_cutoff == 4
    assert design.usable_count >= design.budget.minimum_frequency_samples


def test_amplitude_weighted_harmonics_mark_unusable_high_threshold() -> None:
    rows = weighted_harmonics(
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=4,
        relative_cutoff=2.0,
        samples=1024,
    )
    assert rows
    assert not any(row.usable for row in rows)


def test_physical_detectability_noise_scales_by_harmonic() -> None:
    assert harmonic_noise_sigma(1.0e-3, 4, harmonic_noise_slope=0.5) == 2.0e-3


def test_physical_detectability_splits_low_and_high_signal_scale() -> None:
    low = classify_physical_design(
        "low-scale",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=6,
        amplitude_scale=0.1,
        base_noise_sigma=1.0e-3,
        snr_threshold=5.0,
        polynomial_order=1,
        projection_nuisance=1,
        samples=2048,
    )
    high = classify_physical_design(
        "high-scale",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=6,
        amplitude_scale=1.0,
        base_noise_sigma=1.0e-3,
        snr_threshold=5.0,
        polynomial_order=1,
        projection_nuisance=1,
        samples=2048,
    )
    assert low.verdict == "snr-underbudget-no-go"
    assert high.verdict == "physically-budget-breaking"
    assert high.detectable_count >= high.budget.minimum_frequency_samples


def test_physical_detectability_minimum_scale_examples() -> None:
    acceleration = minimum_physical_amplitude_scale(
        "minimum",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=6,
        base_noise_sigma=1.0e-3,
        snr_threshold=5.0,
        polynomial_order=1,
        projection_nuisance=1,
        samples=2048,
    )
    range_channel = minimum_physical_amplitude_scale(
        "minimum",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="range",
        harmonic_cutoff=6,
        base_noise_sigma=1.0e-6,
        snr_threshold=5.0,
        polynomial_order=1,
        projection_nuisance=2,
        samples=2048,
    )
    assert acceleration is not None
    assert range_channel is not None
    assert acceleration.verdict == "physically-budget-breaking"
    assert range_channel.verdict == "physically-budget-breaking"
    assert 0.8 < acceleration.amplitude_scale < 1.0
    assert 0.2 < range_channel.amplitude_scale < 0.4


def test_physical_detectability_default_rows_include_no_go_and_positive() -> None:
    rows = default_physical_design_rows()
    assert len(rows) == 6
    assert any(row.verdict == "snr-underbudget-no-go" for row in rows)
    assert any(row.verdict == "physically-budget-breaking" for row in rows)


def test_nonlinear_second_order_transfers_share_denominator() -> None:
    syms = nonlinear_second_order_symbols()
    z = syms["z"]
    denominator = nonlinear_second_order_denominator(z)
    expected_internal = syms["alpha"] * syms["omega_chi_sq"] / denominator
    expected_observable = syms["c_Y"] + syms["c_chi"] * expected_internal
    expected_sideband = syms["c_chi"] * syms["beta_F2"] / denominator
    assert sp.simplify(nonlinear_second_order_internal_transfer(z) - expected_internal) == 0
    assert sp.simplify(nonlinear_second_order_observable_transfer(z) - expected_observable) == 0
    assert sp.simplify(nonlinear_drive_sideband_transfer(z) - expected_sideband) == 0


def test_nonlinear_second_order_readout_sideband_law() -> None:
    syms = nonlinear_second_order_symbols()
    u = syms["u"]
    v = syms["v"]
    h_u = nonlinear_second_order_internal_transfer(u)
    h_v = nonlinear_second_order_internal_transfer(v)
    expected = syms["lambda_Fchi"] * (h_u + h_v) + syms["lambda_chi2"] * h_u * h_v
    assert sp.factor(nonlinear_readout_sideband_transfer(u, v) - expected) == 0


def test_nonlinear_second_order_shared_denominator_residual_vanishes() -> None:
    assert shared_denominator_residual() == 0


def test_nonlinear_second_order_interpolation_obstruction_boundaries() -> None:
    syms = nonlinear_second_order_symbols()
    obstruction = drive_sideband_interpolation_obstruction(order=1)
    assert obstruction.order == 1
    assert len(obstruction.nodes) == 2
    assert sp.simplify(obstruction.residual.subs(syms["beta_F2"], 0)) == 0
    assert sp.simplify(obstruction.residual.subs(syms["c_chi"], 0)) == 0
    assert sp.simplify(obstruction.residual.subs(obstruction.target, obstruction.nodes[0])) == 0
    assert sp.simplify(
        obstruction.residual.subs({
            syms["mu_chi"]: 0,
            syms["gamma_chi"]: 0,
        })
    ) == 0


def test_nonlinear_second_order_budget_and_rows() -> None:
    budget = nonlinear_sideband_budget(polynomial_order=1, projection_nuisance=1)
    assert budget.minimum_frequency_samples == 4
    rows = nonlinear_second_order_rows()
    assert len(rows) == 4
    assert rows[0].status == "Proven"
    assert "not unique" in rows[0].verdict
    assert any("shared" in row.target for row in rows)


def test_nonlinear_detectability_uses_exact_f_squared_coefficients() -> None:
    generated = generated_forcing_coefficient(
        p=2.0,
        eccentricity=0.3,
        harmonic=3,
        samples=1024,
    )
    direct = forcing_cosine_coefficient(
        p=4.0,
        eccentricity=0.3,
        harmonic=3,
        samples=1024,
    )
    assert abs(generated - direct) < 1.0e-15


def test_nonlinear_detectability_transfer_has_resonant_shape() -> None:
    transfer = nonlinear_drive_transfer_magnitude(
        harmonic=2,
        rho=Fraction(3, 2),
        damping=0.2,
    )
    expected = 1 / (((2.25 - 4.0) ** 2 + 0.4**2) ** 0.5)
    assert abs(transfer - expected) < 1.0e-15


def test_nonlinear_detectability_relative_audit_splits_eccentricity() -> None:
    low = classify_nonlinear_detectability(
        "low-e",
        p=2.0,
        eccentricity=0.1,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        generated_cutoff=6,
        criterion="relative",
        relative_cutoff=1.0e-3,
        sideband_order=1,
        projection_nuisance=1,
        samples=2048,
    )
    moderate = classify_nonlinear_detectability(
        "moderate-e",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        generated_cutoff=6,
        criterion="relative",
        relative_cutoff=1.0e-3,
        sideband_order=1,
        projection_nuisance=1,
        samples=2048,
    )
    assert low.verdict == "generated-underbudget-no-go"
    assert moderate.verdict == "nonlinear-generated-budget-breaking"
    assert moderate.usable_count >= moderate.budget.minimum_frequency_samples


def test_nonlinear_detectability_snr_minimum_scales() -> None:
    acceleration = minimum_generated_amplitude_scale(
        "minimum",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        generated_cutoff=6,
        base_noise_sigma=1.0e-3,
        snr_threshold=5.0,
        sideband_order=1,
        projection_nuisance=1,
        samples=2048,
    )
    range_channel = minimum_generated_amplitude_scale(
        "minimum",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="range",
        generated_cutoff=6,
        base_noise_sigma=1.0e-6,
        snr_threshold=5.0,
        sideband_order=1,
        projection_nuisance=2,
        samples=2048,
    )
    assert acceleration is not None
    assert range_channel is not None
    assert acceleration.verdict == "nonlinear-generated-budget-breaking"
    assert range_channel.verdict == "nonlinear-generated-budget-breaking"
    assert 0.3 < acceleration.amplitude_scale < 0.6
    assert 0.09 < range_channel.amplitude_scale < 0.2


def test_nonlinear_detectability_default_rows_include_no_go_and_positive() -> None:
    relative_rows = nonlinear_detectability_relative_rows()
    snr_rows = nonlinear_detectability_snr_rows()
    assert len(relative_rows) == 6
    assert len(snr_rows) == 6
    assert any(row.verdict == "generated-underbudget-no-go" for row in relative_rows)
    assert any(row.verdict == "nonlinear-generated-budget-breaking" for row in relative_rows)
    assert any(row.verdict == "generated-underbudget-no-go" for row in snr_rows)
    assert any(row.verdict == "nonlinear-generated-budget-breaking" for row in snr_rows)


def test_nonlinear_detectability_generated_rows_mark_unusable_high_cutoff() -> None:
    rows = generated_line_rows(
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        generated_cutoff=4,
        criterion="relative",
        relative_cutoff=2.0,
        samples=1024,
    )
    assert rows
    assert not any(row.usable for row in rows)


def test_component_separability_uses_same_usable_harmonics() -> None:
    harmonics = separability_usable_harmonics(
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=6,
        relative_cutoff=1.0e-3,
        samples=2048,
    )
    assert harmonics == (1, 2, 3, 4, 5)


def test_component_separability_jacobian_generated_column_adds_rank() -> None:
    design = classify_separability_design(
        "moderate-e",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=6,
        relative_cutoff=1.0e-3,
        projection_nuisance=1,
        samples=2048,
    )
    assert design.verdict == "component-separable-and-budget-ready"
    assert design.generated_adds_rank
    assert design.full_rank == design.parameter_count


def test_component_separability_distinguishes_underbudget_and_degenerate_cases() -> None:
    low_acceleration = classify_separability_design(
        "low-e",
        p=2.0,
        eccentricity=0.1,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        harmonic_cutoff=6,
        relative_cutoff=1.0e-3,
        projection_nuisance=1,
        samples=2048,
    )
    low_range_free_kappa = classify_separability_design(
        "low-e range",
        p=2.0,
        eccentricity=0.1,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="range",
        harmonic_cutoff=6,
        relative_cutoff=1.0e-3,
        projection_nuisance=2,
        include_range_projection_nuisance=True,
        samples=2048,
    )
    assert low_acceleration.verdict == "component-separable-but-budget-underdesigned"
    assert low_range_free_kappa.verdict == "generated-component-degenerate"


def test_component_separability_minimum_cutoffs() -> None:
    acceleration = minimum_separable_cutoff(
        "minimum",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        projection_nuisance=1,
        max_harmonic_cutoff=8,
        samples=2048,
    )
    range_channel = minimum_separable_cutoff(
        "minimum",
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="range",
        projection_nuisance=2,
        include_range_projection_nuisance=True,
        max_harmonic_cutoff=8,
        samples=2048,
    )
    assert acceleration is not None
    assert range_channel is not None
    assert acceleration.harmonic_cutoff == 4
    assert range_channel.harmonic_cutoff == 5


def test_component_separability_default_rows_cover_positive_and_negative() -> None:
    rows = component_separability_rows()
    assert len(rows) == 9
    assert any(row.verdict == "component-separable-and-budget-ready" for row in rows)
    assert any(row.verdict == "component-separable-but-budget-underdesigned" for row in rows)
    assert any(row.verdict == "generated-component-degenerate" for row in rows)


def test_component_separability_joint_observable_is_complex() -> None:
    value = joint_observable(
        harmonic=2,
        p=2.0,
        eccentricity=0.3,
        channel="acceleration",
        params=JointModelParameters(),
        samples=1024,
    )
    assert abs(value.imag) > 0


def test_component_separability_jacobian_shape() -> None:
    matrix = jacobian_matrix(
        harmonics=(1, 2, 3),
        p=2.0,
        eccentricity=0.3,
        channel="acceleration",
        params=JointModelParameters(),
        samples=1024,
    )
    assert matrix.shape == (6, 5)


def test_nonlinear_robustness_row_requires_budget_and_bracket() -> None:
    low = classify_robustness_row(
        p=2.0,
        eccentricity=0.1,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        projection_mode="acceleration",
        projection_nuisance=1,
        include_projection_nuisance=False,
        samples=1024,
    )
    moderate = classify_robustness_row(
        p=2.0,
        eccentricity=0.3,
        rho=Fraction(3, 2),
        damping=0.2,
        channel="acceleration",
        projection_mode="acceleration",
        projection_nuisance=1,
        include_projection_nuisance=False,
        samples=1024,
    )
    assert low.verdict == "component-separable-but-budget-underdesigned"
    assert low.bracketed
    assert moderate.verdict == "robust-positive"


def test_nonlinear_robustness_default_grid_counts_regions() -> None:
    rows = robustness_grid(samples=512)
    summary = summarize_rows(rows)
    assert summary.total == 162
    assert summary.verdict_counts["robust-positive"] == 107
    assert summary.verdict_counts["component-separable-but-budget-underdesigned"] == 40
    assert summary.verdict_counts["generated-component-degenerate"] == 15
    assert 0.6 < summary.positive_fraction < 0.7


def test_nonlinear_robustness_channel_counts_and_representatives() -> None:
    rows = robustness_grid(samples=512)
    summary = summarize_rows(rows)
    assert summary.channel_counts["acceleration"]["robust-positive"] == 37
    assert summary.channel_counts["range-free-kappa"]["generated-component-degenerate"] == 15
    representatives = representative_rows(rows)
    verdicts = {row.verdict for row in representatives}
    assert {"robust-positive", "generated-component-degenerate"}.issubset(verdicts)


def main() -> None:
    test_symmetric_quadratic_jet()
    test_worldline_force_structure()
    test_basis_enumeration_order_four()
    test_total_derivative_reduction()
    test_lower_order_eom_reduction()
    test_algebraic_reduction()
    test_combined_normal_form_reduction()
    test_contraction_enumeration_counts()
    test_gradient_sector_audit()
    test_mixed_time_derivative_audit()
    test_a_e_grade_labels_are_unique()
    test_survivor_rank_independence()
    test_magnetic_family_attack_finds_new_survivor()
    test_eb_sector_survivor_list()
    test_eb_survivor_rank_correction()
    test_es_sector_survivor_list()
    test_es_survivor_rank_independence()
    test_shift_scalar_sector_survivors()
    test_dynamic_sideband_phase_law_returns_tau()
    test_single_sideband_static_mimic_is_exact()
    test_sideband_pole_interpolation_residual_collapse_limits()
    test_two_tone_ratio_has_shared_sideband_pole()
    test_orbital_three_n_ratio_has_shared_sideband_pole()
    test_nonlinear_comparator_audit_rows_classify_boundaries()
    test_shared_tau_ratio_diagnostic_residual_vanishes()
    test_shared_tau_ratio_collapse_limits()
    test_shared_tau_ratio_matches_definition()
    test_static_ratio_template_and_budget_counts()
    test_shared_tau_ratio_audit_rows_record_no_go_boundary()
    test_sample_budget_minimum_for_realistic_first_derivative_case()
    test_sample_budget_projection_nuisance_increases_sideband_requirement()
    test_orbital_design_is_underbudget_for_trivial_sideband_comparator()
    test_two_tone_design_beats_only_trivial_comparator()
    test_sample_budget_classification_rows_include_current_cases()
    test_clean_orbital_pair_count_reproduces_current_case()
    test_richer_orbital_minimum_for_first_derivative_linear_sideband()
    test_richer_orbital_projection_nuisance_minimum()
    test_exact_in_e_has_no_clean_new_line_count()
    test_second_order_transfer_and_first_order_limit()
    test_second_order_phase_tangent_and_resonance()
    test_second_order_low_frequency_series_coefficients()
    test_second_order_cos_sin_components_have_common_denominator()
    test_second_order_projection_diagnostics_deproject_to_internal_transfer()
    test_second_order_projection_nuisance_budgets_and_rows()
    test_resonant_line_shape_budget_counts_projection_nuisance()
    test_resonance_condition_and_phase_wrap_pole()
    test_exact_polynomial_collapse_residual_is_not_identity_generically()
    test_resonant_comparator_rows_record_budget_and_bracket_conditions()
    test_exact_in_e_harmonic_bracket_examples()
    test_exact_in_e_minimum_cutoffs_for_budget_and_bracket()
    test_exact_in_e_design_classification_rows()
    test_amplitude_weighted_exact_in_e_coefficient_formula_and_circular_limit()
    test_amplitude_weighted_design_splits_low_and_moderate_eccentricity()
    test_amplitude_weighted_minimum_cutoff_example()
    test_amplitude_weighted_harmonics_mark_unusable_high_threshold()
    test_physical_detectability_noise_scales_by_harmonic()
    test_physical_detectability_splits_low_and_high_signal_scale()
    test_physical_detectability_minimum_scale_examples()
    test_physical_detectability_default_rows_include_no_go_and_positive()
    test_nonlinear_second_order_transfers_share_denominator()
    test_nonlinear_second_order_readout_sideband_law()
    test_nonlinear_second_order_shared_denominator_residual_vanishes()
    test_nonlinear_second_order_interpolation_obstruction_boundaries()
    test_nonlinear_second_order_budget_and_rows()
    test_nonlinear_detectability_uses_exact_f_squared_coefficients()
    test_nonlinear_detectability_transfer_has_resonant_shape()
    test_nonlinear_detectability_relative_audit_splits_eccentricity()
    test_nonlinear_detectability_snr_minimum_scales()
    test_nonlinear_detectability_default_rows_include_no_go_and_positive()
    test_nonlinear_detectability_generated_rows_mark_unusable_high_cutoff()
    test_component_separability_uses_same_usable_harmonics()
    test_component_separability_jacobian_generated_column_adds_rank()
    test_component_separability_distinguishes_underbudget_and_degenerate_cases()
    test_component_separability_minimum_cutoffs()
    test_component_separability_default_rows_cover_positive_and_negative()
    test_component_separability_joint_observable_is_complex()
    test_component_separability_jacobian_shape()
    test_nonlinear_robustness_row_requires_budget_and_bracket()
    test_nonlinear_robustness_default_grid_counts_regions()
    test_nonlinear_robustness_channel_counts_and_representatives()
    print("symbolic checks passed")


if __name__ == "__main__":
    main()
