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
from family_envelope_census import family_envelope_summary
from family_witness_map import family_witness_summary
from fixed_family_operator_count import fixed_family_operator_count_summary
from high_rank_diff_report import high_rank_diff_report
from high_rank_family_enumerator import high_rank_audit_summary
from irrep_family_census import irrep_family_summary
from nonanalytic_jet_demo import nonanalytic_jet_summary
from hereditary_kernel_demo import hereditary_kernel_summary
from composition_attack_delta4 import composition_summary
from mixed_witness_map import mixed_witness_summary
from primitive_family_attack import primitive_attack_summary
from r1_sector_delta4 import r1_summary
from r1_survivor_rank_check import r1_rank_summary
from r3_sector_delta4 import r3_summary
from r3_survivor_rank_check import r3_rank_summary
from r4_sector_delta4 import r4_summary
from r4_survivor_rank_check import r4_rank_summary
from r5_sector_delta4 import r5_summary
from r5_survivor_rank_check import r5_rank_summary
from r6_sector_delta4 import r6_summary
from r6_survivor_rank_check import r6_rank_summary
from shift_scalar_sector_delta4 import shift_scalar_summary
from stateful_counterexample_demo import stateful_counterexample_summary
from stf_rankL_pattern_check import stf_rank_pattern_summary
from stf_self_witness_check import stf_self_witness_summary
from threshold_formula_check import threshold_formula_summary
from weight_spectrum_demo import weight_spectrum_summary
from witness_threshold_map import witness_threshold_summary
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


def test_family_witness_map() -> None:
    summary = family_witness_summary()
    assert summary.any_harmless is False
    assert tuple(entry.family_class for entry in summary.entries) == (
        "rank2_stf",
        "rank0_scalar_unsuppressed",
        "rank0_scalar_derivative_only",
        "rank1_vector",
        "rank3_tensor_stf",
        "rank4_tensor_stf",
        "rank5_tensor_stf",
        "rank6_tensor_stf",
    )
    assert summary.entries[0].smallest_surviving_operator == "X2"
    assert summary.entries[0].audited_instance == "B2"
    assert summary.entries[1].smallest_surviving_operator == "S"
    assert summary.entries[2].smallest_surviving_operator == "dotS2"
    assert summary.entries[2].weight == 4
    assert summary.entries[3].smallest_surviving_operator == "V2"
    assert summary.entries[3].weight == 2
    assert summary.entries[4].smallest_surviving_operator == "T2"
    assert summary.entries[4].weight == 2
    assert summary.entries[5].smallest_surviving_operator == "Q2"
    assert summary.entries[5].weight == 2
    assert summary.entries[6].smallest_surviving_operator == "U2"
    assert summary.entries[6].weight == 2
    assert summary.entries[7].smallest_surviving_operator == "Z2"
    assert summary.entries[7].weight == 2


def test_witness_threshold_map() -> None:
    summary = witness_threshold_summary()
    assert summary.delta_max == 4
    assert tuple(entry.family_class for entry in summary.entries) == (
        "rank2_stf",
        "rank0_scalar_unsuppressed",
        "rank0_scalar_derivative_only",
        "rank1_vector",
        "rank3_tensor_stf",
        "rank4_tensor_stf",
        "rank5_tensor_stf",
        "rank6_tensor_stf",
    )
    assert summary.entries[0].self_only_lower_bound_at_delta4 == "w_X >= 3"
    assert summary.entries[0].necessary_threshold_at_delta4.startswith("w_X >= 4")
    assert summary.entries[0].threshold_type == "mixed-aware"
    assert summary.entries[1].necessary_threshold_at_delta4.startswith("w_S >= 5")
    assert summary.entries[1].threshold_type == "self-only"
    assert summary.entries[2].necessary_threshold_at_delta4.startswith("w_D >= 3")
    assert summary.entries[2].threshold_type == "tied-sharp"
    assert summary.entries[3].necessary_threshold_at_delta4.startswith("w_V >= 3")
    assert summary.entries[3].threshold_type == "self-only"
    assert summary.entries[4].necessary_threshold_at_delta4.startswith("w_T >= 3")
    assert summary.entries[4].threshold_type == "self-only"
    assert summary.entries[5].necessary_threshold_at_delta4.startswith("w_Q >= 3")
    assert summary.entries[5].threshold_type == "self-only"
    assert summary.entries[6].necessary_threshold_at_delta4.startswith("w_U >= 3")
    assert summary.entries[6].threshold_type == "self-only"
    assert summary.entries[7].necessary_threshold_at_delta4.startswith("w_Z >= 3")
    assert summary.entries[7].threshold_type == "self-only"
    assert all(entry.sufficient_for_uniqueness is False for entry in summary.entries)


def test_mixed_witness_map() -> None:
    summary = mixed_witness_summary()
    assert summary.delta_max == 4
    assert summary.mixed_dominant_family_class is None
    assert tuple(entry.family_class for entry in summary.entries) == (
        "rank2_stf",
        "rank0_scalar_unsuppressed",
        "rank0_scalar_derivative_only",
        "rank1_vector",
        "rank3_tensor_stf",
        "rank4_tensor_stf",
        "rank5_tensor_stf",
        "rank6_tensor_stf",
    )
    rank2, bare_scalar, derivative_only, vector, rank3, rank4, rank5, rank6 = summary.entries
    assert rank2.first_self_witness == "X2"
    assert rank2.first_mixed_witness == "EX"
    assert rank2.self_weight == rank2.mixed_weight == rank2.w_min == 2
    assert "tied" in rank2.sharpness_status
    assert bare_scalar.first_self_witness == "S"
    assert bare_scalar.first_mixed_witness == "SE2"
    assert bare_scalar.self_weight == bare_scalar.w_min == 1
    assert bare_scalar.mixed_weight == 3
    assert bare_scalar.sharpness_status == "self-only"
    assert derivative_only.first_self_witness == "dotS2"
    assert derivative_only.first_mixed_witness == "DtS_E2"
    assert derivative_only.self_weight == derivative_only.mixed_weight == derivative_only.w_min == 4
    assert derivative_only.sharpness_status == "tied-sharp"
    assert vector.first_self_witness == "V2"
    assert vector.first_mixed_witness == "EVV"
    assert vector.self_weight == vector.w_min == 2
    assert vector.mixed_weight == 3
    assert vector.sharpness_status == "self-only"
    assert rank3.first_self_witness == "T2"
    assert rank3.first_mixed_witness == "ETT"
    assert rank3.self_weight == rank3.w_min == 2
    assert rank3.mixed_weight == 3
    assert rank3.sharpness_status == "self-only"
    assert rank4.first_self_witness == "Q2"
    assert rank4.first_mixed_witness == "EEQ and EQQ"
    assert rank4.self_weight == rank4.w_min == 2
    assert rank4.mixed_weight == 3
    assert rank4.sharpness_status == "self-only"
    assert rank5.first_self_witness == "U2"
    assert rank5.first_mixed_witness == "EUU"
    assert rank5.self_weight == rank5.w_min == 2
    assert rank5.mixed_weight == 3
    assert rank5.sharpness_status == "self-only"
    assert rank6.first_self_witness == "Z2"
    assert rank6.first_mixed_witness == "EZZ"
    assert rank6.self_weight == rank6.w_min == 2
    assert rank6.mixed_weight == 3
    assert rank6.sharpness_status == "self-only"


def test_threshold_formula_check() -> None:
    summary = threshold_formula_summary()
    assert summary.delta_max == 4
    assert tuple(entry.family_class for entry in summary.entries) == (
        "rank2_stf",
        "rank0_scalar_unsuppressed",
        "rank0_scalar_derivative_only",
        "rank1_vector",
        "rank3_tensor_stf",
        "rank4_tensor_stf",
        "rank5_tensor_stf",
        "rank6_tensor_stf",
    )
    rank2, bare_scalar, derivative_only, vector, rank3, rank4, rank5, rank6 = summary.entries
    assert rank2.self_formula == "2*w_X"
    assert rank2.mixed_formula == "w_X + 1"
    assert rank2.threshold_type == "mixed-aware"
    assert "tied" in rank2.sharpness_status
    assert bare_scalar.w_min_formula == "w_S"
    assert bare_scalar.threshold_type == "self-only"
    assert derivative_only.self_formula == "2*w_D"
    assert derivative_only.mixed_formula == "w_D + 2"
    assert derivative_only.threshold_type == "tied-sharp"
    assert vector.self_formula == "2*w_V"
    assert vector.mixed_formula == "2*w_V + 1"
    assert vector.w_min_formula == "2*w_V"
    assert vector.threshold_type == "self-only"
    assert rank3.self_formula == "2*w_T"
    assert rank3.mixed_formula == "2*w_T + 1"
    assert rank3.w_min_formula == "2*w_T"
    assert rank3.threshold_type == "self-only"
    assert rank4.self_formula == "2*w_Q"
    assert rank4.mixed_formula == "min(w_Q + 2, 2*w_Q + 1)"
    assert rank4.w_min_formula == "min(2*w_Q, w_Q + 2)"
    assert rank4.threshold_type == "self-only"
    assert rank5.self_formula == "2*w_U"
    assert rank5.mixed_formula == "2*w_U + 1"
    assert rank5.w_min_formula == "2*w_U"
    assert rank5.threshold_type == "self-only"
    assert rank6.self_formula == "2*w_Z"
    assert rank6.mixed_formula == "2*w_Z + 1"
    assert rank6.w_min_formula == "2*w_Z"
    assert rank6.threshold_type == "self-only"


def test_composition_attack_delta4() -> None:
    summary = composition_summary()
    assert summary.baseline_survivors == (
        "E2",
        "E3",
        "dotE2",
        "E2^2",
        "divE2",
        "gradE2",
        "mixedGradE2",
    )
    assert summary.pre_r1_audited_set_joint_sufficient is True
    assert summary.pre_r1_smallest_surviving_cross_family is None
    assert summary.pre_rodd_audited_set_joint_sufficient is True
    assert summary.pre_rodd_smallest_surviving_cross_family is None
    assert summary.pre_reven_audited_set_joint_sufficient is True
    assert summary.pre_reven_smallest_surviving_cross_family is None
    assert summary.pre_rodd5_audited_set_joint_sufficient is True
    assert summary.pre_rodd5_smallest_surviving_cross_family is None
    assert summary.pre_reven6_audited_set_joint_sufficient is True
    assert summary.pre_reven6_smallest_surviving_cross_family is None
    assert summary.post_reven6_audited_set_joint_sufficient is True
    assert summary.post_reven6_smallest_surviving_cross_family is None
    assert len(summary.cases) == 247
    scope_counts = {
        "pre_r1_audited_set": 0,
        "pre_rodd_audited_set": 0,
        "pre_reven_audited_set": 0,
        "pre_rodd5_audited_set": 0,
        "pre_reven6_audited_set": 0,
        "post_reven6_audited_set": 0,
    }
    kind_counts = {
        "pairwise": 0,
        "triple": 0,
        "quadruple": 0,
        "quintuple": 0,
        "six-family": 0,
        "seven-family": 0,
        "full-set": 0,
    }
    case_by_id = {case.combo_id: case for case in summary.cases}
    for case in summary.cases:
        scope_counts[case.case_scope] += 1
        kind_counts[case.combination_kind] += 1
    assert scope_counts == {
        "pre_r1_audited_set": 4,
        "pre_rodd_audited_set": 7,
        "pre_reven_audited_set": 15,
        "pre_rodd5_audited_set": 31,
        "pre_reven6_audited_set": 63,
        "post_reven6_audited_set": 127,
    }
    assert kind_counts == {
        "pairwise": 28,
        "triple": 56,
        "quadruple": 70,
        "quintuple": 55,
        "six-family": 27,
        "seven-family": 7,
        "full-set": 4,
    }
    assert case_by_id["Reven4++R2"].case_scope == "pre_rodd5_audited_set"
    assert case_by_id["Reven4++R2"].combination_kind == "pairwise"
    assert case_by_id["Reven4++R2+R0a"].combination_kind == "triple"
    assert case_by_id["Reven4++R2+R0a+R0b"].combination_kind == "quadruple"
    assert case_by_id["Reven4++R2+R0a+R0b+R1"].combination_kind == "quintuple"
    assert (
        case_by_id["R2+R0a+R0b+R1+Rodd++Reven4+"].combination_kind
        == "full-set"
    )
    assert case_by_id["Rodd5++R2"].case_scope == "pre_reven6_audited_set"
    assert case_by_id["Rodd5++R2"].combination_kind == "pairwise"
    assert case_by_id["Rodd5++R2+R0a"].combination_kind == "triple"
    assert case_by_id["Rodd5++R2+R0a+R0b"].combination_kind == "quadruple"
    assert case_by_id["Rodd5++R2+R0a+R0b+R1"].combination_kind == "quintuple"
    assert case_by_id["Rodd5++R2+R0a+R0b+R1+Rodd+"].combination_kind == "six-family"
    assert (
        case_by_id["R2+R0a+R0b+R1+Rodd++Reven4++Rodd5+"].combination_kind
        == "full-set"
    )
    assert case_by_id["Reven6++R2"].case_scope == "post_reven6_audited_set"
    assert case_by_id["Reven6++R2"].combination_kind == "pairwise"
    assert case_by_id["Reven6++R2+R0a"].combination_kind == "triple"
    assert case_by_id["Reven6++R2+R0a+R0b"].combination_kind == "quadruple"
    assert case_by_id["Reven6++R2+R0a+R0b+R1"].combination_kind == "quintuple"
    assert case_by_id["Reven6++R2+R0a+R0b+R1+Rodd+"].combination_kind == "six-family"
    assert (
        case_by_id["Reven6++R2+R0a+R0b+R1+Rodd++Reven4+"].combination_kind
        == "seven-family"
    )
    assert (
        case_by_id["R2+R0a+R0b+R1+Rodd++Reven4++Rodd5++Reven6+"].combination_kind
        == "full-set"
    )
    assert all(case.sufficient for case in summary.cases)
    assert all(case.new_surviving_labels == () for case in summary.cases)
    assert all(case.smallest_surviving_cross_family is None for case in summary.cases)
    assert all(case.smallest_surviving_cross_family_weight is None for case in summary.cases)


def test_family_envelope_census() -> None:
    summary = family_envelope_summary()
    assert summary.delta_max == 4
    assert (
        summary.package_status
        == "irreducible scalar/vector/STF family-envelope closure is resolved positively inside the current theorem domain"
    )
    assert summary.envelope_closed is True
    assert summary.first_open_class is None
    entries = {entry.class_id: entry for entry in summary.entries}
    assert entries["Scalar"].envelope_state == "audited"
    assert entries["Scalar"].smallest_expected_witness_type == "S / dotS2"
    assert entries["Vector"].envelope_state == "audited"
    assert entries["Vector"].smallest_expected_witness_type == "V2 / EVV"
    assert entries["STF2"].envelope_state == "audited"
    assert entries["STF2"].smallest_expected_witness_type == "X2 / EX"
    assert entries["STFge3"].envelope_state == "audited"
    assert entries["STFge3"].smallest_expected_witness_type == "Y2 / class-limited mixed pattern"
    assert entries["TraceDesc"].envelope_state == "absorbed by trace reduction"
    assert entries["MixedEvenDual"].envelope_state == "absorbed by trace reduction"
    assert entries["PseudoOdd"].envelope_state == "excluded by explicit assumption"
    assert entries["State"].envelope_state == "excluded by explicit assumption"
    assert entries["Nonlocal"].envelope_state == "excluded by explicit assumption"


def test_irrep_family_census() -> None:
    summary = irrep_family_summary()
    assert summary.delta_max == 4
    assert (
        summary.package_status
        == "irreducible scalar/vector/STF family-envelope closure is resolved positively inside the current theorem domain"
    )
    assert summary.envelope_closes_on_audited_classes is True
    assert summary.first_open_nonstf_family is None
    entries = {entry.class_id: entry for entry in summary.entries}
    assert entries["Scalar"].resolution_state == "already audited"
    assert entries["Vector"].resolution_state == "already audited"
    assert entries["STF2"].resolution_state == "already audited"
    assert entries["STFge3"].resolution_state == "already audited"
    assert entries["TraceDesc"].resolution_state == "absorbed by trace reduction"
    assert entries["MixedEvenDual"].resolution_state == "absorbed by trace reduction"
    assert entries["PseudoOdd"].resolution_state == "excluded by parity/nonspin assumptions"


def test_fixed_family_operator_count_summary() -> None:
    summary = fixed_family_operator_count_summary()
    assert summary.delta_max == 4
    assert summary.theorem_domain_classes == ("Scalar", "Vector", "STF2", "STFge3")
    assert summary.baseline_sector == "electric parity-even free-fall scalar sector"
    assert summary.envelope_closed is True
    assert summary.candidate_operator_space_finite is True
    assert summary.reduced_operator_space_finite is True
    assert summary.reduction_layers_applied == (
        "irreducible family-envelope closure",
        "trace-descendant absorption",
        "total derivatives",
        "lower-order EOM",
        "explicit algebraic identities",
        "sector-specific linear-dependence quotient when present",
    )
    entries = {entry.catalog_id: entry for entry in summary.entries}
    assert entries["electric_exact_current_set"].candidate_operator_count == 21
    assert entries["electric_exact_current_set"].reduced_operator_count == 7
    assert entries["rank2_stf_special_case"].candidate_operator_count == 42
    assert entries["rank2_stf_special_case"].reduced_operator_count == 18
    assert entries["scalar_bare_source_extension"].candidate_operator_count == 72
    assert entries["scalar_bare_source_extension"].reduced_operator_count == 33
    assert entries["scalar_derivative_only_extension"].candidate_operator_count == 52
    assert entries["scalar_derivative_only_extension"].reduced_operator_count == 23
    assert entries["vector_representative_sector"].candidate_operator_count == 39
    assert entries["vector_representative_sector"].reduced_operator_count == 17
    assert entries["stf_rank3_representative_sector"].candidate_operator_count == 43
    assert entries["stf_rank3_representative_sector"].reduced_operator_count == 19
    assert entries["stf_rank4_representative_sector"].candidate_operator_count == 50
    assert entries["stf_rank4_representative_sector"].reduced_operator_count == 25
    assert entries["stf_rank5_representative_sector"].candidate_operator_count == 45
    assert entries["stf_rank5_representative_sector"].reduced_operator_count == 19
    assert entries["stf_rank6_representative_sector"].candidate_operator_count == 43
    assert entries["stf_rank6_representative_sector"].reduced_operator_count == 23


def test_nonanalytic_jet_demo() -> None:
    summary = nonanalytic_jet_summary()
    assert summary.delta_max == 4
    assert summary.locality_kept_in_all_cases is True
    assert summary.finite_family_operator_closure_kept_in_all_cases is True
    assert summary.smallest_local_nonanalytic_counterexample == "smooth_flat_single_coordinate"
    assert summary.broken_layer == "analytic monopole jet collapse (Lemma 55 / A5)"
    cases = {case.case_id: case for case in summary.cases}
    analytic = cases["analytic_control_quadratic"]
    assert analytic.finite_taylor_jet_valid is True
    assert analytic.theorem_layer_broken is None
    assert analytic.jet_value_at_sample == analytic.model_value_at_sample
    flat = cases["smooth_flat_single_coordinate"]
    assert flat.locality_kept is True
    assert flat.finite_family_operator_closure_kept is True
    assert flat.finite_taylor_jet_valid is False
    assert flat.canonical_counterexample is True
    assert "all Taylor coefficients" in (flat.first_exact_failure_mode or "")
    assert flat.jet_value_at_sample == 1.0
    assert flat.model_value_at_sample > flat.jet_value_at_sample
    threshold = cases["threshold_sqrt_activation"]
    assert threshold.locality_kept is True
    assert threshold.finite_family_operator_closure_kept is True
    assert threshold.finite_taylor_jet_valid is False
    assert "branch point" in (threshold.first_exact_failure_mode or "")
    assert threshold.jet_value_at_sample is None


def test_hereditary_kernel_demo() -> None:
    summary = hereditary_kernel_summary()
    assert summary.delta_max == 4
    assert summary.finite_primitive_family_envelope_kept_in_all_cases is True
    assert summary.analyticity_kept_in_all_cases is True
    assert summary.smallest_genuinely_hereditary_counterexample == "power_law_single_coordinate"
    assert "A3" in summary.broken_layer
    cases = {case.case_id: case for case in summary.cases}
    local_control = cases["local_analytic_control"]
    assert local_control.locality_kept is True
    assert local_control.finite_local_jet_valid is True
    assert local_control.theorem_layer_first_broken is None
    assert local_control.same_point_response_difference == 0.0
    exponential = cases["exponential_memory_control"]
    assert exponential.locality_kept is False
    assert exponential.analyticity_kept is True
    assert exponential.finite_local_jet_valid is False
    assert exponential.finite_state_markovianizable is True
    assert exponential.same_point_response_difference > 0.0
    assert "A4-type extension" in (exponential.theorem_layer_first_broken or "")
    power_law = cases["power_law_single_coordinate"]
    assert power_law.locality_kept is False
    assert power_law.analyticity_kept is True
    assert power_law.finite_local_jet_valid is False
    assert power_law.finite_state_markovianizable is False
    assert power_law.canonical_counterexample is True
    assert power_law.same_point_response_difference > 0.0
    assert "Lemmas 55 and 56" in (power_law.theorem_layer_first_broken or "")


def test_stateful_counterexample_demo() -> None:
    summary = stateful_counterexample_summary()
    assert summary.delta_max == 4
    assert summary.locality_kept_in_all_cases is True
    assert summary.analyticity_kept_in_all_cases is True
    assert summary.finite_primitive_family_envelope_kept_in_all_cases is True
    assert summary.smallest_local_finite_state_counterexample == "dynamical_one_state_chi"
    assert "A4" in summary.broken_layer
    assert summary.finite_state_augmented_collapse_survives is True
    cases = {case.case_id: case for case in summary.cases}
    local_control = cases["local_analytic_no_state_control"]
    assert local_control.y_only_finite_jet_valid is True
    assert local_control.finite_state_augmented_bookkeeping_valid is True
    assert local_control.theorem_layer_first_broken is None
    assert local_control.same_point_response_difference == 0.0
    slaved = cases["adiabatic_slaved_local_state"]
    assert slaved.y_only_finite_jet_valid is True
    assert slaved.finite_state_augmented_bookkeeping_valid is True
    assert slaved.theorem_layer_first_broken is None
    assert slaved.same_point_response_difference == 0.0
    dynamic = cases["dynamical_one_state_chi"]
    assert dynamic.locality_kept is True
    assert dynamic.analyticity_kept is True
    assert dynamic.y_only_finite_jet_valid is False
    assert dynamic.finite_state_augmented_bookkeeping_valid is True
    assert dynamic.canonical_counterexample is True
    assert dynamic.same_point_response_difference > 0.0
    assert "Y-only form" in (dynamic.theorem_layer_first_broken or "")


def test_weight_spectrum_demo() -> None:
    summary = weight_spectrum_summary()
    assert summary.delta_max == 4
    assert summary.local_weight_spectrum_finiteness_suffices is True
    assert summary.a8_stronger_than_necessary is True
    assert summary.sharp_a8_failure_mode == "infinite_low_weight_stf_tower"
    cases = {case.case_id: case for case in summary.cases}
    finite_control = cases["finite_family_catalog_control"]
    assert finite_control.finite_primitive_family_content is True
    assert finite_control.locally_finite_below_delta_max is True
    assert finite_control.candidate_operator_space_finite is True
    assert finite_control.normal_form_quotient_finite is True
    assert finite_control.theorem_layer_first_broken is None
    local_infinite = cases["infinite_but_locally_finite_weight_spectrum"]
    assert local_infinite.finite_primitive_family_content is False
    assert local_infinite.locally_finite_below_delta_max is True
    assert local_infinite.candidate_operator_space_finite is True
    assert local_infinite.normal_form_quotient_finite is True
    assert local_infinite.theorem_layer_first_broken is None
    low_weight = cases["infinite_low_weight_stf_tower"]
    assert low_weight.finite_primitive_family_content is False
    assert low_weight.locally_finite_below_delta_max is False
    assert low_weight.candidate_operator_space_finite is False
    assert low_weight.normal_form_quotient_finite is False
    assert "Lemma 53" in (low_weight.theorem_layer_first_broken or "")


def test_high_rank_exhaustiveness_audit() -> None:
    rank3 = high_rank_audit_summary(3)
    assert rank3.generated_count == 14
    assert rank3.manual_count == 14
    assert rank3.matched_count == 14
    assert rank3.omitted_from_manual == ()
    assert rank3.manual_only == ()

    # Post-completion truth (r4 survivor list completed to 25): the audit's
    # historical omission set is closed. EEDtQ trades against EDtEQ through
    # the total derivative Dt(EEQ) once DtE is in the operator set, and only
    # one EQQQ contraction class survives reduction (EQ3); EDtEQ itself uses
    # the DtE block, which lies outside this enumerator's candidate
    # universe, so it reports as manual_only.
    rank4 = high_rank_audit_summary(4)
    assert rank4.generated_count == 22
    assert rank4.manual_count == 21
    assert rank4.matched_count == 20
    assert rank4.omitted_from_manual == (
        "EEDtQ",
        "EQQQ_1",
    )
    assert rank4.manual_only == ("EDtEQ",)

    rank5 = high_rank_audit_summary(5)
    assert rank5.generated_count == 16
    assert rank5.manual_count == 16
    assert rank5.matched_count == 16
    assert rank5.omitted_from_manual == ()
    assert rank5.manual_only == ()


def test_high_rank_diff_report_mentions_eeq() -> None:
    # EEQ (and QQQ) were the audit's original omission finds; after the
    # sector completion they must report as matched, the residual
    # generated-only classes are EEDtQ/EQQQ_1, and the out-of-universe
    # EDtEQ surfaces as manual_only.
    report = high_rank_diff_report()
    assert "\tQ\tmatched\tEEQ\t3\tE,E,Q\t" in report
    assert "\tQ\tmatched\tQQQ\t3\tQ,Q,Q\t" in report
    assert "\tQ\tomitted_from_manual\tEEDtQ\t4\tE,E,DtQ\t" in report
    assert "\tQ\tmanual_only\tEDtEQ\t4\tE,DtE,Q\tout-of-universe" in report
    assert "\tT\tmatched\tETDtT\t4\tE,T,DtT\t" in report
    assert "\tU\tmatched\tEUDtU\t4\tE,U,DtU\t" in report


def test_stf_rank_pattern_summary() -> None:
    summary = stf_rank_pattern_summary()
    assert summary.attempted_tower_theorem_holds is False
    assert summary.failure_rank == 4
    assert "EEQ" in (summary.failure_reason or "")
    assert tuple(entry.rank for entry in summary.entries) == (3, 4, 5, 6)
    assert summary.entries[0].first_self_witness == "T2"
    assert summary.entries[0].first_mixed_witness_layer == "ETT"
    assert summary.entries[1].first_self_witness == "Q2"
    assert summary.entries[1].first_mixed_witness_layer == "EEQ and EQQ"
    assert "exception" in summary.entries[1].current_audited_status
    assert summary.entries[2].first_self_witness == "U2"
    assert summary.entries[2].first_mixed_witness_layer == "EUU"
    assert summary.entries[3].first_self_witness == "Z2"
    assert summary.entries[3].first_mixed_witness_layer == "EZZ"


def test_stf_self_witness_summary() -> None:
    summary = stf_self_witness_summary()
    assert summary.theorem_holds is True
    assert summary.failure_rank is None
    assert summary.failure_reason is None
    assert tuple(entry.rank for entry in summary.entries) == (3, 4, 5, 6)
    for entry in summary.entries:
        assert entry.linear_scalar_exists is False
        assert entry.mixed_quadratic_witness_exists is False
        assert entry.threshold_at_delta4 == "w_Y >= 3"
    assert summary.entries[0].first_self_witness == "T2"
    assert summary.entries[0].first_mixed_witness_layer == "ETT"
    assert summary.entries[1].first_self_witness == "Q2"
    assert summary.entries[1].first_mixed_witness_layer == "EEQ and EQQ"
    assert summary.entries[2].first_self_witness == "U2"
    assert summary.entries[2].first_mixed_witness_layer == "EUU"
    assert summary.entries[3].first_self_witness == "Z2"
    assert summary.entries[3].first_mixed_witness_layer == "EZZ"


def test_r1_sector_delta4() -> None:
    summary = r1_summary()
    assert summary.total_classes == 39
    assert summary.first_self_witness == "V2"
    assert summary.first_mixed_witness == "EVV"
    assert summary.smallest_new_witness == "V2"
    assert summary.new_surviving_labels == (
        "V2",
        "EVV",
        "dotV2",
        "EVDtV",
        "E2V2",
        "EVEV",
        "divV2",
        "gradV2",
        "mixedGradV2",
        "V2^2",
    )


def test_r1_survivor_rank_check() -> None:
    summary = r1_rank_summary()
    assert summary.rank == 17
    assert summary.count == 17
    assert summary.nullity == 0
    assert summary.new_labels == (
        "V2",
        "EVV",
        "dotV2",
        "EVDtV",
        "E2V2",
        "EVEV",
        "divV2",
        "gradV2",
        "mixedGradV2",
        "V2^2",
    )


def test_r3_sector_delta4() -> None:
    summary = r3_summary()
    assert summary.total_classes == 43
    assert summary.first_self_witness == "T2"
    assert summary.first_mixed_witness == "ETT"
    assert summary.smallest_new_witness == "T2"
    assert summary.new_surviving_labels == (
        "T2",
        "ETT",
        "dotT2",
        "ETDtT",
        "E2T2",
        "E2T2_mixed_1",
        "E2T2_mixed_2",
        "E2T2_mixed_3",
        "divT2",
        "gradT2",
        "mixedGradT2",
        "T2^2",
        "T4_chain",
        "T4_tetra",
    )


def test_r3_survivor_rank_check() -> None:
    summary = r3_rank_summary()
    assert summary.rank == 19
    assert summary.count == 21
    assert summary.nullity == 2
    assert "E2T2_mixed_1" in str(summary.null_relation)
    assert "E2T2_mixed_2" in str(summary.null_relation)
    assert summary.new_labels == (
        "T2",
        "ETT",
        "dotT2",
        "ETDtT",
        "E2T2",
        "E2T2_mixed_1",
        "E2T2_mixed_2",
        "E2T2_mixed_3",
        "divT2",
        "gradT2",
        "mixedGradT2",
        "T2^2",
        "T4_chain",
        "T4_tetra",
    )


def test_r4_sector_delta4() -> None:
    summary = r4_summary()
    assert summary.total_classes == 50
    assert summary.first_self_witness == "Q2"
    assert summary.first_mixed_witness == "EEQ"
    assert summary.smallest_new_witness == "Q2"
    assert summary.new_surviving_labels == (
        "Q2",
        "EEQ",
        "EQQ",
        "QQQ",
        "dotQ2",
        "EQDtQ",
        "EDtEQ",
        "E2Q2",
        "E2Q2_mixed_1",
        "E2Q2_mixed_2",
        "E2Q2_mixed_3",
        "EQ3",
        "GradEGradQ",
        "divQ2",
        "gradQ2",
        "mixedGradQ2",
        "E3Q",
        "Q2^2",
        "Q4_bridge",
        "Q4_chain",
        "Q4_tetra",
    )


def test_r4_survivor_rank_check() -> None:
    summary = r4_rank_summary()
    assert summary.rank == 25
    assert summary.count == 28
    assert summary.nullity == 3
    assert summary.new_rank == 18
    assert summary.new_count == 21
    assert "E2Q2_mixed_1" in str(summary.first_null_relation)
    assert "Q4_bridge" in str(summary.additional_null_relations[0])
    assert summary.new_labels == (
        "Q2",
        "EEQ",
        "EQQ",
        "QQQ",
        "dotQ2",
        "EQDtQ",
        "EDtEQ",
        "E2Q2",
        "E2Q2_mixed_1",
        "E2Q2_mixed_2",
        "E2Q2_mixed_3",
        "EQ3",
        "GradEGradQ",
        "divQ2",
        "gradQ2",
        "mixedGradQ2",
        "E3Q",
        "Q2^2",
        "Q4_bridge",
        "Q4_chain",
        "Q4_tetra",
    )


def test_r5_sector_delta4() -> None:
    summary = r5_summary()
    assert summary.total_classes == 45
    assert summary.first_self_witness == "U2"
    assert summary.first_mixed_witness == "EUU"
    assert summary.smallest_new_witness == "U2"
    assert summary.new_surviving_labels == (
        "U2",
        "EUU",
        "dotU2",
        "EUDtU",
        "E2U2",
        "E2U2_mixed_1",
        "E2U2_mixed_2",
        "E2U2_mixed_3",
        "divU2",
        "gradU2",
        "mixedGradU2",
        "U2^2",
        "U4_balanced",
        "U4_bridge",
        "U4_chain",
        "U4_tetra",
    )


def test_r5_survivor_rank_check() -> None:
    summary = r5_rank_summary()
    assert summary.rank == 19
    assert summary.count == 23
    assert summary.nullity == 4
    assert summary.new_rank == 12
    assert summary.new_count == 16
    assert "E2U2_mixed_1" in str(summary.first_null_relation)
    assert "U4_bridge" in str(summary.additional_null_relations[0])
    assert summary.new_labels == (
        "U2",
        "EUU",
        "dotU2",
        "EUDtU",
        "E2U2",
        "E2U2_mixed_1",
        "E2U2_mixed_2",
        "E2U2_mixed_3",
        "divU2",
        "gradU2",
        "mixedGradU2",
        "U2^2",
        "U4_balanced",
        "U4_bridge",
        "U4_chain",
        "U4_tetra",
    )


def test_r6_sector_delta4() -> None:
    summary = r6_summary()
    assert summary.total_classes == 43
    assert summary.total_new_family_classes == 22
    assert summary.first_self_witness == "Z2"
    assert summary.first_mixed_witness == "EZZ"
    assert summary.first_mixed_layer_labels == ("EZZ",)
    assert summary.even_rank_additional_first_mixed_labels == ()
    assert summary.smallest_new_witness == "Z2"
    assert summary.new_surviving_labels == (
        "Z2",
        "EZZ",
        "Z3",
        "dotZ2",
        "EEEZ",
        "E2Z2",
        "E2Z2_mixed_1",
        "E2Z2_mixed_2",
        "E2Z2_mixed_3",
        "EZDtZ",
        "EZZZ_1",
        "EZZZ_2",
        "divZ2",
        "gradZ2",
        "mixedGradZ2",
        "Z4_1",
        "Z4_2",
        "Z4_3",
        "Z4_4",
        "Z4_5",
        "Z4_6",
        "Z4_7",
    )


def test_r6_survivor_rank_check() -> None:
    summary = r6_rank_summary()
    assert summary.rank == 23
    assert summary.count == 29
    assert summary.nullity == 6
    assert summary.new_rank == 16
    assert summary.new_count == 22
    assert "E2Z2_mixed_1" in str(summary.first_null_relation)
    signature_ranks = {group.signature: group for group in summary.signature_ranks}
    assert signature_ranks[("E", "E", "Z", "Z")].rank == 3
    assert signature_ranks[("E", "E", "Z", "Z")].count == 4
    assert signature_ranks[("E", "Z", "Z", "Z")].rank == 1
    assert signature_ranks[("Z", "Z", "Z", "Z")].rank == 3
    assert summary.new_labels == (
        "Z2",
        "EZZ",
        "Z3",
        "dotZ2",
        "EEEZ",
        "E2Z2",
        "E2Z2_mixed_1",
        "E2Z2_mixed_2",
        "E2Z2_mixed_3",
        "EZDtZ",
        "EZZZ_1",
        "EZZZ_2",
        "divZ2",
        "gradZ2",
        "mixedGradZ2",
        "Z4_1",
        "Z4_2",
        "Z4_3",
        "Z4_4",
        "Z4_5",
        "Z4_6",
        "Z4_7",
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
    test_family_witness_map()
    test_witness_threshold_map()
    test_mixed_witness_map()
    test_threshold_formula_check()
    test_composition_attack_delta4()
    test_family_envelope_census()
    test_irrep_family_census()
    test_fixed_family_operator_count_summary()
    test_nonanalytic_jet_demo()
    test_hereditary_kernel_demo()
    test_stateful_counterexample_demo()
    test_weight_spectrum_demo()
    test_high_rank_exhaustiveness_audit()
    test_high_rank_diff_report_mentions_eeq()
    test_stf_rank_pattern_summary()
    test_stf_self_witness_summary()
    test_r1_sector_delta4()
    test_r1_survivor_rank_check()
    test_r3_sector_delta4()
    test_r3_survivor_rank_check()
    test_r4_sector_delta4()
    test_r4_survivor_rank_check()
    test_r5_sector_delta4()
    test_r5_survivor_rank_check()
    test_r6_sector_delta4()
    test_r6_survivor_rank_check()
    print("symbolic checks passed (theorem + dynamic-chi)")


if __name__ == "__main__":
    main()
