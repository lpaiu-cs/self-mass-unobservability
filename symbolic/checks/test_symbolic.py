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
    assert entries["stf_rank4_representative_sector"].candidate_operator_count == 44
    assert entries["stf_rank4_representative_sector"].reduced_operator_count == 19
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

    rank4 = high_rank_audit_summary(4)
    assert rank4.generated_count == 22
    assert rank4.manual_count == 15
    assert rank4.matched_count == 15
    assert rank4.omitted_from_manual == (
        "EEDtQ",
        "EEEQ",
        "EEQ",
        "EQQQ_1",
        "EQQQ_2",
        "GradEGradQ",
        "Q3",
    )
    assert rank4.manual_only == ()

    rank5 = high_rank_audit_summary(5)
    assert rank5.generated_count == 16
    assert rank5.manual_count == 16
    assert rank5.matched_count == 16
    assert rank5.omitted_from_manual == ()
    assert rank5.manual_only == ()


def test_high_rank_diff_report_mentions_eeq() -> None:
    report = high_rank_diff_report()
    assert "\tQ\tomitted_from_manual\tEEQ\t3\tE,E,Q\t" in report
    assert "\tQ\tomitted_from_manual\tQ3\t3\tQ,Q,Q\t" in report
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
    assert summary.total_classes == 44
    assert summary.first_self_witness == "Q2"
    assert summary.first_mixed_witness == "EQQ"
    assert summary.smallest_new_witness == "Q2"
    assert summary.new_surviving_labels == (
        "Q2",
        "EQQ",
        "dotQ2",
        "EQDtQ",
        "E2Q2",
        "E2Q2_mixed_1",
        "E2Q2_mixed_2",
        "E2Q2_mixed_3",
        "divQ2",
        "gradQ2",
        "mixedGradQ2",
        "Q2^2",
        "Q4_bridge",
        "Q4_chain",
        "Q4_tetra",
    )


def test_r4_survivor_rank_check() -> None:
    summary = r4_rank_summary()
    assert summary.rank == 19
    assert summary.count == 22
    assert summary.nullity == 3
    assert summary.new_rank == 12
    assert summary.new_count == 15
    assert "E2Q2_mixed_1" in str(summary.first_null_relation)
    assert "Q4_bridge" in str(summary.additional_null_relations[0])
    assert summary.new_labels == (
        "Q2",
        "EQQ",
        "dotQ2",
        "EQDtQ",
        "E2Q2",
        "E2Q2_mixed_1",
        "E2Q2_mixed_2",
        "E2Q2_mixed_3",
        "divQ2",
        "gradQ2",
        "mixedGradQ2",
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
    print("symbolic checks passed")


if __name__ == "__main__":
    main()
