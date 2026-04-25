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
from composition_attack_delta4 import composition_summary
from mixed_witness_map import mixed_witness_summary
from primitive_family_attack import primitive_attack_summary
from r1_sector_delta4 import r1_summary
from r1_survivor_rank_check import r1_rank_summary
from shift_scalar_sector_delta4 import shift_scalar_summary
from threshold_formula_check import threshold_formula_summary
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
    )
    assert summary.entries[0].smallest_surviving_operator == "X2"
    assert summary.entries[0].audited_instance == "B2"
    assert summary.entries[1].smallest_surviving_operator == "S"
    assert summary.entries[2].smallest_surviving_operator == "dotS2"
    assert summary.entries[2].weight == 4
    assert summary.entries[3].smallest_surviving_operator == "V2"
    assert summary.entries[3].weight == 2


def test_witness_threshold_map() -> None:
    summary = witness_threshold_summary()
    assert summary.delta_max == 4
    assert tuple(entry.family_class for entry in summary.entries) == (
        "rank2_stf",
        "rank0_scalar_unsuppressed",
        "rank0_scalar_derivative_only",
        "rank1_vector",
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
    )
    rank2, bare_scalar, derivative_only, vector = summary.entries
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


def test_threshold_formula_check() -> None:
    summary = threshold_formula_summary()
    assert summary.delta_max == 4
    assert tuple(entry.family_class for entry in summary.entries) == (
        "rank2_stf",
        "rank0_scalar_unsuppressed",
        "rank0_scalar_derivative_only",
        "rank1_vector",
    )
    rank2, bare_scalar, derivative_only, vector = summary.entries
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
    assert tuple(case.combo_id for case in summary.cases) == (
        "R2+R0a",
        "R2+R0b",
        "R0a+R0b",
        "R2+R0a+R0b",
    )
    assert all(case.sufficient for case in summary.cases)
    assert all(case.new_surviving_labels == () for case in summary.cases)
    assert all(case.smallest_surviving_cross_family is None for case in summary.cases)


def test_family_envelope_census() -> None:
    summary = family_envelope_summary()
    assert summary.delta_max == 4
    assert (
        summary.live_bottleneck
        == "family-envelope completeness or the next smallest unaudited family obstruction"
    )
    assert summary.envelope_closed is False
    assert summary.smallest_unaudited_class == "Rodd+"
    entries = {entry.class_id: entry for entry in summary.entries}
    assert entries["R2"].envelope_state == "audited"
    assert entries["R0a"].envelope_state == "audited"
    assert entries["R0b"].envelope_state == "audited"
    assert entries["R1"].envelope_state == "audited"
    assert entries["R1"].smallest_expected_witness_type == "V2 / EVV"
    assert entries["Rodd+"].envelope_state == "still unaudited"
    assert entries["Podd"].envelope_state == "excluded by explicit assumption"


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
    test_r1_sector_delta4()
    test_r1_survivor_rank_check()
    print("symbolic checks passed")


if __name__ == "__main__":
    main()
