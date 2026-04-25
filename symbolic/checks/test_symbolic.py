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
    print("symbolic checks passed")


if __name__ == "__main__":
    main()
