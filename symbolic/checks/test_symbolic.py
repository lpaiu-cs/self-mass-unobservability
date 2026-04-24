from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

SYMBOLIC_ROOT = Path(__file__).resolve().parents[1]
if str(SYMBOLIC_ROOT) not in sys.path:
    sys.path.insert(0, str(SYMBOLIC_ROOT))

from enumerate_basis import enumerate_minimal_scalar_monomials
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


def main() -> None:
    test_symmetric_quadratic_jet()
    test_worldline_force_structure()
    test_basis_enumeration_order_four()
    print("symbolic checks passed")


if __name__ == "__main__":
    main()
