from __future__ import annotations

import sympy as sp

from sensitivity_expand import SensitivityJet, make_quadratic_jet


def build_worldline_model() -> dict[str, object]:
    x = sp.Symbol("x", real=True)
    y1 = sp.Function("Y1")(x)
    y2 = sp.Function("Y2")(x)
    e2 = sp.Function("E2")(x)
    m0 = sp.Symbol("m0", positive=True)
    lambda_e = sp.Symbol("lambda_E")

    jet = make_quadratic_jet("s", (y1, y2), mass_scale=m0)
    monopole = jet.mass_shift()
    lagrangian = sp.expand(-monopole - sp.Rational(1, 2) * lambda_e * e2)
    force = sp.expand(-sp.diff(lagrangian, x))

    return {
        "x": x,
        "invariants": (y1, y2),
        "tidal_scalar": e2,
        "m0": m0,
        "lambda_E": lambda_e,
        "jet": jet,
        "monopole": monopole,
        "lagrangian": lagrangian,
        "force": force,
    }


def monopole_force(model: dict[str, object]) -> sp.Expr:
    lagrangian = model["lagrangian"]
    x = model["x"]
    lambda_e = model["lambda_E"]
    return sp.expand(-sp.diff(lagrangian.subs(lambda_e, 0), x))


def describe_model(model: dict[str, object]) -> str:
    jet = model["jet"]
    if not isinstance(jet, SensitivityJet):
        raise TypeError("worldline model is missing a SensitivityJet")

    lines = [
        "Monopole mass jet:",
        sp.sstr(model["monopole"]),
        "",
        "Worldline Lagrangian:",
        sp.sstr(model["lagrangian"]),
        "",
        "Free-fall force:",
        sp.sstr(model["force"]),
        "",
        "Sensitivity coordinates:",
        ", ".join(str(symbol) for symbol in jet.coordinates()),
        "Higher-multipole Wilson coefficient:",
        str(model["lambda_E"]),
    ]
    return "\n".join(lines)


def main() -> None:
    model = build_worldline_model()
    print(describe_model(model))


if __name__ == "__main__":
    main()
