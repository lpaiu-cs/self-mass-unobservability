from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class ReductionExample:
    label: str
    weight: int
    reduction_channel: str
    expression: sp.Expr


def operator_symbols() -> dict[str, sp.Expr]:
    names = (
        "E2",
        "E3",
        "E4",
        "dotE2",
        "gradE2",
        "E_DtE",
        "Dt2_E2",
        "E_Dt2E",
        "a2",
        "aEa",
    )
    return {name: sp.Symbol(name) for name in names}


def reduction_examples() -> tuple[ReductionExample, ...]:
    ops = operator_symbols()
    return (
        ReductionExample("E2", 2, "already normal form", ops["E2"]),
        ReductionExample("E3", 3, "already normal form", ops["E3"]),
        ReductionExample("E2^2", 4, "already normal form", ops["E2"] ** 2),
        ReductionExample("dotE2", 4, "already normal form", ops["dotE2"]),
        ReductionExample("gradE2", 4, "already normal form", ops["gradE2"]),
        ReductionExample("E_DtE", 3, "total derivative", ops["E_DtE"]),
        ReductionExample("Dt2_E2", 4, "total derivative", ops["Dt2_E2"]),
        ReductionExample("E_Dt2E", 4, "total derivative", ops["E_Dt2E"]),
        ReductionExample("E4", 4, "algebraic identity", ops["E4"]),
        ReductionExample("a2", 2, "lower-order EOM", ops["a2"]),
        ReductionExample("aEa", 3, "lower-order EOM", ops["aEa"]),
    )


def reduce_total_derivatives(expr: sp.Expr) -> sp.Expr:
    ops = operator_symbols()
    rules = {
        ops["E_DtE"]: sp.Integer(0),
        ops["Dt2_E2"]: sp.Integer(0),
        ops["E_Dt2E"]: -ops["dotE2"],
    }
    return sp.expand(expr.subs(rules))


def reduce_lower_order_eom(expr: sp.Expr) -> sp.Expr:
    ops = operator_symbols()
    rules = {
        ops["a2"]: sp.Integer(0),
        ops["aEa"]: sp.Integer(0),
    }
    return sp.expand(expr.subs(rules))


def reduce_algebraic_identities(expr: sp.Expr) -> sp.Expr:
    ops = operator_symbols()
    rules = {
        ops["E4"]: sp.Rational(1, 2) * ops["E2"] ** 2,
    }
    return sp.expand(expr.subs(rules))


def reduce_to_normal_form(expr: sp.Expr) -> sp.Expr:
    reduced = reduce_total_derivatives(expr)
    reduced = reduce_lower_order_eom(reduced)
    reduced = reduce_algebraic_identities(reduced)
    return sp.expand(reduced)


def normal_form_basis() -> tuple[sp.Expr, ...]:
    ops = operator_symbols()
    return (
        ops["E2"],
        ops["E3"],
        ops["E2"] ** 2,
        ops["dotE2"],
        ops["gradE2"],
    )


def reduction_report() -> str:
    lines = [
        "Delta<=4 normal-form reduction table",
        "",
        "Normal-form target basis:",
        "- " + ", ".join(sp.sstr(term) for term in normal_form_basis()),
        "",
        "Catalog reductions:",
    ]

    for example in reduction_examples():
        reduced = reduce_to_normal_form(example.expression)
        lines.append(
            f"- {example.label} [weight {example.weight}; {example.reduction_channel}]"
            f" -> {sp.sstr(reduced)}"
        )

    lines.extend(
        [
            "",
            "Status note:",
            "- No explicit obstruction was found inside the exact Delta<=4 catalog handled here.",
            "- This script does not prove catalog exhaustiveness.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(reduction_report())


if __name__ == "__main__":
    main()
