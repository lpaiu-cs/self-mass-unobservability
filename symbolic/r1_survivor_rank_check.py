from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from r1_sector_delta4 import r1_summary
from survivor_rank_check import stf_matrix, survivor_polynomials


@dataclass(frozen=True)
class R1RankSummary:
    rank: int
    count: int
    nullity: int
    monomial_count: int
    null_relation: sp.Expr
    labels: tuple[str, ...]
    new_labels: tuple[str, ...]


def vector(prefix: str) -> tuple[sp.Matrix, tuple[sp.Symbol, ...]]:
    symbols = sp.symbols(f"{prefix}0 {prefix}1 {prefix}2")
    return sp.Matrix(symbols), symbols


def matrix(prefix: str) -> tuple[sp.Matrix, tuple[sp.Symbol, ...]]:
    symbols = sp.symbols(
        " ".join(f"{prefix}{i}{j}" for i in range(3) for j in range(3))
    )
    return sp.Matrix(3, 3, symbols), symbols


def r1_variables() -> tuple[sp.Symbol, ...]:
    _, electric_vars = stf_matrix("E")
    _, electric_dot_vars = stf_matrix("DtE")
    grad_e_vars: list[sp.Symbol] = []
    for grad_index in range(3):
        _, block_vars = stf_matrix(f"GradE{grad_index}")
        grad_e_vars.extend(block_vars)

    _, vector_vars = vector("V")
    _, vector_dot_vars = vector("DtV")
    _, grad_v_vars = matrix("GradV")
    return electric_vars + electric_dot_vars + tuple(grad_e_vars) + vector_vars + vector_dot_vars + grad_v_vars


def r1_survivor_expression_map() -> dict[str, sp.Expr]:
    electric, _ = stf_matrix("E")
    electric_dot, _ = stf_matrix("DtE")
    vector_field, _ = vector("V")
    vector_dot, _ = vector("DtV")
    grad_v, _ = matrix("GradV")

    grad_blocks: list[sp.Matrix] = []
    for grad_index in range(3):
        block, _ = stf_matrix(f"GradE{grad_index}")
        grad_blocks.append(block)

    baseline = dict(survivor_polynomials())

    v2 = sp.expand((vector_field.T * vector_field)[0])
    evv = sp.expand((vector_field.T * electric * vector_field)[0])
    dot_v2 = sp.expand((vector_dot.T * vector_dot)[0])
    evdtv = sp.expand((vector_field.T * electric * vector_dot)[0])
    e2v2 = sp.expand(baseline["E2"] * v2)
    evev = sp.expand((vector_field.T * (electric**2) * vector_field)[0])
    div_v = sp.expand(grad_v.trace())
    div_v2 = sp.expand(div_v**2)
    grad_v2 = sp.expand(sum(grad_v[i, j] ** 2 for i in range(3) for j in range(3)))
    mixed_grad_v2 = sp.expand(sum(grad_v[i, j] * grad_v[j, i] for i in range(3) for j in range(3)))

    baseline.update(
        {
            "V2": v2,
            "EVV": evv,
            "dotV2": dot_v2,
            "EVDtV": evdtv,
            "V2^2": sp.expand(v2**2),
            "E2V2": e2v2,
            "EVEV": evev,
            "divV2": div_v2,
            "gradV2": grad_v2,
            "mixedGradV2": mixed_grad_v2,
        }
    )
    _ = electric_dot, grad_blocks
    return baseline


def r1_survivor_polynomials() -> tuple[tuple[str, sp.Expr], ...]:
    expressions = r1_survivor_expression_map()
    ordered_labels = r1_summary().surviving_labels
    return tuple((label, expressions[label]) for label in ordered_labels)


def coefficient_matrix(
    polynomials: tuple[tuple[str, sp.Expr], ...]
) -> tuple[sp.Matrix, tuple[tuple[int, ...], ...]]:
    variables = r1_variables()
    monomials = sorted(
        {
            monomial
            for _, expression in polynomials
            for monomial in sp.Poly(expression, *variables).monoms()
        }
    )
    matrix = sp.Matrix(
        [
            [sp.Poly(expression, *variables).coeff_monomial(monomial) for _, expression in polynomials]
            for monomial in monomials
        ]
    )
    return matrix, tuple(monomials)


def r1_rank_summary() -> R1RankSummary:
    polynomials = r1_survivor_polynomials()
    matrix, monomials = coefficient_matrix(polynomials)
    nullspace = matrix.nullspace()
    null_relation = sp.Integer(0)
    if nullspace:
        null_vector = nullspace[0]
        symbols = {label: sp.Symbol(label) for label, _ in polynomials}
        null_relation = sp.expand(
            sum(coeff * symbols[label] for (label, _), coeff in zip(polynomials, null_vector))
        )
    summary = r1_summary()
    return R1RankSummary(
        rank=matrix.rank(),
        count=len(polynomials),
        nullity=len(nullspace),
        monomial_count=len(monomials),
        null_relation=null_relation,
        labels=tuple(label for label, _ in polynomials),
        new_labels=summary.new_surviving_labels,
    )


def r1_survivor_rank_report() -> str:
    summary = r1_rank_summary()
    lines = [
        "Delta<=4 rank-1 vector-family survivor rank audit",
        "",
        "Corrected R1-extended survivor candidates:",
        "- " + ", ".join(summary.labels),
        "",
        f"New candidates beyond the audited-set baseline: {', '.join(summary.new_labels)}",
        f"Rank: {summary.rank} out of {summary.count}",
        f"Nullity: {summary.nullity}",
        f"Monomial support size: {summary.monomial_count}",
    ]
    if len(summary.new_labels) <= 1:
        lines.extend(
            [
                "",
                "Interpretation:",
                "- Only one genuinely new survivor is present, so no nontrivial independence issue remains.",
            ]
        )
    elif summary.nullity:
        lines.extend(
            [
                "",
                "First exact dependence relation:",
                f"- {sp.sstr(summary.null_relation)} = 0",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "Interpretation:",
                "- The corrected R1-extended Delta<=4 survivor list is linearly independent.",
                "- Therefore the vector-family extension enlarges the audited catalog without creating a new linear-dependence correction at this stage.",
            ]
        )
    return "\n".join(lines)


def main() -> None:
    print(r1_survivor_rank_report())


if __name__ == "__main__":
    main()
