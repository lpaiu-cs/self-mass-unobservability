from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from eb_survivor_rank_check import (
    frobenius_pair,
    square_norm,
    stf_matrix,
)
from es_sector_delta4 import es_summary


@dataclass(frozen=True)
class ESRankSummary:
    rank: int
    count: int
    nullity: int
    monomial_count: int
    null_relation: sp.Expr
    labels: tuple[str, ...]


def es_variables() -> tuple[sp.Symbol, ...]:
    _, electric_vars = stf_matrix("E")
    _, magnetic_vars = stf_matrix("B")
    _, electric_dot_vars = stf_matrix("DtE")
    _, magnetic_dot_vars = stf_matrix("DtB")
    grad_vars: list[sp.Symbol] = []
    for grad_index in range(3):
        _, block_vars = stf_matrix(f"GradE{grad_index}")
        grad_vars.extend(block_vars)
    for grad_index in range(3):
        _, block_vars = stf_matrix(f"GradB{grad_index}")
        grad_vars.extend(block_vars)
    scalar_vars = sp.symbols("S DtS")
    grad_s_vars = sp.symbols("GradS0 GradS1 GradS2")
    return (
        electric_vars
        + magnetic_vars
        + electric_dot_vars
        + magnetic_dot_vars
        + tuple(grad_vars)
        + scalar_vars
        + grad_s_vars
    )


def es_survivor_expression_map() -> dict[str, sp.Expr]:
    electric, _ = stf_matrix("E")
    magnetic, _ = stf_matrix("B")
    electric_dot, _ = stf_matrix("DtE")
    magnetic_dot, _ = stf_matrix("DtB")
    scalar, scalar_dot = sp.symbols("S DtS")
    grad_s = sp.symbols("GradS0 GradS1 GradS2")

    electric_grad_blocks: list[sp.Matrix] = []
    magnetic_grad_blocks: list[sp.Matrix] = []
    for grad_index in range(3):
        block, _ = stf_matrix(f"GradE{grad_index}")
        electric_grad_blocks.append(block)
    for grad_index in range(3):
        block, _ = stf_matrix(f"GradB{grad_index}")
        magnetic_grad_blocks.append(block)

    e2 = square_norm(electric)
    b2 = square_norm(magnetic)
    e3 = sp.expand((electric**3).trace())
    eb2 = sp.expand((electric * (magnetic**2)).trace())
    e2_sq = sp.expand(e2**2)
    b2_sq = sp.expand(b2**2)
    dot_e2 = square_norm(electric_dot)
    dot_b2 = square_norm(magnetic_dot)
    ebdtb = sp.expand((electric * magnetic * magnetic_dot).trace())
    e2b2 = sp.expand(e2 * b2)
    eb_sq = sp.expand(frobenius_pair(electric, magnetic) ** 2)
    tr_e2b2 = sp.expand((electric**2 * magnetic**2).trace())
    grad_e2 = sp.expand(
        sum(
            electric_grad_blocks[k][i, j] ** 2
            for k in range(3)
            for i in range(3)
            for j in range(3)
        )
    )
    div_e = [sp.expand(sum(electric_grad_blocks[i][i, j] for i in range(3))) for j in range(3)]
    div_e2 = sp.expand(sum(entry**2 for entry in div_e))
    mixed_grad_e2 = sp.expand(
        sum(
            electric_grad_blocks[k][i, j] * electric_grad_blocks[i][k, j]
            for i in range(3)
            for j in range(3)
            for k in range(3)
        )
    )
    grad_b2 = sp.expand(
        sum(
            magnetic_grad_blocks[k][i, j] ** 2
            for k in range(3)
            for i in range(3)
            for j in range(3)
        )
    )
    div_b = [sp.expand(sum(magnetic_grad_blocks[i][i, j] for i in range(3))) for j in range(3)]
    div_b2 = sp.expand(sum(entry**2 for entry in div_b))
    mixed_grad_b2 = sp.expand(
        sum(
            magnetic_grad_blocks[k][i, j] * magnetic_grad_blocks[i][k, j]
            for i in range(3)
            for j in range(3)
            for k in range(3)
        )
    )
    grad_s2 = sp.expand(sum(component**2 for component in grad_s))
    div_e_grad_s = sp.expand(sum(div_e[index] * grad_s[index] for index in range(3)))

    return {
        "S": scalar,
        "B2": b2,
        "E2": e2,
        "S2": sp.expand(scalar**2),
        "SB2": sp.expand(scalar * b2),
        "SE2": sp.expand(scalar * e2),
        "S3": sp.expand(scalar**3),
        "EB2": eb2,
        "E3": e3,
        "B2^2": b2_sq,
        "DtS_B2": sp.expand(scalar_dot * b2),
        "DtS_E2": sp.expand(scalar_dot * e2),
        "E2B2": e2b2,
        "EB_sq": eb_sq,
        "TrE2B2": tr_e2b2,
        "SE3": sp.expand(scalar * e3),
        "SEB2": sp.expand(scalar * eb2),
        "dotB2": dot_b2,
        "dotE2": dot_e2,
        "dotS2": sp.expand(scalar_dot**2),
        "S2B2": sp.expand((scalar**2) * b2),
        "S2E2": sp.expand((scalar**2) * e2),
        "S4": sp.expand(scalar**4),
        "E2^2": e2_sq,
        "divB2": div_b2,
        "divE2": div_e2,
        "divEGradS": div_e_grad_s,
        "EBDtB": ebdtb,
        "gradB2": grad_b2,
        "gradE2": grad_e2,
        "gradS2": grad_s2,
        "mixedGradB2": mixed_grad_b2,
        "mixedGradE2": mixed_grad_e2,
    }


def es_survivor_polynomials() -> tuple[tuple[str, sp.Expr], ...]:
    expressions = es_survivor_expression_map()
    ordered_labels = es_summary().surviving_labels
    return tuple((label, expressions[label]) for label in ordered_labels)


def coefficient_matrix(
    polynomials: tuple[tuple[str, sp.Expr], ...]
) -> tuple[sp.Matrix, tuple[tuple[int, ...], ...]]:
    variables = es_variables()
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


def es_rank_summary() -> ESRankSummary:
    polynomials = es_survivor_polynomials()
    matrix, monomials = coefficient_matrix(polynomials)
    nullspace = matrix.nullspace()
    null_relation = sp.Integer(0)
    if nullspace:
        null_vector = nullspace[0]
        symbols = {label: sp.Symbol(label) for label, _ in polynomials}
        null_relation = sp.expand(
            sum(coeff * symbols[label] for (label, _), coeff in zip(polynomials, null_vector))
        )
    return ESRankSummary(
        rank=matrix.rank(),
        count=len(polynomials),
        nullity=len(nullspace),
        monomial_count=len(monomials),
        null_relation=null_relation,
        labels=tuple(label for label, _ in polynomials),
    )


def es_survivor_rank_report() -> str:
    summary = es_rank_summary()
    lines = [
        "Delta<=4 E/B+scalar survivor rank audit",
        "",
        "Corrected E/B+scalar survivor candidates:",
        "- " + ", ".join(summary.labels),
        "",
        f"Rank: {summary.rank} out of {summary.count}",
        f"Nullity: {summary.nullity}",
        f"Monomial support size: {summary.monomial_count}",
    ]
    if summary.nullity:
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
                "- No new linear dependence was found in the corrected E/B+scalar Delta<=4 survivor list.",
                "- The scalar-family extension enlarges the basis but does not presently break fixed-order independence.",
            ]
        )
    return "\n".join(lines)


def main() -> None:
    print(es_survivor_rank_report())


if __name__ == "__main__":
    main()
