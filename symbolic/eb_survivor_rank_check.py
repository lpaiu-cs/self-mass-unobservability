from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class EBRankSummary:
    raw_rank: int
    raw_count: int
    corrected_rank: int
    corrected_count: int
    nullity: int
    monomial_count: int
    null_relation: sp.Expr
    raw_labels: tuple[str, ...]
    corrected_labels: tuple[str, ...]


def stf_matrix(prefix: str) -> tuple[sp.Matrix, tuple[sp.Symbol, ...]]:
    xx, xy, xz, yy, yz = sp.symbols(f"{prefix}_xx {prefix}_xy {prefix}_xz {prefix}_yy {prefix}_yz")
    matrix = sp.Matrix(
        [
            [xx, xy, xz],
            [xy, yy, yz],
            [xz, yz, -xx - yy],
        ]
    )
    return matrix, (xx, xy, xz, yy, yz)


def frobenius_pair(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand(sum(left[i, j] * right[i, j] for i in range(3) for j in range(3)))


def square_norm(matrix: sp.Matrix) -> sp.Expr:
    return sp.expand(sum(matrix[i, j] ** 2 for i in range(3) for j in range(3)))


def eb_variables() -> tuple[sp.Symbol, ...]:
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
    return electric_vars + magnetic_vars + electric_dot_vars + magnetic_dot_vars + tuple(grad_vars)


def raw_eb_survivor_polynomials() -> tuple[tuple[str, sp.Expr], ...]:
    electric, _ = stf_matrix("E")
    magnetic, _ = stf_matrix("B")
    electric_dot, _ = stf_matrix("DtE")
    magnetic_dot, _ = stf_matrix("DtB")

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
    ebeb = sp.expand((electric * magnetic * electric * magnetic).trace())
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

    return (
        ("E2", e2),
        ("B2", b2),
        ("E3", e3),
        ("EB2", eb2),
        ("E2^2", e2_sq),
        ("B2^2", b2_sq),
        ("dotE2", dot_e2),
        ("dotB2", dot_b2),
        ("EBDtB", ebdtb),
        ("E2B2", e2b2),
        ("EB_sq", eb_sq),
        ("TrE2B2", tr_e2b2),
        ("EBEB", ebeb),
        ("gradE2", grad_e2),
        ("divE2", div_e2),
        ("mixedGradE2", mixed_grad_e2),
        ("gradB2", grad_b2),
        ("divB2", div_b2),
        ("mixedGradB2", mixed_grad_b2),
    )


def corrected_eb_survivor_polynomials() -> tuple[tuple[str, sp.Expr], ...]:
    return tuple(item for item in raw_eb_survivor_polynomials() if item[0] != "EBEB")


def coefficient_matrix(
    polynomials: tuple[tuple[str, sp.Expr], ...]
) -> tuple[sp.Matrix, tuple[tuple[int, ...], ...]]:
    variables = eb_variables()
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


def eb_rank_summary() -> EBRankSummary:
    raw = raw_eb_survivor_polynomials()
    corrected = corrected_eb_survivor_polynomials()
    raw_matrix, monomials = coefficient_matrix(raw)
    corrected_matrix, _ = coefficient_matrix(corrected)
    nullspace = raw_matrix.nullspace()
    labels = [label for label, _ in raw]
    null_relation = sp.Integer(0)
    if nullspace:
        null_vector = nullspace[0]
        symbols = {label: sp.Symbol(label) for label in labels}
        null_relation = sp.expand(
            sum(coeff * symbols[label] for label, coeff in zip(labels, null_vector))
        )
    return EBRankSummary(
        raw_rank=raw_matrix.rank(),
        raw_count=len(raw),
        corrected_rank=corrected_matrix.rank(),
        corrected_count=len(corrected),
        nullity=len(nullspace),
        monomial_count=len(monomials),
        null_relation=null_relation,
        raw_labels=tuple(label for label, _ in raw),
        corrected_labels=tuple(label for label, _ in corrected),
    )


def eb_survivor_rank_report() -> str:
    summary = eb_rank_summary()
    lines = [
        "Delta<=4 E/B survivor rank audit",
        "",
        "Raw E/B survivor candidates:",
        "- " + ", ".join(summary.raw_labels),
        "",
        f"Raw rank: {summary.raw_rank} out of {summary.raw_count}",
        f"Corrected rank: {summary.corrected_rank} out of {summary.corrected_count}",
        f"Nullity of the raw set: {summary.nullity}",
        f"Monomial support size: {summary.monomial_count}",
        "",
        "First exact dependence relation:",
        f"- {sp.sstr(summary.null_relation)} = 0",
        "",
        "Corrected E/B basis:",
        "- " + ", ".join(summary.corrected_labels),
        "",
        "Interpretation:",
        "- The raw 19-element E/B survivor list is not linearly independent.",
        "- The first explicit dependence sits in the mixed quartic algebraic sector.",
        "- After removing EBEB via the explicit quartic STF identity, the corrected 18-element E/B basis is linearly independent.",
    ]
    return "\n".join(lines)


def main() -> None:
    print(eb_survivor_rank_report())


if __name__ == "__main__":
    main()
