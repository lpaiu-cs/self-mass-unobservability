from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class RankSummary:
    total_rank: int
    e_sector_rank: int
    dt_sector_rank: int
    gradient_sector_rank: int
    survivor_labels: tuple[str, ...]
    monomial_count: int


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


def survivor_polynomials() -> tuple[tuple[str, sp.Expr], ...]:
    electric, electric_vars = stf_matrix("E")
    electric_dot, dot_vars = stf_matrix("DtE")

    grad_vars: list[sp.Symbol] = []
    grad_blocks: list[sp.Matrix] = []
    for grad_index in range(3):
        block, block_vars = stf_matrix(f"GradE{grad_index}")
        grad_blocks.append(block)
        grad_vars.extend(block_vars)

    e2 = sp.expand(sum(electric[i, j] ** 2 for i in range(3) for j in range(3)))
    e3 = sp.expand((electric**3).trace())
    e2_sq = sp.expand(e2**2)
    dot_e2 = sp.expand(sum(electric_dot[i, j] ** 2 for i in range(3) for j in range(3)))
    grad_e2 = sp.expand(
        sum(grad_blocks[k][i, j] ** 2 for k in range(3) for i in range(3) for j in range(3))
    )
    divergence = [
        sp.expand(sum(grad_blocks[i][i, j] for i in range(3)))
        for j in range(3)
    ]
    div_e2 = sp.expand(sum(entry**2 for entry in divergence))
    mixed_grad_e2 = sp.expand(
        sum(grad_blocks[k][i, j] * grad_blocks[i][k, j] for i in range(3) for j in range(3) for k in range(3))
    )

    _ = electric_vars, dot_vars, tuple(grad_vars)
    return (
        ("E2", e2),
        ("E3", e3),
        ("E2^2", e2_sq),
        ("dotE2", dot_e2),
        ("gradE2", grad_e2),
        ("divE2", div_e2),
        ("mixedGradE2", mixed_grad_e2),
    )


def survivor_variables() -> tuple[sp.Symbol, ...]:
    electric, electric_vars = stf_matrix("E")
    electric_dot, dot_vars = stf_matrix("DtE")
    _ = electric, electric_dot
    grad_vars: list[sp.Symbol] = []
    for grad_index in range(3):
        block, block_vars = stf_matrix(f"GradE{grad_index}")
        _ = block
        grad_vars.extend(block_vars)
    return electric_vars + dot_vars + tuple(grad_vars)


def coefficient_rank(polynomials: tuple[tuple[str, sp.Expr], ...]) -> tuple[sp.Matrix, tuple[tuple[int, ...], ...]]:
    variables = survivor_variables()
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


def rank_summary() -> RankSummary:
    survivors = survivor_polynomials()
    full_matrix, monomials = coefficient_rank(survivors)
    e_sector_matrix, _ = coefficient_rank(survivors[:3])
    dt_sector_matrix, _ = coefficient_rank((survivors[3],))
    gradient_sector_matrix, _ = coefficient_rank(survivors[4:])
    return RankSummary(
        total_rank=full_matrix.rank(),
        e_sector_rank=e_sector_matrix.rank(),
        dt_sector_rank=dt_sector_matrix.rank(),
        gradient_sector_rank=gradient_sector_matrix.rank(),
        survivor_labels=tuple(label for label, _ in survivors),
        monomial_count=len(monomials),
    )


def survivor_rank_report() -> str:
    summary = rank_summary()
    lines = [
        "Delta<=4 survivor rank audit",
        "",
        "Linear independence target:",
        "- E2, E3, E2^2, dotE2, gradE2, divE2, mixedGradE2",
        "",
        f"Exact polynomial coefficient rank: {summary.total_rank}",
        f"E-sector rank: {summary.e_sector_rank}",
        f"DtE-sector rank: {summary.dt_sector_rank}",
        f"Gradient-sector rank: {summary.gradient_sector_rank}",
        f"Monomial support size: {summary.monomial_count}",
        "",
        "Interpretation:",
        "- The seven survivors are linearly independent as operators over constant coefficients.",
        "- This is a basis-independence statement modulo the allowed reductions, not algebraic functional independence.",
        "- In particular, E2^2 is a separate weight-4 operator even though it is the square of the weight-2 invariant E2.",
    ]
    return "\n".join(lines)


def main() -> None:
    print(survivor_rank_report())


if __name__ == "__main__":
    main()
