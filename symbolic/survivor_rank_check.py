from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations

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


def stf3_tensor(prefix: str) -> tuple[list[list[list[sp.Expr]]], tuple[sp.Symbol, ...]]:
    """A totally symmetric, trace-free rank-3 tensor in 3D (an STF-3 octupole).

    grad E = \\partial_k\\partial_i\\partial_j\\Phi is symmetric in ALL three indices
    by equality of mixed partials (Schwarz) and, in the external vacuum
    (\\nabla^2\\Phi = 0, the same condition that makes E traceless), trace-free on
    every pair.  It therefore carries 7 independent components, not the 15 of a
    generic (STF-2)\\otimes vector object.
    """
    xxx, xxy, xxz, xyy, xyz, yyy, yyz = sp.symbols(
        f"{prefix}_xxx {prefix}_xxy {prefix}_xxz {prefix}_xyy {prefix}_xyz {prefix}_yyy {prefix}_yyz"
    )
    # trace-free conditions sum_a G_aaj = 0 fix the remaining three components
    xzz = -(xxx + xyy)
    yzz = -(xxy + yyy)
    zzz = -(xxz + yyz)
    seed = {
        (0, 0, 0): xxx, (0, 0, 1): xxy, (0, 0, 2): xxz,
        (0, 1, 1): xyy, (0, 1, 2): xyz, (0, 2, 2): xzz,
        (1, 1, 1): yyy, (1, 1, 2): yyz, (1, 2, 2): yzz,
        (2, 2, 2): zzz,
    }
    comp: dict[tuple[int, int, int], sp.Expr] = {}
    for key, val in seed.items():
        for perm in set(permutations(key)):
            comp[perm] = val
    tensor = [[[comp[(a, b, c)] for c in range(3)] for b in range(3)] for a in range(3)]
    return tensor, (xxx, xxy, xxz, xyy, xyz, yyy, yyz)


def survivor_polynomials() -> tuple[tuple[str, sp.Expr], ...]:
    electric, electric_vars = stf_matrix("E")
    electric_dot, dot_vars = stf_matrix("DtE")
    grad, grad_vars = stf3_tensor("GradE")

    e2 = sp.expand(sum(electric[i, j] ** 2 for i in range(3) for j in range(3)))
    e3 = sp.expand((electric**3).trace())
    e2_sq = sp.expand(e2**2)
    dot_e2 = sp.expand(sum(electric_dot[i, j] ** 2 for i in range(3) for j in range(3)))
    grad_e2 = sp.expand(
        sum(grad[k][i][j] ** 2 for k in range(3) for i in range(3) for j in range(3))
    )

    _ = electric_vars, dot_vars, grad_vars
    # divE2 and mixedGradE2 are no longer independent survivors: with grad E an
    # STF-3 octupole, divE2 == 0 (trace-free) and mixedGradE2 == gradE2 (total
    # symmetry).  The corrected Delta<=4 electric survivor set is the following five.
    return (
        ("E2", e2),
        ("E3", e3),
        ("E2^2", e2_sq),
        ("dotE2", dot_e2),
        ("gradE2", grad_e2),
    )


def survivor_variables() -> tuple[sp.Symbol, ...]:
    electric, electric_vars = stf_matrix("E")
    electric_dot, dot_vars = stf_matrix("DtE")
    grad, grad_vars = stf3_tensor("GradE")
    _ = electric, electric_dot, grad
    return electric_vars + dot_vars + grad_vars


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
        "- E2, E3, E2^2, dotE2, gradE2",
        "",
        f"Exact polynomial coefficient rank: {summary.total_rank}",
        f"E-sector rank: {summary.e_sector_rank}",
        f"DtE-sector rank: {summary.dt_sector_rank}",
        f"Gradient-sector rank: {summary.gradient_sector_rank}",
        f"Monomial support size: {summary.monomial_count}",
        "",
        "Interpretation:",
        "- The five survivors are linearly independent as operators over constant coefficients.",
        "- The gradient sector is one-dimensional (gradE2): with grad E an STF-3 octupole,",
        "  divE2 vanishes (trace-free) and mixedGradE2 coincides with gradE2 (total symmetry).",
        "- This is a basis-independence statement modulo the allowed reductions, not algebraic functional independence.",
        "- In particular, E2^2 is a separate weight-4 operator even though it is the square of the weight-2 invariant E2.",
    ]
    return "\n".join(lines)


def main() -> None:
    print(survivor_rank_report())


if __name__ == "__main__":
    main()
