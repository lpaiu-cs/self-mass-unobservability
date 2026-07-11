from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import sympy as sp

from r3_sector_delta4 import r3_summary, r3_surviving_classes
from survivor_rank_check import stf_matrix, survivor_polynomials, survivor_variables


@dataclass(frozen=True)
class R3RankSummary:
    rank: int
    count: int
    nullity: int
    monomial_count: int
    null_relation: sp.Expr
    labels: tuple[str, ...]
    new_labels: tuple[str, ...]


def _stf_rank3_array(prefix: str) -> tuple[sp.MutableDenseNDimArray, tuple[sp.Symbol, ...]]:
    c000, c001, c002, c011, c012, c111, c112 = sp.symbols(
        f"{prefix}_000 {prefix}_001 {prefix}_002 {prefix}_011 {prefix}_012 {prefix}_111 {prefix}_112"
    )
    components = {
        (0, 0, 0): c000,
        (0, 0, 1): c001,
        (0, 0, 2): c002,
        (0, 1, 1): c011,
        (0, 1, 2): c012,
        (0, 2, 2): -c000 - c011,
        (1, 1, 1): c111,
        (1, 1, 2): c112,
        (1, 2, 2): -c001 - c111,
        (2, 2, 2): -c002 - c112,
    }
    array = sp.MutableDenseNDimArray.zeros(3, 3, 3)
    for indices, value in components.items():
        for perm in set(
            (
                indices,
                (indices[0], indices[2], indices[1]),
                (indices[1], indices[0], indices[2]),
                (indices[1], indices[2], indices[0]),
                (indices[2], indices[0], indices[1]),
                (indices[2], indices[1], indices[0]),
            )
        ):
            array[perm] = value
    return array, (c000, c001, c002, c011, c012, c111, c112)


def _grad_stf_rank3_array(prefix: str) -> tuple[sp.MutableDenseNDimArray, tuple[sp.Symbol, ...]]:
    array = sp.MutableDenseNDimArray.zeros(3, 3, 3, 3)
    variables: list[sp.Symbol] = []
    for grad_index in range(3):
        block, block_vars = _stf_rank3_array(f"{prefix}{grad_index}")
        variables.extend(block_vars)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    array[grad_index, i, j, k] = block[i, j, k]
    return array, tuple(variables)


def _tensor_map() -> tuple[dict[str, object], tuple[sp.Symbol, ...]]:
    electric, electric_vars = stf_matrix("E")
    electric_dot, electric_dot_vars = stf_matrix("DtE")
    tidal, tidal_vars = _stf_rank3_array("T")
    tidal_dot, tidal_dot_vars = _stf_rank3_array("DtT")
    tidal_grad, tidal_grad_vars = _grad_stf_rank3_array("GradT")
    tensor_map: dict[str, object] = {
        "E": electric,
        "DtE": electric_dot,
        "T": tidal,
        "DtT": tidal_dot,
        "GradT": tidal_grad,
    }
    variables = survivor_variables() + tidal_vars + tidal_dot_vars + tidal_grad_vars
    _ = electric_vars, electric_dot_vars
    return tensor_map, variables


def _component(tensor: object, indices: tuple[int, ...]) -> sp.Expr:
    if isinstance(tensor, sp.MatrixBase):
        return tensor[indices[0], indices[1]]
    if isinstance(tensor, sp.NDimArray):
        return tensor[indices]
    raise TypeError(f"Unsupported tensor container {type(tensor)}")


def _expression_from_representative(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> sp.Expr:
    tensor_map, _ = _tensor_map()
    edge_ids = {edge: edge_index for edge_index, edge in enumerate(representative)}
    slot_to_edge: dict[tuple[int, int], int] = {}
    for edge, edge_index in edge_ids.items():
        slot_to_edge[edge[0]] = edge_index
        slot_to_edge[edge[1]] = edge_index

    expression = sp.Integer(0)
    for assignment in product(range(3), repeat=len(representative)):
        term = sp.Integer(1)
        for instance_index, block_name in enumerate(signature):
            rank = 4 if block_name == "GradT" else (3 if "T" in block_name else 2)
            indices = tuple(assignment[slot_to_edge[(instance_index, slot)]] for slot in range(rank))
            term *= _component(tensor_map[block_name], indices)
        expression += term
    return sp.expand(expression)


def r3_survivor_expression_map() -> dict[str, sp.Expr]:
    expressions = dict(survivor_polynomials())
    for item in r3_surviving_classes():
        if item.label in expressions:
            continue
        expressions[item.label] = _expression_from_representative(item.signature, item.representative)
    return expressions


def r3_survivor_polynomials() -> tuple[tuple[str, sp.Expr], ...]:
    expressions = r3_survivor_expression_map()
    ordered_labels = r3_summary().surviving_labels
    return tuple((label, expressions[label]) for label in ordered_labels)


def coefficient_matrix(
    polynomials: tuple[tuple[str, sp.Expr], ...]
) -> tuple[sp.Matrix, tuple[tuple[int, ...], ...]]:
    _, variables = _tensor_map()
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


def r3_rank_summary() -> R3RankSummary:
    polynomials = r3_survivor_polynomials()
    matrix, monomials = coefficient_matrix(polynomials)
    nullspace = matrix.nullspace()
    null_relation = sp.Integer(0)
    if nullspace:
        null_vector = nullspace[0]
        symbols = {label: sp.Symbol(label) for label, _ in polynomials}
        null_relation = sp.expand(
            sum(coeff * symbols[label] for (label, _), coeff in zip(polynomials, null_vector))
        )
    summary = r3_summary()
    return R3RankSummary(
        rank=matrix.rank(),
        count=len(polynomials),
        nullity=len(nullspace),
        monomial_count=len(monomials),
        null_relation=null_relation,
        labels=tuple(label for label, _ in polynomials),
        new_labels=summary.new_surviving_labels,
    )


def r3_survivor_rank_report() -> str:
    summary = r3_rank_summary()
    lines = [
        "Delta<=4 rank-3 family survivor rank audit",
        "",
        "Corrected R3-extended survivor candidates:",
        "- " + ", ".join(summary.labels),
        "",
        f"New candidates beyond the enlarged audited-set baseline: {', '.join(summary.new_labels)}",
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
                "- The corrected R3-extended Delta<=4 survivor list is linearly independent.",
                "- Therefore the genuine rank-3 family enlarges the audited catalog without requiring a linear-dependence correction at this stage.",
            ]
        )
    return "\n".join(lines)


def main() -> None:
    print(r3_survivor_rank_report())


if __name__ == "__main__":
    main()
