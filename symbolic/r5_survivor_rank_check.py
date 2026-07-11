from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import random

import sympy as sp

from survivor_rank_check import rank_summary
from r5_sector_delta4 import r5_summary


SAMPLE_COUNT = 8
RNG_SEED = 525252

PAIR_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


@dataclass(frozen=True)
class R5RankSummary:
    rank: int
    count: int
    nullity: int
    baseline_rank: int
    new_rank: int
    new_count: int
    sample_count: int
    first_null_relation: sp.Expr
    additional_null_relations: tuple[sp.Expr, ...]
    labels: tuple[str, ...]
    new_labels: tuple[str, ...]


def _count_patterns(rank: int) -> tuple[tuple[int, int, int], ...]:
    patterns: list[tuple[int, int, int]] = []
    for count0 in range(rank + 1):
        for count1 in range(rank - count0 + 1):
            patterns.append((count0, count1, rank - count0 - count1))
    return tuple(patterns)


def _counts(indices: tuple[int, ...]) -> tuple[int, int, int]:
    return (indices.count(0), indices.count(1), indices.count(2))


def _rank5_entry_order() -> tuple[tuple[int, int, int, int, int], ...]:
    return tuple(
        (i, j, k, l, m)
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for l in range(3)
        for m in range(3)
    )


@lru_cache(maxsize=1)
def _stf_rank5_basis() -> tuple[sp.MutableDenseNDimArray, ...]:
    patterns5 = _count_patterns(5)
    pattern_index5 = {pattern: idx for idx, pattern in enumerate(patterns5)}
    patterns3 = _count_patterns(3)
    constraint_rows: list[list[int]] = []
    for pattern3 in patterns3:
        ijk = tuple([0] * pattern3[0] + [1] * pattern3[1] + [2] * pattern3[2])
        row = [0] * len(patterns5)
        for trace_index in range(3):
            row[pattern_index5[_counts((trace_index, trace_index) + ijk)]] += 1
        constraint_rows.append(row)
    nullspace = sp.Matrix(constraint_rows).nullspace()
    entries5 = _rank5_entry_order()
    basis: list[sp.MutableDenseNDimArray] = []
    for vector in nullspace:
        array = sp.MutableDenseNDimArray.zeros(3, 3, 3, 3, 3)
        for indices in entries5:
            array[indices] = vector[pattern_index5[_counts(indices)]]
        basis.append(array)
    return tuple(basis)


def _random_stf2(rng: random.Random) -> sp.Matrix:
    a, b, c, d, e = [sp.Integer(rng.randint(-3, 3)) for _ in range(5)]
    return sp.Matrix([[a, b, c], [b, d, e], [c, e, -a - d]])


@lru_cache(maxsize=1)
def _rank5_index() -> dict[tuple[int, int, int, int, int], int]:
    return {indices: idx for idx, indices in enumerate(_rank5_entry_order())}


def _random_stf5(rng: random.Random) -> list[int]:
    basis = _stf_rank5_basis()
    entries5 = _rank5_entry_order()
    coefficients = [rng.randint(-2, 2) for _ in range(len(basis))]
    while all(value == 0 for value in coefficients):
        coefficients = [rng.randint(-2, 2) for _ in range(len(basis))]
    values = [0] * len(entries5)
    for coeff, basis_array in zip(coefficients, basis):
        for index, indices in enumerate(entries5):
            values[index] += coeff * int(basis_array[indices])
    return values


def _random_grad_stf5(rng: random.Random) -> list[list[int]]:
    return [_random_stf5(rng) for _ in range(3)]


def _multigraph_representative(
    degrees: tuple[int, ...],
    multiplicities: tuple[int, int, int, int, int, int],
) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    slot_cursors = [0] * len(degrees)
    representative: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for (left, right), multiplicity in zip(PAIR_EDGES, multiplicities):
        for _ in range(multiplicity):
            representative.append(
                ((left, slot_cursors[left]), (right, slot_cursors[right]))
            )
            slot_cursors[left] += 1
            slot_cursors[right] += 1
    if tuple(slot_cursors) != degrees:
        raise ValueError(
            f"Representative construction drifted: got {tuple(slot_cursors)} expected {degrees}."
        )
    return tuple(representative)


def _r5_new_specs() -> dict[str, tuple[tuple[str, ...], tuple[tuple[tuple[int, int], tuple[int, int]], ...]]]:
    return {
        "U2": (
            ("U", "U"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
            ),
        ),
        "EUU": (
            ("E", "U", "U"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
                ((1, 3), (2, 3)),
                ((1, 4), (2, 4)),
            ),
        ),
        "dotU2": (
            ("DtU", "DtU"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
            ),
        ),
        "EUDtU": (
            ("DtU", "E", "U"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((0, 2), (2, 1)),
                ((0, 3), (2, 2)),
                ((0, 4), (2, 3)),
                ((1, 1), (2, 4)),
            ),
        ),
        "E2U2": (
            ("E", "E", "U", "U"),
            _multigraph_representative((2, 2, 5, 5), (2, 0, 0, 0, 0, 5)),
        ),
        "E2U2_mixed_1": (
            ("E", "E", "U", "U"),
            _multigraph_representative((2, 2, 5, 5), (1, 0, 1, 1, 0, 4)),
        ),
        "E2U2_mixed_2": (
            ("E", "E", "U", "U"),
            _multigraph_representative((2, 2, 5, 5), (0, 0, 2, 2, 0, 3)),
        ),
        "E2U2_mixed_3": (
            ("E", "E", "U", "U"),
            _multigraph_representative((2, 2, 5, 5), (0, 1, 1, 1, 1, 3)),
        ),
        "divU2": (
            ("GradU", "GradU"),
            (
                ((0, 0), (0, 1)),
                ((0, 2), (1, 1)),
                ((0, 3), (1, 2)),
                ((0, 4), (1, 3)),
                ((0, 5), (1, 4)),
                ((1, 0), (1, 5)),
            ),
        ),
        "gradU2": (
            ("GradU", "GradU"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
                ((0, 5), (1, 5)),
            ),
        ),
        "mixedGradU2": (
            ("GradU", "GradU"),
            (
                ((0, 0), (1, 1)),
                ((0, 1), (1, 0)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
                ((0, 5), (1, 5)),
            ),
        ),
        "U2^2": (
            ("U", "U", "U", "U"),
            _multigraph_representative((5, 5, 5, 5), (0, 0, 5, 5, 0, 0)),
        ),
        "U4_chain": (
            ("U", "U", "U", "U"),
            _multigraph_representative((5, 5, 5, 5), (0, 1, 4, 4, 1, 0)),
        ),
        "U4_bridge": (
            ("U", "U", "U", "U"),
            _multigraph_representative((5, 5, 5, 5), (0, 2, 3, 3, 2, 0)),
        ),
        "U4_tetra": (
            ("U", "U", "U", "U"),
            _multigraph_representative((5, 5, 5, 5), (1, 1, 3, 3, 1, 1)),
        ),
        "U4_balanced": (
            ("U", "U", "U", "U"),
            _multigraph_representative((5, 5, 5, 5), (1, 2, 2, 2, 2, 1)),
        ),
    }


def _numeric_spec_meta() -> dict[str, tuple[tuple[str, ...], dict[tuple[int, int], int], int]]:
    out: dict[str, tuple[tuple[str, ...], dict[tuple[int, int], int], int]] = {}
    for label, (signature, representative) in _r5_new_specs().items():
        slot_to_edge: dict[tuple[int, int], int] = {}
        for edge_index, edge in enumerate(representative):
            slot_to_edge[edge[0]] = edge_index
            slot_to_edge[edge[1]] = edge_index
        out[label] = (signature, slot_to_edge, len(representative))
    return out


def _u_component(values: list[int], indices: tuple[int, int, int, int, int]) -> int:
    return values[_rank5_index()[indices]]


def _eval_new_label(
    sample: dict[str, object],
    label: str,
    metadata: dict[str, tuple[tuple[str, ...], dict[tuple[int, int], int], int]],
) -> int:
    signature, slot_to_edge, edge_count = metadata[label]
    total = 0
    for assignment in product(range(3), repeat=edge_count):
        term = 1
        for instance_index, block_name in enumerate(signature):
            if block_name == "E":
                indices = tuple(assignment[slot_to_edge[(instance_index, slot)]] for slot in range(2))
                matrix = sample["E"]
                assert isinstance(matrix, sp.MatrixBase)
                term *= int(matrix[indices[0], indices[1]])
            elif block_name in {"U", "DtU"}:
                indices = tuple(
                    assignment[slot_to_edge[(instance_index, slot)]] for slot in range(5)
                )
                values = sample[block_name]
                assert isinstance(values, list)
                term *= _u_component(values, indices)
            elif block_name == "GradU":
                indices = tuple(
                    assignment[slot_to_edge[(instance_index, slot)]] for slot in range(6)
                )
                values = sample["GradU"]
                assert isinstance(values, list)
                term *= values[indices[0]][_rank5_index()[indices[1:]]]
            else:
                raise ValueError(f"Unsupported block name {block_name}")
        total += term
    return total


def _normalize_nullspace_vector(
    labels: tuple[str, ...],
    vector: sp.Matrix,
) -> tuple[tuple[int, ...], sp.Expr]:
    denominators = [
        value.q for value in vector if isinstance(value, sp.Rational) and value.q != 1
    ]
    if not denominators:
        scale = 1
    elif len(denominators) == 1:
        scale = denominators[0]
    else:
        scale = sp.ilcm(*denominators)
    integers = [sp.Integer(scale) * value for value in vector]
    gcd = sp.igcd(*[abs(int(value)) for value in integers if value != 0]) if any(
        value != 0 for value in integers
    ) else 1
    integers = [sp.Integer(value // gcd) for value in integers]
    first_nonzero = next((value for value in integers if value != 0), sp.Integer(1))
    if first_nonzero < 0:
        integers = [-value for value in integers]
    relation = sp.expand(
        sum(coefficient * sp.Symbol(label) for label, coefficient in zip(labels, integers))
    )
    return tuple(int(value) for value in integers), relation


def _sampled_null_relations(
    labels: tuple[str, ...],
    rows: list[list[int]],
) -> tuple[tuple[tuple[int, ...], sp.Expr], ...]:
    matrix = sp.Matrix(rows)
    relations: list[sp.Expr] = []
    for vector in matrix.nullspace():
        relations.append(_normalize_nullspace_vector(labels, vector))
    return tuple(relations)


def _sample_rows_for_labels(
    labels: tuple[str, ...],
    sample_count: int,
    seed: int,
) -> list[list[int]]:
    metadata = _numeric_spec_meta()
    rows: list[list[int]] = []
    rng = random.Random(seed)
    for _ in range(sample_count):
        sample = {
            "E": _random_stf2(rng),
            "U": _random_stf5(rng),
            "DtU": _random_stf5(rng),
            "GradU": _random_grad_stf5(rng),
        }
        rows.append([_eval_new_label(sample, label, metadata) for label in labels])
    return rows


def _sector_rank(labels: tuple[str, ...], sample_count: int, seed: int) -> int:
    return sp.Matrix(_sample_rows_for_labels(labels, sample_count=sample_count, seed=seed)).rank()


def _verify_relations_on_rows(
    relations: tuple[tuple[tuple[int, ...], sp.Expr], ...],
    labels: tuple[str, ...],
    rows: list[list[int]],
) -> tuple[sp.Expr, ...]:
    verified: list[sp.Expr] = []
    for coefficients, relation in relations:
        for row in rows:
            if sum(coefficient * value for coefficient, value in zip(coefficients, row)) != 0:
                raise ValueError(
                    f"Rank-5 sample-based null relation failed revalidation for labels {labels}."
                )
        verified.append(relation)
    return tuple(verified)


def _verified_null_relations() -> tuple[sp.Expr, ...]:
    mixed_labels = ("E2U2", "E2U2_mixed_1", "E2U2_mixed_2", "E2U2_mixed_3")
    mixed_relations = _sampled_null_relations(
        mixed_labels,
        _sample_rows_for_labels(mixed_labels, sample_count=6, seed=RNG_SEED + 11),
    )
    verified_mixed = _verify_relations_on_rows(
        mixed_relations,
        mixed_labels,
        _sample_rows_for_labels(mixed_labels, sample_count=4, seed=RNG_SEED + 17),
    )

    quartic_labels = ("U2^2", "U4_chain", "U4_bridge", "U4_tetra", "U4_balanced")
    quartic_relations = _sampled_null_relations(
        quartic_labels,
        _sample_rows_for_labels(quartic_labels, sample_count=8, seed=RNG_SEED + 29),
    )
    verified_quartic = _verify_relations_on_rows(
        quartic_relations,
        quartic_labels,
        _sample_rows_for_labels(quartic_labels, sample_count=4, seed=RNG_SEED + 41),
    )
    return verified_mixed + verified_quartic


def r5_rank_summary() -> R5RankSummary:
    summary = r5_summary()
    baseline = rank_summary()
    null_relations = _verified_null_relations()
    witness_labels = ("U2", "EUU", "dotU2", "EUDtU")
    mixed_quadratic_labels = ("E2U2", "E2U2_mixed_1", "E2U2_mixed_2", "E2U2_mixed_3")
    gradient_labels = ("divU2", "gradU2", "mixedGradU2")
    quartic_labels = ("U2^2", "U4_chain", "U4_bridge", "U4_tetra", "U4_balanced")
    sector_ranks = (
        _sector_rank(witness_labels[:1], sample_count=4, seed=RNG_SEED + 51),
        _sector_rank(witness_labels[1:2], sample_count=4, seed=RNG_SEED + 53),
        _sector_rank(witness_labels[2:3], sample_count=4, seed=RNG_SEED + 55),
        _sector_rank(witness_labels[3:], sample_count=4, seed=RNG_SEED + 57),
        _sector_rank(mixed_quadratic_labels, sample_count=6, seed=RNG_SEED + 59),
        _sector_rank(gradient_labels, sample_count=5, seed=RNG_SEED + 61),
        _sector_rank(quartic_labels, sample_count=8, seed=RNG_SEED + 63),
    )
    new_rank_lower_bound = sum(sector_ranks)
    new_count = len(summary.new_surviving_labels)
    new_nullity = len(null_relations)
    if new_rank_lower_bound != new_count - new_nullity:
        raise ValueError(
            "Rank-5 new-sector lower bound does not match the verified null-relation upper bound."
        )
    total_rank = baseline.total_rank + new_rank_lower_bound
    total_count = len(summary.surviving_labels)
    total_nullity = total_count - total_rank
    return R5RankSummary(
        rank=total_rank,
        count=total_count,
        nullity=total_nullity,
        baseline_rank=baseline.total_rank,
        new_rank=new_rank_lower_bound,
        new_count=new_count,
        sample_count=SAMPLE_COUNT,
        first_null_relation=null_relations[0],
        additional_null_relations=null_relations[1:],
        labels=summary.surviving_labels,
        new_labels=summary.new_surviving_labels,
    )


def r5_survivor_rank_report() -> str:
    summary = r5_rank_summary()
    lines = [
        "Delta<=4 rank-5 family survivor rank audit",
        "",
        "Raw R5-extended survivor candidates:",
        "- " + ", ".join(summary.labels),
        "",
        "New candidates beyond the enlarged audited-set baseline:",
        "- " + ", ".join(summary.new_labels),
        "",
        f"Baseline electric rank: {summary.baseline_rank}",
        f"New rank-5-sector rank: {summary.new_rank} out of {summary.new_count}",
        f"Total rank: {summary.rank} out of {summary.count}",
        f"Nullity: {summary.nullity}",
        f"Deterministic evaluation sample count: {summary.sample_count}",
        "",
        "First exact dependence relation:",
        f"- {sp.sstr(summary.first_null_relation)} = 0",
    ]
    if summary.additional_null_relations:
        lines.extend(
            [
                "",
                "Additional exact dependence relations:",
                *[f"- {sp.sstr(relation)} = 0" for relation in summary.additional_null_relations],
            ]
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "- The raw R5-extended survivor list is not a corrected basis statement.",
            "- The rank result uses two ingredients: sample-extracted null relations revalidated on an independent exact-integer batch, and a deterministic component-evaluation lower bound.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(r5_survivor_rank_report())


if __name__ == "__main__":
    main()
