from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import random

import sympy as sp

from r6_sector_delta4 import enumerate_r6_family_classes, r6_summary
from survivor_rank_check import rank_summary


RNG_SEED = 626262
VERIFICATION_ROWS = 4


@dataclass(frozen=True)
class R6SectorSignatureRank:
    signature: tuple[str, ...]
    labels: tuple[str, ...]
    sample_count: int
    rank: int
    count: int
    nullity: int
    first_null_relation: sp.Expr | None


@dataclass(frozen=True)
class R6RankSummary:
    rank: int
    count: int
    nullity: int
    baseline_rank: int
    new_rank: int
    new_count: int
    signature_ranks: tuple[R6SectorSignatureRank, ...]
    first_null_relation: sp.Expr | None
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


@lru_cache(maxsize=None)
def _entry_order(rank: int) -> tuple[tuple[int, ...], ...]:
    if rank == 2:
        return tuple((i, j) for i in range(3) for j in range(3))
    if rank == 6:
        return tuple(
            (i, j, k, l, m, n)
            for i in range(3)
            for j in range(3)
            for k in range(3)
            for l in range(3)
            for m in range(3)
            for n in range(3)
        )
    raise ValueError(f"Unsupported STF rank {rank}")


@lru_cache(maxsize=None)
def _index_map(rank: int) -> dict[tuple[int, ...], int]:
    return {indices: idx for idx, indices in enumerate(_entry_order(rank))}


@lru_cache(maxsize=None)
def _stf_basis(rank: int) -> tuple[sp.MutableDenseNDimArray, ...]:
    patterns = _count_patterns(rank)
    pattern_index = {pattern: idx for idx, pattern in enumerate(patterns)}
    lower_patterns = _count_patterns(rank - 2)
    constraint_rows: list[list[int]] = []
    for lower_pattern in lower_patterns:
        suffix = tuple([0] * lower_pattern[0] + [1] * lower_pattern[1] + [2] * lower_pattern[2])
        row = [0] * len(patterns)
        for trace_index in range(3):
            row[pattern_index[_counts((trace_index, trace_index) + suffix)]] += 1
        constraint_rows.append(row)

    nullspace = sp.Matrix(constraint_rows).nullspace()
    basis: list[sp.MutableDenseNDimArray] = []
    for vector in nullspace:
        array = sp.MutableDenseNDimArray.zeros(*((3,) * rank))
        for indices in _entry_order(rank):
            array[indices] = vector[pattern_index[_counts(indices)]]
        basis.append(array)
    return tuple(basis)


def _random_stf(rank: int, rng: random.Random) -> list[int]:
    basis = _stf_basis(rank)
    entries = _entry_order(rank)
    coefficients = [rng.randint(-2, 2) for _ in range(len(basis))]
    while all(value == 0 for value in coefficients):
        coefficients = [rng.randint(-2, 2) for _ in range(len(basis))]
    values = [0] * len(entries)
    for coeff, basis_array in zip(coefficients, basis):
        for index, indices in enumerate(entries):
            values[index] += coeff * int(basis_array[indices])
    return values


def _random_grad_stf(rank: int, rng: random.Random) -> list[list[int]]:
    return [_random_stf(rank, rng) for _ in range(3)]


def _random_stf2(rng: random.Random) -> sp.Matrix:
    a, b, c, d, e = [sp.Integer(rng.randint(-3, 3)) for _ in range(5)]
    return sp.Matrix([[a, b, c], [b, d, e], [c, e, -a - d]])


def _random_grad_e(rng: random.Random) -> list[sp.Matrix]:
    """Random integer STF-3 octupole (GradE), as its three k-slices.

    GradE = partial^3 Phi is totally symmetric (Schwarz) and vacuum
    trace-free: 7 free components, the rest fixed by the trace conditions.
    (Corrected 2026-07-12 from three independent STF-2 slices.)
    """
    c000, c001, c002, c011, c012, c111, c112 = [
        sp.Integer(rng.randint(-3, 3)) for _ in range(7)
    ]
    comp = {
        (0, 0, 0): c000, (0, 0, 1): c001, (0, 0, 2): c002,
        (0, 1, 1): c011, (0, 1, 2): c012, (0, 2, 2): -c000 - c011,
        (1, 1, 1): c111, (1, 1, 2): c112, (1, 2, 2): -c001 - c111,
        (2, 2, 2): -c002 - c112,
    }
    def entry(a: int, b: int, c: int) -> sp.Integer:
        return comp[tuple(sorted((a, b, c)))]
    return [
        sp.Matrix(3, 3, lambda i, j, k=k: entry(k, i, j))
        for k in range(3)
    ]


def _component(values: object, rank: int, indices: tuple[int, ...]) -> int:
    if rank == 2:
        if isinstance(values, sp.MatrixBase):
            return int(values[indices[0], indices[1]])
        if isinstance(values, list):
            return int(values[_index_map(2)[indices]])
    if rank == 6 and isinstance(values, list):
        return int(values[_index_map(6)[indices]])
    raise TypeError(f"Unsupported component lookup for rank {rank} and value type {type(values)}")


def _evaluate_label(
    sample: dict[str, object],
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> int:
    slot_to_edge: dict[tuple[int, int], int] = {}
    for edge_index, edge in enumerate(representative):
        slot_to_edge[edge[0]] = edge_index
        slot_to_edge[edge[1]] = edge_index

    total = 0
    for assignment in product(range(3), repeat=len(representative)):
        term = 1
        for instance_index, block_name in enumerate(signature):
            if block_name == "E":
                indices = tuple(assignment[slot_to_edge[(instance_index, slot)]] for slot in range(2))
                term *= _component(sample["E"], 2, indices)
            elif block_name in {"Z", "DtZ"}:
                indices = tuple(assignment[slot_to_edge[(instance_index, slot)]] for slot in range(6))
                term *= _component(sample[block_name], 6, indices)
            elif block_name == "GradE":
                derivative_index = assignment[slot_to_edge[(instance_index, 0)]]
                grad_e = sample["GradE"]
                assert isinstance(grad_e, list)
                indices = tuple(assignment[slot_to_edge[(instance_index, slot)]] for slot in range(1, 3))
                term *= _component(grad_e[derivative_index], 2, indices)
            elif block_name == "GradZ":
                derivative_index = assignment[slot_to_edge[(instance_index, 0)]]
                grad_z = sample["GradZ"]
                assert isinstance(grad_z, list)
                indices = tuple(assignment[slot_to_edge[(instance_index, slot)]] for slot in range(1, 7))
                term *= _component(grad_z[derivative_index], 6, indices)
            else:
                raise ValueError(f"Unsupported block name {block_name}")
        total += term
    return total


def _normalize_nullspace_vector(
    labels: tuple[str, ...],
    vector: sp.Matrix,
) -> tuple[tuple[int, ...], sp.Expr]:
    denominators = [value.q for value in vector if isinstance(value, sp.Rational) and value.q != 1]
    if not denominators:
        scale = 1
    elif len(denominators) == 1:
        scale = denominators[0]
    else:
        scale = sp.ilcm(*denominators)
    integers = [sp.Integer(scale) * value for value in vector]
    nonzero = [abs(int(value)) for value in integers if value != 0]
    gcd = sp.igcd(*nonzero) if nonzero else 1
    integers = [sp.Integer(value // gcd) for value in integers]
    first_nonzero = next((value for value in integers if value != 0), sp.Integer(1))
    if first_nonzero < 0:
        integers = [-value for value in integers]
    relation = sp.expand(sum(coefficient * sp.Symbol(label) for label, coefficient in zip(labels, integers)))
    return tuple(int(value) for value in integers), relation


def _sample_rows(
    group: tuple[tuple[str, ContractionClass], ...],
    sample_count: int,
    seed: int,
) -> list[list[int]]:
    rng = random.Random(seed)
    rows: list[list[int]] = []
    for _ in range(sample_count):
        sample = {
            "E": _random_stf2(rng),
            "GradE": _random_grad_e(rng),
            "Z": _random_stf(6, rng),
            "DtZ": _random_stf(6, rng),
            "GradZ": _random_grad_stf(6, rng),
        }
        rows.append(
            [
                _evaluate_label(sample, item.signature, item.representative)
                for _, item in group
            ]
        )
    return rows


def _verify_relations(
    group: tuple[tuple[str, ContractionClass], ...],
    relations: tuple[tuple[tuple[int, ...], sp.Expr], ...],
    seed: int,
) -> tuple[sp.Expr, ...]:
    rows = _sample_rows(group, sample_count=VERIFICATION_ROWS, seed=seed)
    verified: list[sp.Expr] = []
    for coefficients, relation in relations:
        for row in rows:
            if sum(coefficient * value for coefficient, value in zip(coefficients, row)) != 0:
                raise ValueError("Rank-6 sample-based null relation failed revalidation.")
        verified.append(relation)
    return tuple(verified)


def _signature_groups() -> tuple[tuple[tuple[str, ContractionClass], ...], ...]:
    groups: dict[tuple[str, ...], list[tuple[str, ContractionClass]]] = {}
    for item in enumerate_r6_family_classes():
        groups.setdefault(item.signature, []).append((item.label, item))
    return tuple(
        tuple(sorted(entries, key=lambda pair: pair[0]))
        for _, entries in sorted(groups.items(), key=lambda pair: (len(pair[0]), pair[0]))
    )


def r6_rank_summary() -> R6RankSummary:
    summary = r6_summary()
    baseline = rank_summary()
    signature_summaries: list[R6SectorSignatureRank] = []
    null_relations: list[sp.Expr] = []
    new_rank = 0
    group_seed = RNG_SEED

    for group in _signature_groups():
        labels = tuple(label for label, _ in group)
        sample_count = max(len(labels) + 3, 6)
        rows = _sample_rows(group, sample_count=sample_count, seed=group_seed)
        group_seed += 97
        matrix = sp.Matrix(rows)
        group_rank = matrix.rank()
        relations = tuple(
            _normalize_nullspace_vector(labels, vector)
            for vector in matrix.nullspace()
        )
        verified_relations = _verify_relations(group, relations, seed=group_seed)
        group_seed += 97
        new_rank += group_rank
        null_relations.extend(verified_relations)
        signature_summaries.append(
            R6SectorSignatureRank(
                signature=group[0][1].signature,
                labels=labels,
                sample_count=sample_count,
                rank=group_rank,
                count=len(labels),
                nullity=len(labels) - group_rank,
                first_null_relation=None if not verified_relations else verified_relations[0],
            )
        )

    total_rank = baseline.total_rank + new_rank
    total_count = len(summary.surviving_labels)
    total_nullity = total_count - total_rank
    return R6RankSummary(
        rank=total_rank,
        count=total_count,
        nullity=total_nullity,
        baseline_rank=baseline.total_rank,
        new_rank=new_rank,
        new_count=len(summary.new_surviving_labels),
        signature_ranks=tuple(signature_summaries),
        first_null_relation=None if not null_relations else null_relations[0],
        additional_null_relations=tuple(null_relations[1:]),
        labels=summary.surviving_labels,
        new_labels=summary.new_surviving_labels,
    )


def r6_survivor_rank_report() -> str:
    summary = r6_rank_summary()
    lines = [
        "Delta<=4 rank-6 family survivor rank audit",
        "",
        "Raw R6-extended survivor candidates:",
        "- " + ", ".join(summary.labels),
        "",
        "New rank-6 candidates beyond the enlarged audited-set baseline:",
        "- " + ", ".join(summary.new_labels),
        "",
        f"Baseline electric rank: {summary.baseline_rank}",
        f"New rank-6-sector sample-stable rank: {summary.new_rank} out of {summary.new_count}",
        f"Total sample-stable rank: {summary.rank} out of {summary.count}",
        f"Sample-stable nullity: {summary.nullity}",
        "",
        "Signature-group ranks:",
    ]
    for group in summary.signature_ranks:
        lines.append(
            f"- {group.signature}: rank {group.rank} / {group.count} using {group.sample_count} samples"
        )
        if group.first_null_relation is not None:
            lines.append(f"  first relation: {sp.sstr(group.first_null_relation)} = 0")

    if summary.first_null_relation is not None:
        lines.extend(
            [
                "",
                "First revalidated null relation:",
                f"- {sp.sstr(summary.first_null_relation)} = 0",
            ]
        )

    lines.extend(
        [
            "",
            "Interpretation:",
            "- The raw R6-extended candidate list is not being presented as a corrected basis statement.",
            "- The rank result is sample-stable rather than a closed symbolic basis theorem; it is used only as bookkeeping for raw-list dependence.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(r6_survivor_rank_report())


if __name__ == "__main__":
    main()
