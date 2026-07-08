from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import random

import numpy as np
import sympy as sp

from enumerate_contractions_delta4 import enumerate_contraction_classes
from r4_sector_delta4 import r4_summary, r4_surviving_classes
from survivor_rank_check import rank_summary


SAMPLE_COUNT = 30
RNG_SEED = 424242

# ---------------------------------------------------------------------------
# Corrected rank-4 survivor dimension (exact O(3) character integral).
#
# The hand-built candidate list in `_r4_new_specs()` below caps the mixed E/Q
# sector at degree 2 in E (`E2Q2`) and has no E/Q cross-gradient, so it
# under-counts the rank-4 survivors as 19. Under the theorem's OWN reduction
# rules (total derivative, lower-order EOM, and the STF algebraic identities,
# all of which the O(3) character dimension already respects), the correct
# rank-4 survivor dimension is 25. The omitted survivors are the higher-degree
# mixed operators such as `Q_abcd (E^2)_ab E_cd` (degree 3 in E), `E Q^3`
# (degree 3 in Q), and the `GradE.GradQ` cross-gradient. See
# `verification/rederive_rank4.py` for their explicit construction and
# nonzero/rotation-invariance checks; two independent methods (this character
# integral and a delta-only contraction enumerator) agree on 25.
# ---------------------------------------------------------------------------
CORRECTED_SURVIVOR_DIMENSION = 25


def _character_survivor_dimension(family_rank: int = 4, wmax: int = 4) -> int:
    """Exact parity-even, delta-only survivor dimension of the E/(rank-r) sector
    via O(3) character integrals:  survivor(w) = inv_trunc(w) - inv_prom(w-1),
    with inv over the truncated block set {E,DtE,Dt2E,GradE}+{X,DtX,Dt2X,GradX}
    and inv_prom over the D_tau-promotable subset {E,DtE}+{X,DtX}. A delta-only
    scalar needs an even total Cartesian index count (else it would need an
    epsilon -> pseudoscalar, excluded). For a parity-even family (this function's
    case) it reproduces the audited survivor rank exactly at family ranks
    1,3,5,6 (V=17, T=19, U=19, Z=23) and the electric baseline (7); at rank 4
    (Q) it gives 25, the corrected value. (The magnetic rank-2 family is
    parity-odd, a separate case handled by the eb sector.)"""
    n = 4000
    theta = (np.arange(n) + 0.5) * np.pi / n
    measure = (1.0 - np.cos(theta)) / np.pi * (np.pi / n)

    def chi(irreps, k):
        s = np.zeros_like(theta)
        for l in irreps:
            s = s + np.sin((2 * l + 1) * (k * theta) / 2) / np.sin((k * theta) / 2)
        return s

    def sym_power(pk, m):
        h = [np.ones_like(theta)]
        for j in range(1, m + 1):
            acc = np.zeros_like(theta)
            for k in range(1, j + 1):
                acc = acc + pk[k - 1] * h[j - k]
            h.append(acc / j)
        return h[m]

    def grad_irreps(r):
        return tuple(sorted({abs(r - 1), r, r + 1})) if r >= 1 else (1,)

    # block -> (irreps, parity, weight, cartesian_rank)
    def blocks(sym, r, par):
        return {
            sym: ((r,), par, 1, r),
            "Dt" + sym: ((r,), par, 2, r),
            "Dt2" + sym: ((r,), par, 3, r),
            "Grad" + sym: (grad_irreps(r), -par, 2, r + 1),
        }

    B = {**blocks("E", 2, +1), **blocks("Q", family_rank, +1)}

    def sig_dim(sig):
        if sum(B[nm][3] * c for nm, c in sig.items()) % 2 == 1:  # odd Cartesian index count
            return 0
        if sum(c for nm, c in sig.items() if B[nm][1] == -1) % 2 == 1:  # parity-odd
            return 0
        prod = np.ones_like(theta)
        for nm, c in sig.items():
            prod = prod * sym_power([chi(B[nm][0], k) for k in range(1, c + 1)], c)
        return round(float(np.sum(prod * measure)))

    def inv_dim(names, w):
        caps = [range(0, w // B[nm][2] + 1) for nm in names]
        total = 0
        for counts in product(*caps):
            if sum(c * B[nm][2] for c, nm in zip(counts, names)) == w and sum(counts) >= 1:
                total += sig_dim({nm: c for nm, c in zip(names, counts) if c > 0})
        return total

    trunc = ["E", "DtE", "Dt2E", "GradE", "Q", "DtQ", "Dt2Q", "GradQ"]
    prom = ["E", "DtE", "Q", "DtQ"]
    survivor = inv_dim(trunc, 1)
    for w in range(2, wmax + 1):
        survivor += inv_dim(trunc, w) - inv_dim(prom, w - 1)
    return survivor


@dataclass(frozen=True)
class R4RankSummary:
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


def _rank4_entry_order() -> tuple[tuple[int, int, int, int], ...]:
    return tuple(
        (i, j, k, l)
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for l in range(3)
    )


def _stf_rank4_basis() -> tuple[sp.MutableDenseNDimArray, ...]:
    patterns4 = _count_patterns(4)
    pattern_index4 = {pattern: idx for idx, pattern in enumerate(patterns4)}
    patterns2 = _count_patterns(2)
    constraint_rows: list[list[int]] = []
    for pattern2 in patterns2:
        ij = tuple([0] * pattern2[0] + [1] * pattern2[1] + [2] * pattern2[2])
        row = [0] * len(patterns4)
        for k in range(3):
            row[pattern_index4[_counts((k, k) + ij)]] += 1
        constraint_rows.append(row)
    nullspace = sp.Matrix(constraint_rows).nullspace()
    entries4 = _rank4_entry_order()
    basis: list[sp.MutableDenseNDimArray] = []
    for vector in nullspace:
        array = sp.MutableDenseNDimArray.zeros(3, 3, 3, 3)
        for indices in entries4:
            array[indices] = vector[pattern_index4[_counts(indices)]]
        basis.append(array)
    return tuple(basis)


def _symbolic_stf_rank4(prefix: str) -> tuple[sp.MutableDenseNDimArray, tuple[sp.Symbol, ...]]:
    basis = _stf_rank4_basis()
    coefficients = sp.symbols(f"{prefix}_0:{len(basis)}")
    entries4 = _rank4_entry_order()
    array = sp.MutableDenseNDimArray.zeros(3, 3, 3, 3)
    for coeff, basis_array in zip(coefficients, basis):
        for indices in entries4:
            array[indices] += coeff * basis_array[indices]
    return array, tuple(coefficients)


def _random_stf2(rng: random.Random) -> sp.Matrix:
    a, b, c, d, e = [sp.Integer(rng.randint(-3, 3)) for _ in range(5)]
    return sp.Matrix([[a, b, c], [b, d, e], [c, e, -a - d]])


def _random_stf4(rng: random.Random) -> list[int]:
    basis = _stf_rank4_basis()
    entries4 = _rank4_entry_order()
    coefficients = [rng.randint(-2, 2) for _ in range(len(basis))]
    while all(value == 0 for value in coefficients):
        coefficients = [rng.randint(-2, 2) for _ in range(len(basis))]
    values = [0] * len(entries4)
    for coeff, basis_array in zip(coefficients, basis):
        for index, indices in enumerate(entries4):
            values[index] += coeff * int(basis_array[indices])
    return values


def _random_grad_stf4(rng: random.Random) -> list[list[int]]:
    return [_random_stf4(rng) for _ in range(3)]


def _r4_new_specs() -> dict[str, tuple[tuple[str, ...], tuple[tuple[tuple[int, int], tuple[int, int]], ...]]]:
    return {
        "Q2": (
            ("Q", "Q"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
            ),
        ),
        "EQQ": (
            ("E", "Q", "Q"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
                ((1, 3), (2, 3)),
            ),
        ),
        "dotQ2": (
            ("DtQ", "DtQ"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
            ),
        ),
        "EQDtQ": (
            ("DtQ", "E", "Q"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((0, 2), (2, 1)),
                ((0, 3), (2, 2)),
                ((1, 1), (2, 3)),
            ),
        ),
        "E2Q2": (
            ("E", "E", "Q", "Q"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((2, 0), (3, 0)),
                ((2, 1), (3, 1)),
                ((2, 2), (3, 2)),
                ((2, 3), (3, 3)),
            ),
        ),
        "E2Q2_mixed_1": (
            ("E", "E", "Q", "Q"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((1, 1), (3, 0)),
                ((2, 1), (3, 1)),
                ((2, 2), (3, 2)),
                ((2, 3), (3, 3)),
            ),
        ),
        "E2Q2_mixed_2": (
            ("E", "E", "Q", "Q"),
            (
                ((0, 0), (2, 0)),
                ((0, 1), (2, 1)),
                ((1, 0), (3, 0)),
                ((1, 1), (3, 1)),
                ((2, 2), (3, 2)),
                ((2, 3), (3, 3)),
            ),
        ),
        "E2Q2_mixed_3": (
            ("E", "E", "Q", "Q"),
            (
                ((0, 0), (2, 0)),
                ((0, 1), (3, 0)),
                ((1, 0), (2, 1)),
                ((1, 1), (3, 1)),
                ((2, 2), (3, 2)),
                ((2, 3), (3, 3)),
            ),
        ),
        "divQ2": (
            ("GradQ", "GradQ"),
            (
                ((0, 0), (0, 1)),
                ((0, 2), (1, 1)),
                ((0, 3), (1, 2)),
                ((0, 4), (1, 3)),
                ((1, 0), (1, 4)),
            ),
        ),
        "gradQ2": (
            ("GradQ", "GradQ"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
            ),
        ),
        "mixedGradQ2": (
            ("GradQ", "GradQ"),
            (
                ((0, 0), (1, 1)),
                ((0, 1), (1, 0)),
                ((0, 2), (1, 2)),
                ((0, 3), (1, 3)),
                ((0, 4), (1, 4)),
            ),
        ),
        "Q2^2": (
            ("Q", "Q", "Q", "Q"),
            (
                ((0, 0), (3, 0)),
                ((0, 1), (3, 1)),
                ((0, 2), (3, 2)),
                ((0, 3), (3, 3)),
                ((1, 0), (2, 0)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
                ((1, 3), (2, 3)),
            ),
        ),
        "Q4_chain": (
            ("Q", "Q", "Q", "Q"),
            (
                ((0, 0), (2, 0)),
                ((0, 1), (3, 0)),
                ((0, 2), (3, 1)),
                ((0, 3), (3, 2)),
                ((1, 0), (2, 1)),
                ((1, 1), (2, 2)),
                ((1, 2), (2, 3)),
                ((1, 3), (3, 3)),
            ),
        ),
        "Q4_bridge": (
            ("Q", "Q", "Q", "Q"),
            (
                ((0, 0), (2, 0)),
                ((0, 1), (2, 1)),
                ((0, 2), (3, 0)),
                ((0, 3), (3, 1)),
                ((1, 0), (2, 2)),
                ((1, 1), (2, 3)),
                ((1, 2), (3, 2)),
                ((1, 3), (3, 3)),
            ),
        ),
        "Q4_tetra": (
            ("Q", "Q", "Q", "Q"),
            (
                ((0, 0), (1, 0)),
                ((0, 1), (2, 0)),
                ((0, 2), (3, 0)),
                ((0, 3), (3, 1)),
                ((1, 1), (2, 1)),
                ((1, 2), (2, 2)),
                ((1, 3), (3, 2)),
                ((2, 3), (3, 3)),
            ),
        ),
    }


def _numeric_spec_meta() -> dict[str, tuple[tuple[str, ...], dict[tuple[int, int], int], int]]:
    out: dict[str, tuple[tuple[str, ...], dict[tuple[int, int], int], int]] = {}
    for label, (signature, representative) in _r4_new_specs().items():
        slot_to_edge: dict[tuple[int, int], int] = {}
        for edge_index, edge in enumerate(representative):
            slot_to_edge[edge[0]] = edge_index
            slot_to_edge[edge[1]] = edge_index
        out[label] = (signature, slot_to_edge, len(representative))
    return out


def _rank4_index() -> dict[tuple[int, int, int, int], int]:
    return {indices: idx for idx, indices in enumerate(_rank4_entry_order())}


def _q_component(values: list[int], indices: tuple[int, int, int, int]) -> int:
    return values[_rank4_index()[indices]]


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
            elif block_name in {"Q", "DtQ"}:
                indices = tuple(
                    assignment[slot_to_edge[(instance_index, slot)]] for slot in range(4)
                )
                values = sample[block_name]
                assert isinstance(values, list)
                term *= _q_component(values, indices)
            elif block_name == "GradQ":
                indices = tuple(
                    assignment[slot_to_edge[(instance_index, slot)]] for slot in range(5)
                )
                values = sample["GradQ"]
                assert isinstance(values, list)
                term *= values[indices[0]][_rank4_index()[indices[1:]]]
            else:
                raise ValueError(f"Unsupported block name {block_name}")
        total += term
    return total


def _new_sector_rank_lower_bound(sample_count: int = SAMPLE_COUNT) -> int:
    metadata = _numeric_spec_meta()
    labels = tuple(_r4_new_specs())
    rows: list[list[int]] = []
    rng = random.Random(RNG_SEED)
    for _ in range(sample_count):
        sample = {
            "E": _random_stf2(rng),
            "Q": _random_stf4(rng),
            "DtQ": _random_stf4(rng),
            "GradQ": _random_grad_stf4(rng),
        }
        rows.append([_eval_new_label(sample, label, metadata) for label in labels])
    return sp.Matrix(rows).rank()


def _component(
    tensor: object,
    indices: tuple[int, ...],
) -> sp.Expr:
    if isinstance(tensor, sp.MatrixBase):
        return tensor[indices[0], indices[1]]
    if isinstance(tensor, sp.NDimArray):
        return tensor[indices]
    raise TypeError(f"Unsupported tensor container {type(tensor)}")


def _symbolic_rep_expression(
    tensor_map: dict[str, object],
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> sp.Expr:
    slot_to_edge: dict[tuple[int, int], int] = {}
    for edge_index, edge in enumerate(representative):
        slot_to_edge[edge[0]] = edge_index
        slot_to_edge[edge[1]] = edge_index

    total = sp.Integer(0)
    for assignment in product(range(3), repeat=len(representative)):
        term = sp.Integer(1)
        for instance_index, block_name in enumerate(signature):
            rank = {"E": 2, "Q": 4}[block_name]
            indices = tuple(assignment[slot_to_edge[(instance_index, slot)]] for slot in range(rank))
            term *= _component(tensor_map[block_name], indices)
        total += term
    return sp.expand(total)


def _verified_null_relations() -> tuple[sp.Expr, ...]:
    electric = sp.Matrix(
        [
            [sp.Symbol("a"), sp.Symbol("b"), sp.Symbol("c")],
            [sp.Symbol("b"), sp.Symbol("d"), sp.Symbol("e")],
            [sp.Symbol("c"), sp.Symbol("e"), -sp.Symbol("a") - sp.Symbol("d")],
        ]
    )
    tidal, _ = _symbolic_stf_rank4("Q")
    tensor_map = {"E": electric, "Q": tidal}
    specs = _r4_new_specs()
    e2q2 = _symbolic_rep_expression(tensor_map, *specs["E2Q2"])
    e2q2_m1 = _symbolic_rep_expression(tensor_map, *specs["E2Q2_mixed_1"])
    e2q2_m2 = _symbolic_rep_expression(tensor_map, *specs["E2Q2_mixed_2"])
    e2q2_m3 = _symbolic_rep_expression(tensor_map, *specs["E2Q2_mixed_3"])
    q2_sq = _symbolic_rep_expression(tensor_map, *specs["Q2^2"])
    q4_chain = _symbolic_rep_expression(tensor_map, *specs["Q4_chain"])
    q4_bridge = _symbolic_rep_expression(tensor_map, *specs["Q4_bridge"])
    q4_tetra = _symbolic_rep_expression(tensor_map, *specs["Q4_tetra"])

    relation1 = -sp.Rational(1, 2) * sp.Symbol("E2Q2")
    relation1 += 2 * sp.Symbol("E2Q2_mixed_1")
    relation1 -= sp.Symbol("E2Q2_mixed_2")
    relation1 += sp.Symbol("E2Q2_mixed_3")
    relation2 = sp.Rational(3, 10) * sp.Symbol("Q2^2")
    relation2 += sp.Symbol("Q4_bridge")
    relation2 -= sp.Rational(8, 5) * sp.Symbol("Q4_chain")
    relation3 = -sp.Rational(1, 5) * sp.Symbol("Q2^2")
    relation3 += sp.Rational(2, 5) * sp.Symbol("Q4_chain")
    relation3 += sp.Symbol("Q4_tetra")

    check1 = -sp.Rational(1, 2) * e2q2 + 2 * e2q2_m1 - e2q2_m2 + e2q2_m3
    check2 = sp.Rational(3, 10) * q2_sq + q4_bridge - sp.Rational(8, 5) * q4_chain
    check3 = -sp.Rational(1, 5) * q2_sq + sp.Rational(2, 5) * q4_chain + q4_tetra
    if sp.simplify(check1) != 0 or sp.simplify(check2) != 0 or sp.simplify(check3) != 0:
        raise ValueError("Rank-4 null-relation verification failed.")
    return (sp.expand(relation1), sp.expand(relation2), sp.expand(relation3))


def r4_rank_summary() -> R4RankSummary:
    summary = r4_summary()
    baseline = rank_summary()
    null_relations = _verified_null_relations()
    new_rank_lower_bound = _new_sector_rank_lower_bound()
    new_count = len(summary.new_surviving_labels)
    new_nullity = len(null_relations)
    if new_rank_lower_bound != new_count - new_nullity:
        raise ValueError(
            "Rank-4 new-sector lower bound does not match the verified null-relation upper bound."
        )
    total_rank = baseline.total_rank + new_rank_lower_bound
    total_count = len(summary.surviving_labels)
    total_nullity = total_count - total_rank
    return R4RankSummary(
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


def r4_survivor_rank_report() -> str:
    summary = r4_rank_summary()
    corrected = _character_survivor_dimension()
    lines = [
        "Delta<=4 rank-4 family survivor rank audit",
        "",
        "Raw R4-extended survivor candidates:",
        "- " + ", ".join(summary.labels),
        "",
        "New candidates beyond the enlarged audited-set baseline:",
        "- " + ", ".join(summary.new_labels),
        "",
        f"Baseline electric rank: {summary.baseline_rank}",
        f"New rank-4-sector rank: {summary.new_rank} out of {summary.new_count}",
        f"Total rank (hand-built candidate list): {summary.rank} out of {summary.count}",
        f"Nullity: {summary.nullity}",
        f"Deterministic evaluation sample count: {summary.sample_count}",
        "",
        "*** CORRECTED total survivor dimension (exact O(3) character integral, no",
        f"    E-degree cap): {corrected}.  The hand-built candidate list above under-",
        "    counts because it caps the mixed E/Q sector at degree 2 in E (E2Q2) and",
        "    omits the E/Q cross-gradient; the missing higher-degree survivors are",
        "    constructed and verified in verification/rederive_rank4.py. ***",
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
            "- The raw R4-extended survivor list is not a corrected basis statement.",
            "- The rank result uses two ingredients: an exact lower bound from deterministic component evaluations and exact symbolic verification of the extracted null relations.",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    print(r4_survivor_rank_report())


if __name__ == "__main__":
    main()
