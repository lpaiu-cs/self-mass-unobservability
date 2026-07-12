from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement, permutations, product

from enumerate_contractions_delta4 import ContractionClass
from r3_sector_delta4 import r3_surviving_classes
from r4_sector_delta4 import r4_surviving_classes
from r5_sector_delta4 import r5_surviving_classes


MAX_WEIGHT = 4
# Corrected electric baseline (2026-07-12): the gradient sector is the single
# invariant gradE2 (GradE is an STF-3 octupole; divE2 = 0, mixedGradE2 = gradE2).
BASELINE_SURVIVOR_LABELS = {
    "E2",
    "E3",
    "dotE2",
    "E2^2",
    "gradE2",
}


@dataclass(frozen=True)
class GeneratedHighRankClass:
    rank: int
    family_name: str
    signature: tuple[str, ...]
    weight: int
    key: tuple
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...]
    label: str
    matched_manual_label: str | None


@dataclass(frozen=True)
class HighRankAuditSummary:
    rank: int
    family_name: str
    generated_count: int
    manual_count: int
    matched_count: int
    omitted_from_manual: tuple[str, ...]
    manual_only: tuple[str, ...]


def family_symbol(rank: int) -> str:
    return {3: "T", 4: "Q", 5: "U", 6: "Z"}[rank]


def family_rank(rank: int, name: str) -> int:
    family_name = family_symbol(rank)
    if name == "E":
        return 2
    if name == family_name or name == f"Dt{family_name}":
        return rank
    if name == "GradE":
        return 3
    if name == f"Grad{family_name}":
        return rank + 1
    raise ValueError(f"Unhandled block name {name}")


def block_weight(name: str, rank: int) -> int:
    family_name = family_symbol(rank)
    if name in {"E", family_name}:
        return 1
    if name in {"GradE", f"Dt{family_name}", f"Grad{family_name}"}:
        return 2
    raise ValueError(f"Unhandled block name {name}")


def block_parity(name: str, rank: int) -> int:
    family_name = family_symbol(rank)
    if name in {"GradE", f"Grad{family_name}"}:
        return 1
    return 0


def available_block_names(rank: int) -> tuple[str, ...]:
    family_name = family_symbol(rank)
    return ("E", family_name, f"Dt{family_name}", "GradE", f"Grad{family_name}")


def signature_in_name_order(rank: int, names: tuple[str, ...]) -> tuple[str, ...]:
    order = {name: index for index, name in enumerate(available_block_names(rank))}
    return tuple(sorted(names, key=lambda name: order[name]))


def is_normal_form_candidate_signature(rank: int, signature: tuple[str, ...]) -> bool:
    family_name = family_symbol(rank)
    dt_family = f"Dt{family_name}"
    grad_family = f"Grad{family_name}"

    if family_name not in signature and dt_family not in signature and grad_family not in signature:
        return False

    if sum(block_weight(name, rank) for name in signature) > MAX_WEIGHT:
        return False
    if sum(block_parity(name, rank) for name in signature) % 2:
        return False
    if sum(family_rank(rank, name) for name in signature) % 2:
        return False

    if "GradE" in signature or grad_family in signature:
        return signature in {
            signature_in_name_order(rank, (grad_family, grad_family)),
            signature_in_name_order(rank, ("GradE", grad_family)),
        }

    if dt_family in signature:
        allowed = {
            signature_in_name_order(rank, (dt_family, dt_family)),
            signature_in_name_order(rank, ("E", dt_family, family_name)),
            signature_in_name_order(rank, ("E", "E", dt_family)),
        }
        return signature in allowed

    return all(name in {"E", family_name} for name in signature)


def candidate_signatures(rank: int) -> tuple[tuple[str, ...], ...]:
    names = available_block_names(rank)
    out: list[tuple[str, ...]] = []
    for size in range(1, 5):
        for signature in combinations_with_replacement(names, size):
            if is_normal_form_candidate_signature(rank, signature):
                out.append(signature)
    return tuple(sorted(out, key=lambda item: (sum(block_weight(name, rank) for name in item), item)))


def signature_permutation_to_name_order(rank: int, signature: tuple[str, ...]) -> tuple[int, ...]:
    order = {name: index for index, name in enumerate(available_block_names(rank))}
    return tuple(sorted(range(len(signature)), key=lambda index: (order[signature[index]], index)))


def reorder_signature_and_representative(
    rank: int,
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[tuple[str, ...], tuple[tuple[tuple[int, int], tuple[int, int]], ...]]:
    perm = signature_permutation_to_name_order(rank, signature)
    inverse = {old: new for new, old in enumerate(perm)}
    reordered_signature = tuple(signature[old] for old in perm)
    reordered_edges = []
    for (left_inst, left_slot), (right_inst, right_slot) in representative:
        lhs = (inverse[left_inst], left_slot)
        rhs = (inverse[right_inst], right_slot)
        if rhs < lhs:
            lhs, rhs = rhs, lhs
        reordered_edges.append((lhs, rhs))
    reordered_edges.sort()
    return reordered_signature, tuple(reordered_edges)


def identical_group_permutations(signature: tuple[str, ...]) -> tuple[tuple[int, ...], ...]:
    groups: dict[str, list[int]] = {}
    for index, name in enumerate(signature):
        groups.setdefault(name, []).append(index)
    per_group = [tuple(permutations(indices)) for _, indices in sorted(groups.items())]
    out = []
    for choice in product(*per_group):
        perm = list(range(len(signature)))
        for indices, mapping in zip((groups[name] for name in sorted(groups)), choice):
            for old, new in zip(indices, mapping):
                perm[old] = new
        out.append(tuple(perm))
    return tuple(sorted(set(out)))


def multigraph_solutions(degrees: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    edges = tuple(combinations(range(len(degrees)), 2))
    remaining = list(degrees)
    current = [0] * len(edges)
    out: list[tuple[int, ...]] = []

    def search(edge_index: int) -> None:
        if edge_index == len(edges):
            if all(value == 0 for value in remaining):
                out.append(tuple(current))
            return
        left, right = edges[edge_index]
        upper = min(remaining[left], remaining[right])
        for multiplicity in range(upper + 1):
            remaining[left] -= multiplicity
            remaining[right] -= multiplicity
            current[edge_index] = multiplicity
            search(edge_index + 1)
            current[edge_index] = 0
            remaining[left] += multiplicity
            remaining[right] += multiplicity

    search(0)
    return tuple(out)


def multiplicities_to_matrix(count: int, multiplicities: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    matrix = [[0] * count for _ in range(count)]
    for multiplicity, (left, right) in zip(multiplicities, combinations(range(count), 2)):
        matrix[left][right] = multiplicity
        matrix[right][left] = multiplicity
    return tuple(tuple(row) for row in matrix)


def matrix_to_multiplicities(matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    return tuple(matrix[left][right] for left, right in combinations(range(len(matrix)), 2))


def canonical_pure_stf_key(signature: tuple[str, ...], multiplicities: tuple[int, ...]) -> tuple:
    matrix = multiplicities_to_matrix(len(signature), multiplicities)
    candidates = []
    for perm in identical_group_permutations(signature):
        permuted = tuple(
            tuple(matrix[perm[left]][perm[right]] for right in range(len(signature)))
            for left in range(len(signature))
        )
        candidates.append(matrix_to_multiplicities(permuted))
    return ("pure", signature, min(candidates))


def representative_from_pure_key(
    signature: tuple[str, ...],
    multiplicities: tuple[int, ...],
) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    ranks = tuple(2 if name == "E" else family_rank_from_signature_name(name) for name in signature)
    cursors = [0] * len(signature)
    edges_out: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for multiplicity, (left, right) in zip(multiplicities, combinations(range(len(signature)), 2)):
        for _ in range(multiplicity):
            edges_out.append(((left, cursors[left]), (right, cursors[right])))
            cursors[left] += 1
            cursors[right] += 1
    if tuple(cursors) != ranks:
        raise ValueError(f"Degree bookkeeping drifted for {signature}: {tuple(cursors)} != {ranks}")
    return tuple(edges_out)


def family_rank_from_signature_name(name: str) -> int:
    if name == "E":
        return 2
    if name.startswith("Dt") and not name.startswith("DtE"):
        return {"DtT": 3, "DtQ": 4, "DtU": 5, "DtZ": 6}[name]
    if name in {"T", "Q", "U", "Z"}:
        return {"T": 3, "Q": 4, "U": 5, "Z": 6}[name]
    raise ValueError(f"Unhandled pure STF signature name {name}")


def pure_stf_generated_classes(rank: int, signature: tuple[str, ...]) -> tuple[GeneratedHighRankClass, ...]:
    family_name = family_symbol(rank)
    degrees = tuple(family_rank(rank, name) for name in signature)
    classes: list[GeneratedHighRankClass] = []
    seen: set[tuple] = set()
    manual_map = manual_label_map(rank)
    base_label = signature_base_label(signature, family_name)

    for multiplicities in multigraph_solutions(degrees):
        key = canonical_pure_stf_key(signature, multiplicities)
        if key in seen:
            continue
        seen.add(key)
    sorted_keys = sorted(seen)
    needs_suffix = len(sorted_keys) > 1
    for index, key in enumerate(sorted_keys, start=1):
        _, key_signature, canonical_multiplicities = key
        representative = representative_from_pure_key(key_signature, canonical_multiplicities)
        matched = manual_map.get(key)
        label = matched or (base_label if not needs_suffix else f"{base_label}_{index}")
        classes.append(
            GeneratedHighRankClass(
                rank=rank,
                family_name=family_name,
                signature=signature,
                weight=sum(block_weight(name, rank) for name in signature),
                key=key,
                representative=representative,
                label=label,
                matched_manual_label=matched,
            )
        )
    return tuple(classes)


def gradient_profile(rank: int, name: str) -> tuple[int, bool]:
    family_name = family_symbol(rank)
    if name == "GradE":
        return 2, False
    if name == f"Grad{family_name}":
        return rank, True
    raise ValueError(f"Unhandled gradient name {name}")


def gradient_descriptor(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[int, int, int, int, int, int]:
    internal_left = 0
    internal_right = 0
    dd = 0
    d_left_to_right = 0
    d_right_to_left = 0
    cross_stf = 0
    for left, right in representative:
        if left[0] == right[0] == 0:
            internal_left += 1
        elif left[0] == right[0] == 1:
            internal_right += 1
        elif {left, right} == {((0, 0)), ((1, 0))}:  # type: ignore[arg-type]
            dd += 1
        else:
            a, b = left, right
            if a[0] == 0 and a[1] == 0 and b[0] == 1 and b[1] > 0:
                d_left_to_right += 1
            elif b[0] == 0 and b[1] == 0 and a[0] == 1 and a[1] > 0:
                d_left_to_right += 1
            elif a[0] == 1 and a[1] == 0 and b[0] == 0 and b[1] > 0:
                d_right_to_left += 1
            elif b[0] == 1 and b[1] == 0 and a[0] == 0 and a[1] > 0:
                d_right_to_left += 1
            else:
                cross_stf += 1
    return internal_left, internal_right, dd, d_left_to_right, d_right_to_left, cross_stf


def canonical_gradient_key(signature: tuple[str, ...], descriptor: tuple[int, int, int, int, int, int]) -> tuple:
    if signature[0] == signature[1]:
        il, ir, dd, dlr, drl, s = descriptor
        swapped = (ir, il, dd, drl, dlr, s)
        descriptor = min(descriptor, swapped)
    return ("gradient", signature, descriptor)


def representative_from_gradient_descriptor(
    signature: tuple[str, ...],
    descriptor: tuple[int, int, int, int, int, int],
    rank: int,
) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    left_stf_slots, _ = gradient_profile(rank, signature[0])
    right_stf_slots, _ = gradient_profile(rank, signature[1])
    internal_left, internal_right, dd, d_left_to_right, d_right_to_left, cross_stf = descriptor

    left_available = list(range(1, left_stf_slots + 1))
    right_available = list(range(1, right_stf_slots + 1))
    edges_out: list[tuple[tuple[int, int], tuple[int, int]]] = []

    if internal_left:
        edges_out.append(((0, 0), (0, left_available.pop(0))))
    if internal_right:
        edges_out.append(((1, 0), (1, right_available.pop(0))))
    if dd:
        edges_out.append(((0, 0), (1, 0)))
    if d_left_to_right:
        edges_out.append(((0, 0), (1, right_available.pop(0))))
    if d_right_to_left:
        edges_out.append(((1, 0), (0, left_available.pop(0))))
    for _ in range(cross_stf):
        edges_out.append(((0, left_available.pop(0)), (1, right_available.pop(0))))

    if left_available or right_available:
        raise ValueError(
            f"Gradient descriptor did not consume all STF slots for {signature}: "
            f"{left_available}, {right_available}"
        )
    return tuple(edges_out)


def gradient_descriptors(rank: int, signature: tuple[str, ...]) -> tuple[tuple[int, int, int, int, int, int], ...]:
    left_stf_slots, _ = gradient_profile(rank, signature[0])
    right_stf_slots, _ = gradient_profile(rank, signature[1])
    out = set()
    for internal_left in (0, 1):
        for internal_right in (0, 1):
            for dd in (0, 1):
                for d_left_to_right in (0, 1):
                    for d_right_to_left in (0, 1):
                        if internal_left + dd + d_left_to_right != 1:
                            continue
                        if internal_right + dd + d_right_to_left != 1:
                            continue
                        cross_stf_left = left_stf_slots - internal_left - d_right_to_left
                        cross_stf_right = right_stf_slots - internal_right - d_left_to_right
                        if cross_stf_left < 0 or cross_stf_right < 0:
                            continue
                        if cross_stf_left != cross_stf_right:
                            continue
                        out.add((internal_left, internal_right, dd, d_left_to_right, d_right_to_left, cross_stf_left))
    keys = {canonical_gradient_key(signature, descriptor) for descriptor in out}
    return tuple(sorted(key[2] for key in keys))


def gradient_generated_classes(rank: int, signature: tuple[str, ...]) -> tuple[GeneratedHighRankClass, ...]:
    family_name = family_symbol(rank)
    descriptors = gradient_descriptors(rank, signature)
    manual_map = manual_label_map(rank)
    base_label = signature_base_label(signature, family_name)
    needs_suffix = len(descriptors) > 1
    out: list[GeneratedHighRankClass] = []
    for index, descriptor in enumerate(descriptors, start=1):
        key = canonical_gradient_key(signature, descriptor)
        representative = representative_from_gradient_descriptor(signature, descriptor, rank)
        matched = manual_map.get(key)
        label = matched or (base_label if not needs_suffix else f"{base_label}_{index}")
        out.append(
            GeneratedHighRankClass(
                rank=rank,
                family_name=family_name,
                signature=signature,
                weight=sum(block_weight(name, rank) for name in signature),
                key=key,
                representative=representative,
                label=label,
                matched_manual_label=matched,
            )
        )
    return tuple(out)


def generated_family_classes(rank: int) -> tuple[GeneratedHighRankClass, ...]:
    out: list[GeneratedHighRankClass] = []
    for signature in candidate_signatures(rank):
        if any(name.startswith("Grad") for name in signature):
            out.extend(gradient_generated_classes(rank, signature))
        else:
            out.extend(pure_stf_generated_classes(rank, signature))
    return tuple(sorted(out, key=lambda item: (item.weight, item.signature, item.label)))


def manual_survivor_classes(rank: int) -> tuple[ContractionClass, ...]:
    by_rank = {
        3: r3_surviving_classes,
        4: r4_surviving_classes,
        5: r5_surviving_classes,
    }
    if rank not in by_rank:
        return ()
    return tuple(
        item for item in by_rank[rank]() if item.label not in BASELINE_SURVIVOR_LABELS
    )


def manual_class_key(rank: int, item: ContractionClass) -> tuple:
    signature, representative = reorder_signature_and_representative(rank, item.signature, item.representative)
    if any(name.startswith("Grad") for name in signature):
        return canonical_gradient_key(signature, gradient_descriptor(signature, representative))
    count = len(signature)
    matrix = [[0] * count for _ in range(count)]
    for (left_inst, _), (right_inst, _) in representative:
        matrix[left_inst][right_inst] += 1
        matrix[right_inst][left_inst] += 1
    return canonical_pure_stf_key(signature, matrix_to_multiplicities(tuple(tuple(row) for row in matrix)))


def manual_label_map(rank: int) -> dict[tuple, str]:
    # The completed manual survivor lists may use blocks outside this
    # enumerator's candidate universe (DtE entered at the rank-4
    # completion); such classes can never match a generated class, so they
    # are keyed out here and surface through manual_only in the audit.
    universe = set(available_block_names(rank))
    mapping: dict[tuple, str] = {}
    for item in manual_survivor_classes(rank):
        if not set(item.signature) <= universe:
            continue
        mapping[manual_class_key(rank, item)] = item.label
    return mapping


def signature_base_label(signature: tuple[str, ...], family_name: str) -> str:
    dt_family = f"Dt{family_name}"
    grad_family = f"Grad{family_name}"
    table = {
        (family_name, family_name): f"{family_name}2",
        ("E", family_name, family_name): f"E{family_name}{family_name}",
        ("E", "E", family_name): f"EE{family_name}",
        (family_name, family_name, family_name): f"{family_name}3",
        (dt_family, dt_family): f"dot{family_name}2",
        ("E", family_name, dt_family): f"E{family_name}Dt{family_name}",
        ("E", "E", dt_family): f"EEDt{family_name}",
        ("E", "E", family_name, family_name): f"E2{family_name}2",
        (family_name, family_name, family_name, family_name): f"{family_name}4",
        ("E", family_name, family_name, family_name): f"E{family_name}{family_name}{family_name}",
        ("E", "E", "E", family_name): f"EEE{family_name}",
        (grad_family, grad_family): f"Grad{family_name}2",
        ("GradE", grad_family): f"GradEGrad{family_name}",
    }
    return table.get(signature, "_".join(signature))


def high_rank_audit_summary(rank: int) -> HighRankAuditSummary:
    generated = generated_family_classes(rank)
    generated_labels = {item.label for item in generated}
    manual_labels = {item.label for item in manual_survivor_classes(rank)}
    matched = {item.matched_manual_label for item in generated if item.matched_manual_label is not None}
    return HighRankAuditSummary(
        rank=rank,
        family_name=family_symbol(rank),
        generated_count=len(generated_labels),
        manual_count=len(manual_labels),
        matched_count=len(matched),
        omitted_from_manual=tuple(sorted(generated_labels - manual_labels)),
        manual_only=tuple(sorted(manual_labels - generated_labels)),
    )


def high_rank_enumerator_report() -> str:
    lines = ["rank\tfamily\tgenerated_count\tmanual_count\tmatched_count\tomitted_from_manual\tmanual_only"]
    for rank in (3, 4, 5):
        summary = high_rank_audit_summary(rank)
        lines.append(
            "\t".join(
                (
                    str(summary.rank),
                    summary.family_name,
                    str(summary.generated_count),
                    str(summary.manual_count),
                    str(summary.matched_count),
                    ",".join(summary.omitted_from_manual) if summary.omitted_from_manual else "-",
                    ",".join(summary.manual_only) if summary.manual_only else "-",
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(high_rank_enumerator_report())


if __name__ == "__main__":
    main()
