from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement, permutations, product
from typing import Callable


MAX_WEIGHT = 4


@dataclass(frozen=True)
class BlockType:
    name: str
    weight: int
    rank: int
    parity: int = 0
    sym_groups: tuple[tuple[int, ...], ...] = ()
    tracefree_pairs: tuple[tuple[int, int], ...] = ()


@dataclass(frozen=True)
class ContractionClass:
    signature: tuple[str, ...]
    label: str
    classification: str
    status: str
    weight: int
    reduction_channel: str
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...]


BLOCK_TYPES = (
    BlockType("E", 1, 2, parity=0, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
    BlockType("a", 1, 1, parity=1),
    BlockType("DtE", 2, 2, parity=0, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
    BlockType("GradE", 2, 3, parity=1, sym_groups=((1, 2),), tracefree_pairs=((1, 2),)),
    BlockType("Dt2E", 3, 2, parity=0, sym_groups=((0, 1),), tracefree_pairs=((0, 1),)),
)


def pairings(items: list[tuple[int, int]]) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    if not items:
        return ((),)
    first = items[0]
    out: list[tuple[tuple[tuple[int, int], tuple[int, int]], ...]] = []
    for idx in range(1, len(items)):
        second = items[idx]
        rest = items[1:idx] + items[idx + 1 :]
        for sub in pairings(rest):
            out.append(((first, second),) + sub)
    return tuple(out)


def slot_permutations(block: BlockType) -> tuple[tuple[int, ...], ...]:
    base = tuple(range(block.rank))
    perms = {base}
    for group in block.sym_groups:
        if len(group) == 2:
            perm = list(base)
            i, j = group
            perm[i], perm[j] = perm[j], perm[i]
            perms.add(tuple(perm))
    return tuple(sorted(perms))


def symmetry_maps(types: tuple[BlockType, ...]) -> tuple[dict[tuple[int, int], tuple[int, int]], ...]:
    groups: dict[str, list[int]] = {}
    for idx, block in enumerate(types):
        groups.setdefault(block.name, []).append(idx)

    identical_groups = [groups[name] for name in sorted(groups)]
    identical_perms = [tuple(permutations(indices)) for indices in identical_groups]
    internal_perms = [slot_permutations(block) for block in types]

    maps: list[dict[tuple[int, int], tuple[int, int]]] = []
    for ident_choice in product(*identical_perms):
        instance_map = {idx: idx for idx in range(len(types))}
        for indices, perm in zip(identical_groups, ident_choice):
            for old, new in zip(indices, perm):
                instance_map[old] = new
        for slot_choice in product(*internal_perms):
            mapping: dict[tuple[int, int], tuple[int, int]] = {}
            for old_inst, slot_perm in enumerate(slot_choice):
                new_inst = instance_map[old_inst]
                for src_slot, dst_slot in enumerate(slot_perm):
                    mapping[(old_inst, src_slot)] = (new_inst, dst_slot)
            maps.append(mapping)
    return tuple(maps)


def canonicalize(
    pairing: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
    maps: tuple[dict[tuple[int, int], tuple[int, int]], ...],
) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    candidates = []
    for mapping in maps:
        edges = []
        for left, right in pairing:
            lhs = mapping[left]
            rhs = mapping[right]
            if rhs < lhs:
                lhs, rhs = rhs, lhs
            edges.append((lhs, rhs))
        edges.sort()
        candidates.append(tuple(edges))
    return min(candidates)


def bad_trace(
    pairing: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
    types: tuple[BlockType, ...],
) -> bool:
    for (inst_a, slot_a), (inst_b, slot_b) in pairing:
        if inst_a != inst_b:
            continue
        for trace_pair in types[inst_a].tracefree_pairs:
            if {slot_a, slot_b} == set(trace_pair):
                return True
    return False


def weight_signature(types: tuple[BlockType, ...]) -> tuple[tuple[str, ...], int]:
    signature = tuple(sorted(block.name for block in types))
    weight = sum(block.weight for block in types)
    return signature, weight


def total_parity(types: tuple[BlockType, ...]) -> int:
    return sum(block.parity for block in types) % 2


def classify_and_label(
    signature: tuple[str, ...],
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...],
) -> tuple[str, str, str]:
    internal_pairs = {
        edge
        for edge in representative
        if edge[0][0] == edge[1][0]
    }
    if signature == ("a", "a"):
        return "a2", "Proven reducible", "lower-order EOM"
    if signature == ("E", "E"):
        return "E2", "Surviving candidate", "normal form"
    if signature == ("E", "a", "a"):
        return "aEa", "Proven reducible", "lower-order EOM"
    if signature == ("DtE", "E"):
        return "E_DtE", "Proven reducible", "total derivative"
    if signature == ("E", "E", "E"):
        return "E3", "Surviving candidate", "normal form"
    if signature == ("GradE", "a"):
        return "aDivE", "Proven reducible", "lower-order EOM"
    if signature == ("Dt2E", "E"):
        return "E_Dt2E", "Proven reducible", "total derivative"
    if signature == ("DtE", "DtE"):
        return "dotE2", "Surviving candidate", "normal form"
    if signature == ("DtE", "E", "E"):
        return "TrE2DtE", "Proven reducible", "total derivative"
    if signature == ("DtE", "a", "a"):
        return "aDtEa", "Proven reducible", "lower-order EOM"
    if signature == ("a", "a", "a", "a"):
        return "a4", "Proven reducible", "lower-order EOM"
    if signature == ("E", "E", "E", "E"):
        components = connected_components(representative)
        if len(components) == 2:
            return "E2^2", "Surviving candidate", "normal form"
        return "E4", "Proven reducible", "Cayley-Hamilton"
    if signature == ("E", "E", "a", "a"):
        components = connected_components(representative)
        if len(components) == 2:
            return "a2E2", "Proven reducible", "lower-order EOM"
        return "aE2a", "Proven reducible", "lower-order EOM"
    if signature == ("E", "GradE", "a"):
        if internal_pairs:
            return "aEGradE_1", "Proven reducible", "lower-order EOM"
        accel_edge = next(
            edge for edge in representative if edge[0][0] == 1 or edge[1][0] == 1
        )
        grad_endpoint = accel_edge[0] if accel_edge[0][0] == 2 else accel_edge[1]
        if grad_endpoint[1] == 0:
            return "aEGradE_3", "Proven reducible", "lower-order EOM"
        return "aEGradE_2", "Proven reducible", "lower-order EOM"
    if signature == ("GradE", "GradE"):
        if len(internal_pairs) == 2:
            return "divE2", "Surviving candidate", "survives current rules"
        if representative == (
            ((0, 0), (1, 0)),
            ((0, 1), (1, 1)),
            ((0, 2), (1, 2)),
        ):
            return "gradE2", "Surviving candidate", "normal form"
        return "mixedGradE2", "Surviving candidate", "survives current rules"
    raise ValueError(f"Unhandled signature {signature} with representative {representative}")


def connected_components(
    representative: tuple[tuple[tuple[int, int], tuple[int, int]], ...]
) -> tuple[frozenset[int], ...]:
    adjacency: dict[int, set[int]] = {}
    for left, right in representative:
        adjacency.setdefault(left[0], set()).add(right[0])
        adjacency.setdefault(right[0], set()).add(left[0])
    seen: set[int] = set()
    components: list[frozenset[int]] = []
    for node in sorted(adjacency):
        if node in seen:
            continue
        stack = [node]
        component = set()
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            component.add(current)
            stack.extend(adjacency[current] - seen)
        components.append(frozenset(component))
    return tuple(components)


def enumerate_contraction_classes(
    max_weight: int = MAX_WEIGHT,
    block_types: tuple[BlockType, ...] = BLOCK_TYPES,
    classifier: Callable[
        [tuple[str, ...], tuple[tuple[tuple[int, int], tuple[int, int]], ...]],
        tuple[str, str, str],
    ] = classify_and_label,
    require_even_parity: bool = True,
) -> tuple[ContractionClass, ...]:
    classes: list[ContractionClass] = []
    for count in range(1, 5):
        for combo in combinations_with_replacement(range(len(block_types)), count):
            types = tuple(block_types[index] for index in combo)
            signature, weight = weight_signature(types)
            total_rank = sum(block.rank for block in types)
            if weight > max_weight or total_rank % 2:
                continue
            if require_even_parity and total_parity(types):
                continue
            slots = [(inst, slot) for inst, block in enumerate(types) for slot in range(block.rank)]
            maps = symmetry_maps(types)
            unique: dict[
                tuple[tuple[tuple[int, int], tuple[int, int]], ...],
                tuple[tuple[tuple[int, int], tuple[int, int]], ...],
            ] = {}
            for pairing in pairings(slots):
                if bad_trace(pairing, types):
                    continue
                rep = canonicalize(pairing, maps)
                unique.setdefault(rep, rep)
            for rep in sorted(unique):
                label, classification, reduction_channel = classifier(signature, rep)
                classes.append(
                    ContractionClass(
                        signature=signature,
                        label=label,
                        classification=classification,
                        status="Proven" if classification == "Proven reducible" else "Conjectural",
                        weight=weight,
                        reduction_channel=reduction_channel,
                        representative=rep,
                    )
                )
    classes.sort(key=lambda item: (item.weight, item.signature, item.label))
    return tuple(classes)


def gradient_sector_classes() -> tuple[ContractionClass, ...]:
    return tuple(
        item for item in enumerate_contraction_classes() if item.signature == ("GradE", "GradE")
    )


def mixed_time_derivative_classes() -> tuple[ContractionClass, ...]:
    return tuple(
        item for item in enumerate_contraction_classes() if item.signature == ("DtE", "E", "E")
    )


def contraction_report(max_weight: int = MAX_WEIGHT) -> str:
    classes = enumerate_contraction_classes(max_weight=max_weight)
    lines = [
        f"Delta<=4 scalar contraction classes from the exact current primitive blocks",
        "",
        f"Total classes: {len(classes)}",
        "",
    ]
    current_signature: tuple[str, ...] | None = None
    for item in classes:
        if item.signature != current_signature:
            current_signature = item.signature
            lines.append(f"Signature {current_signature}:")
        lines.append(
            f"- {item.label} [weight {item.weight}; {item.classification}; {item.reduction_channel}]"
        )
    return "\n".join(lines)


def main() -> None:
    print(contraction_report())


if __name__ == "__main__":
    main()
