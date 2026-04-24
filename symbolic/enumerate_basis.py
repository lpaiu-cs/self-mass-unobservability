from __future__ import annotations

from dataclasses import dataclass
from itertools import product


DEFAULT_MAX_WEIGHT = 4


@dataclass(frozen=True)
class Primitive:
    name: str
    weight: int
    description: str


@dataclass(frozen=True)
class Monomial:
    label: str
    weight: int
    exponents: tuple[int, ...]


def minimal_sector_primitives() -> tuple[Primitive, ...]:
    return (
        Primitive("E2", 2, "E_ij E^ij"),
        Primitive("E3", 3, "E_i^j E_j^k E_k^i"),
        Primitive("dotE2", 4, "(D_tau E_ij)(D_tau E^ij)"),
        Primitive("gradE2", 4, "(nabla_k E_ij)(nabla^k E^ij)"),
    )


def format_monomial(
    primitives: tuple[Primitive, ...], exponents: tuple[int, ...]
) -> str:
    factors: list[str] = []
    for primitive, exponent in zip(primitives, exponents):
        if exponent == 0:
            continue
        if exponent == 1:
            factors.append(primitive.name)
        else:
            factors.append(f"{primitive.name}^{exponent}")
    return " * ".join(factors) if factors else "1"


def enumerate_monomials(
    primitives: tuple[Primitive, ...], max_weight: int = DEFAULT_MAX_WEIGHT
) -> tuple[Monomial, ...]:
    if max_weight < 0:
        raise ValueError("max_weight must be non-negative")
    if any(primitive.weight <= 0 for primitive in primitives):
        raise ValueError("primitive weights must be positive")

    bounds = [max_weight // primitive.weight for primitive in primitives]
    monomials: list[Monomial] = []
    for exponents in product(*(range(bound + 1) for bound in bounds)):
        weight = sum(
            exponent * primitive.weight
            for exponent, primitive in zip(exponents, primitives)
        )
        if weight > max_weight:
            continue
        monomials.append(
            Monomial(
                label=format_monomial(primitives, exponents),
                weight=weight,
                exponents=tuple(exponents),
            )
        )

    monomials.sort(
        key=lambda item: (item.weight, 0 if item.label == "1" else 1, item.label)
    )
    return tuple(monomials)


def enumerate_minimal_scalar_monomials(
    max_weight: int = DEFAULT_MAX_WEIGHT,
) -> tuple[Monomial, ...]:
    return enumerate_monomials(minimal_sector_primitives(), max_weight=max_weight)


def basis_report(max_weight: int = DEFAULT_MAX_WEIGHT) -> str:
    primitives = minimal_sector_primitives()
    monomials = enumerate_minimal_scalar_monomials(max_weight=max_weight)

    lines = [
        f"Minimal-sector candidate scalar monomials through Delta <= {max_weight}",
        "",
        "Primitive generators:",
    ]
    for primitive in primitives:
        lines.append(
            f"- {primitive.name} [weight {primitive.weight}]: {primitive.description}"
        )

    lines.extend(["", "Enumerated monomials:"])
    for monomial in monomials:
        lines.append(f"- weight {monomial.weight}: {monomial.label}")
    return "\n".join(lines)


def main() -> None:
    print(basis_report())


if __name__ == "__main__":
    main()
