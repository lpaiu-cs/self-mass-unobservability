from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class SensitivityJet:
    invariants: tuple[sp.Expr, ...]
    linear: sp.Matrix
    quadratic: sp.Matrix
    mass_scale: sp.Symbol

    def mass_shift(self) -> sp.Expr:
        y = sp.Matrix(self.invariants)
        linear_term = (self.linear.T * y)[0]
        quadratic_term = sp.Rational(1, 2) * (y.T * self.quadratic * y)[0]
        return sp.expand(self.mass_scale * (1 + linear_term + quadratic_term))

    def coordinates(self) -> tuple[sp.Symbol, ...]:
        coords: list[sp.Symbol] = [self.mass_scale]
        coords.extend(self.linear)
        for i in range(self.quadratic.rows):
            for j in range(i, self.quadratic.cols):
                coords.append(self.quadratic[i, j])
        return tuple(coords)


def make_quadratic_jet(
    prefix: str,
    invariants: tuple[sp.Expr, ...],
    mass_scale: sp.Symbol | None = None,
) -> SensitivityJet:
    if not invariants:
        raise ValueError("At least one invariant is required.")

    size = len(invariants)
    mass = mass_scale if mass_scale is not None else sp.Symbol(f"{prefix}0")
    linear = sp.Matrix(size, 1, lambda i, _j: sp.Symbol(f"{prefix}_{i + 1}"))
    quadratic = sp.Matrix(
        size,
        size,
        lambda i, j: sp.Symbol(f"{prefix}_{min(i, j) + 1}{max(i, j) + 1}"),
    )
    return SensitivityJet(tuple(invariants), linear, quadratic, mass)


def demo_expression() -> sp.Expr:
    y1, y2 = sp.symbols("Y1 Y2")
    jet = make_quadratic_jet("s", (y1, y2), mass_scale=sp.Symbol("m0"))
    return jet.mass_shift()


def main() -> None:
    expr = demo_expression()
    print("Quadratic sensitivity jet:")
    print(sp.sstr(expr))


if __name__ == "__main__":
    main()
