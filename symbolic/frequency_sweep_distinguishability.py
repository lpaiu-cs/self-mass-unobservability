from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class SweepSymbols:
    z_nodes: tuple[sp.Symbol, ...]
    z_star: sp.Symbol
    z: sp.Symbol
    beta: sp.Symbol
    tau: sp.Symbol


@dataclass(frozen=True)
class RealOddSymbols:
    u_nodes: tuple[sp.Symbol, ...]
    u_star: sp.Symbol
    beta: sp.Symbol
    tau: sp.Symbol


def symbols(order: int) -> SweepSymbols:
    if order < 0:
        raise ValueError("order must be non-negative")
    z_nodes = sp.symbols(f"z0:{order + 1}")
    return SweepSymbols(
        z_nodes=tuple(z_nodes),
        z_star=sp.Symbol("z_star"),
        z=sp.Symbol("z"),
        beta=sp.Symbol("beta", nonzero=True),
        tau=sp.Symbol("tau_chi", nonzero=True),
    )


def real_odd_symbols(order: int) -> RealOddSymbols:
    if order < 0:
        raise ValueError("order must be non-negative")
    node_count = max(0, (order + 1) // 2)
    u_nodes = sp.symbols(f"u0:{node_count}")
    if node_count == 1:
        u_nodes = (u_nodes,)
    return RealOddSymbols(
        u_nodes=tuple(u_nodes),
        u_star=sp.Symbol("u_star"),
        beta=sp.Symbol("beta", nonzero=True, real=True),
        tau=sp.Symbol("tau_chi", positive=True, real=True),
    )


def relaxation_transfer(z: sp.Expr, beta: sp.Expr, tau: sp.Expr) -> sp.Expr:
    return beta / (1 + tau * z)


def polynomial_sample_limit(order: int) -> dict[str, int]:
    if order < 0:
        raise ValueError("order must be non-negative")
    return {
        "degree": order,
        "exact_interpolation_limit": order + 1,
        "first_exact_obstruction_count": order + 2,
    }


def real_coefficient_sample_limit(order: int) -> dict[str, int]:
    if order < 0:
        raise ValueError("order must be non-negative")
    even_degree = order // 2
    odd_degree = (order - 1) // 2
    exact_limit = max(0, (order + 1) // 2)
    return {
        "degree": order,
        "even_channel_degree": even_degree,
        "odd_channel_degree": odd_degree,
        "positive_frequency_exact_interpolation_limit": exact_limit,
        "first_positive_frequency_obstruction_count": exact_limit + 1,
    }


def divided_difference(nodes: tuple[sp.Expr, ...], beta: sp.Expr, tau: sp.Expr) -> sp.Expr:
    if not nodes:
        raise ValueError("at least one node is required")
    order = len(nodes) - 1
    denominator = sp.prod(1 + tau * node for node in nodes)
    return sp.simplify(beta * (-tau) ** order / denominator)


def interpolation_residual(order: int) -> sp.Expr:
    """Residual at the (N+2)-th point after degree-N interpolation."""
    sym = symbols(order)
    all_nodes = sym.z_nodes + (sym.z_star,)
    product = sp.prod(sym.z_star - node for node in sym.z_nodes)
    return sp.factor(divided_difference(all_nodes, sym.beta, sym.tau) * product)


def real_odd_channel_residual(order: int) -> sp.Expr:
    """Residual for the imaginary channel with real derivative coefficients.

    For real d_n and positive frequencies, P_N(i Omega) splits into an even
    real polynomial and i Omega times an odd-channel polynomial in u=Omega^2.
    The relaxation odd channel is -beta*tau/(1+tau^2 u).
    """
    sym = real_odd_symbols(order)
    if not sym.u_nodes:
        return sp.factor(-sym.beta * sym.tau / (1 + sym.tau**2 * sym.u_star))

    all_nodes = sym.u_nodes + (sym.u_star,)
    gamma = -sym.beta * sym.tau
    denominator = sp.prod(1 + sym.tau**2 * node for node in all_nodes)
    divided = gamma * (-sym.tau**2) ** (len(all_nodes) - 1) / denominator
    product = sp.prod(sym.u_star - node for node in sym.u_nodes)
    return sp.factor(sp.simplify(divided * product))


def taylor_polynomial(order: int) -> sp.Expr:
    sym = symbols(order)
    return sp.expand(
        sym.beta * sum((-sym.tau * sym.z) ** n for n in range(order + 1))
    )


def taylor_residual(order: int) -> sp.Expr:
    sym = symbols(order)
    exact = relaxation_transfer(sym.z, sym.beta, sym.tau)
    return sp.factor(sp.simplify(exact - taylor_polynomial(order)))


def taylor_residual_closed_form(order: int) -> sp.Expr:
    sym = symbols(order)
    return sp.factor(
        sym.beta * (-sym.tau * sym.z) ** (order + 1) / (1 + sym.tau * sym.z)
    )


def low_frequency_error_bound(order: int) -> sp.Expr:
    rho, beta_abs = sp.symbols("rho abs_beta", nonnegative=True, real=True)
    return beta_abs * rho ** (order + 1)


def linear_interpolation_residual_example() -> dict[str, sp.Expr]:
    beta = sp.Integer(2)
    tau = sp.Integer(3)
    z0 = sp.I
    z1 = 2 * sp.I
    z_star = 3 * sp.I
    z = sp.Symbol("z")

    f0 = relaxation_transfer(z0, beta, tau)
    f1 = relaxation_transfer(z1, beta, tau)
    interpolant = f0 * (z - z1) / (z0 - z1) + f1 * (z - z0) / (z1 - z0)
    direct = sp.simplify(
        relaxation_transfer(z_star, beta, tau) - interpolant.subs(z, z_star)
    )

    sym = symbols(1)
    formula = interpolation_residual(1).subs(
        {
            sym.beta: beta,
            sym.tau: tau,
            sym.z_nodes[0]: z0,
            sym.z_nodes[1]: z1,
            sym.z_star: z_star,
        }
    )
    return {
        "direct_residual": sp.simplify(direct),
        "formula_residual": sp.simplify(formula),
    }


def summary(order: int = 3) -> dict[str, str | dict[str, int]]:
    return {
        "complex_coefficient_sample_boundary": polynomial_sample_limit(order),
        "real_coefficient_sample_boundary": real_coefficient_sample_limit(order),
        "interpolation_residual": sp.sstr(interpolation_residual(order)),
        "real_odd_channel_residual": sp.sstr(real_odd_channel_residual(order)),
        "taylor_residual": sp.sstr(taylor_residual(order)),
        "low_frequency_error_bound": sp.sstr(low_frequency_error_bound(order)),
    }


def main() -> None:
    data = summary()
    for key, value in data.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
