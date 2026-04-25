from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


I = sp.I


@dataclass(frozen=True)
class InterpolationObstruction:
    order: int
    nodes: tuple[sp.Symbol, ...]
    target: sp.Symbol
    residual: sp.Expr


@dataclass(frozen=True)
class ComparatorAuditRow:
    observable: str
    status: str
    dynamic_law: str
    static_comparator: str
    verdict: str
    collapse_boundary: str


def symbols() -> dict[str, sp.Symbol]:
    tau_chi, nu, n = sp.symbols("tau_chi nu n", nonzero=True)
    omega1, omega2 = sp.symbols("Omega1 Omega2", nonzero=True)
    c_y, beta, c_dyn = sp.symbols("c_Y beta C_dyn")
    z = sp.symbols("z")
    return {
        "tau_chi": tau_chi,
        "nu": nu,
        "n": n,
        "Omega1": omega1,
        "Omega2": omega2,
        "c_Y": c_y,
        "beta": beta,
        "C_dyn": c_dyn,
        "z": z,
    }


def relaxation_pole(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return 1 / (1 + syms["tau_chi"] * z)


def dynamic_sideband_transfer(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return syms["C_dyn"] * relaxation_pole(z)


def dynamic_sideband_cos_sin(nu: sp.Expr | None = None) -> tuple[sp.Expr, sp.Expr]:
    syms = symbols()
    if nu is None:
        nu = syms["nu"]
    x = syms["tau_chi"] * nu
    cos_component = syms["C_dyn"] / (1 + x**2)
    sin_component = syms["C_dyn"] * x / (1 + x**2)
    return sp.simplify(cos_component), sp.simplify(sin_component)


def shared_tau_phase_law(nu: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if nu is None:
        nu = syms["nu"]
    cos_component, sin_component = dynamic_sideband_cos_sin(nu)
    return sp.simplify(sin_component / (nu * cos_component))


def linear_dynamic_transfer(omega: sp.Expr) -> sp.Expr:
    syms = symbols()
    return syms["c_Y"] + syms["beta"] / (1 + I * omega * syms["tau_chi"])


def two_tone_sideband_ratio() -> sp.Expr:
    syms = symbols()
    omega1 = syms["Omega1"]
    omega2 = syms["Omega2"]
    tau_chi = syms["tau_chi"]
    numerator = syms["C_dyn"]
    denominator = (
        (1 + I * (omega1 + omega2) * tau_chi)
        * linear_dynamic_transfer(omega1)
        * linear_dynamic_transfer(omega2)
    )
    return sp.factor(numerator / denominator)


def orbital_three_n_ratio() -> sp.Expr:
    syms = symbols()
    n = syms["n"]
    tau_chi = syms["tau_chi"]
    numerator = syms["C_dyn"]
    denominator = (
        (1 + I * 3 * n * tau_chi)
        * linear_dynamic_transfer(n)
        * linear_dynamic_transfer(2 * n)
    )
    return sp.factor(numerator / denominator)


def static_single_line_mimic_coefficient(nu: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if nu is None:
        nu = syms["nu"]
    return dynamic_sideband_transfer(I * nu)


def sideband_pole_interpolation_obstruction(order: int) -> InterpolationObstruction:
    if order < 0:
        raise ValueError("order must be nonnegative")
    syms = symbols()
    tau_chi = syms["tau_chi"]
    c_dyn = syms["C_dyn"]
    target = sp.Symbol("z_star")
    nodes = sp.symbols(f"z0:{order + 1}")
    node_product = sp.prod(1 + tau_chi * node for node in nodes)
    target_product = sp.prod(target - node for node in nodes)
    residual = sp.simplify(
        c_dyn
        * (-tau_chi) ** (order + 1)
        * target_product
        / ((1 + tau_chi * target) * node_product)
    )
    return InterpolationObstruction(
        order=order,
        nodes=tuple(nodes),
        target=target,
        residual=residual,
    )


def audit_rows() -> tuple[ComparatorAuditRow, ...]:
    return (
        ComparatorAuditRow(
            observable="single generated sideband",
            status="Proven",
            dynamic_law="C_dyn/(1 + tau_chi z)",
            static_comparator="one free nonlinear coefficient at the sampled sideband",
            verdict="not unique to dynamic chi",
            collapse_boundary="a static nonlinear coefficient can match one complex sideband sample",
        ),
        ComparatorAuditRow(
            observable="shared sideband pole over many generated lines",
            status="Counterexample candidate",
            dynamic_law="C_dyn/(1 + tau_chi z) at every generated sideband z",
            static_comparator="shared finite polynomial P_N(z)",
            verdict="distinguishable after the finite interpolation budget is exceeded",
            collapse_boundary="complex P_N absorbs at most N+1 distinct samples; the next distinct sample has the pole residual",
        ),
        ComparatorAuditRow(
            observable="two-tone ratio O(Omega1+Omega2)/(O(Omega1)O(Omega2))",
            status="Counterexample candidate",
            dynamic_law="C_dyn/[(1+i(Omega1+Omega2)tau_chi)G(Omega1)G(Omega2)]",
            static_comparator="local nonlinear derivative terms in F and its derivatives",
            verdict="useful only if the same tau_chi also fits the linear lines",
            collapse_boundary="arbitrary per-sideband nonlinear nuisance destroys the ratio test",
        ),
        ComparatorAuditRow(
            observable="orbital n,2n -> 3n ratio",
            status="Counterexample candidate",
            dynamic_law="C_dyn/[(1+3 i n tau_chi)G(n)G(2n)]",
            static_comparator="static nonlinear F^2 and derivative products can create 3n",
            verdict="sideband existence is insufficient; shared pole phase is the target",
            collapse_boundary="if static nonlinear coefficients are independent at 3n, the sideband alone is degenerate",
        ),
    )


def audit_report() -> str:
    obstruction = sideband_pole_interpolation_obstruction(order=1)
    lines = [
        "Nonlinear comparator audit",
        "",
        "Dynamic sideband transfer:",
        f"- {sp.sstr(dynamic_sideband_transfer())}",
        "Shared tau phase law:",
        f"- sin_component/(nu*cos_component) = {sp.sstr(shared_tau_phase_law())}",
        "Two-tone shared-pole ratio:",
        f"- {sp.sstr(two_tone_sideband_ratio())}",
        "Orbital 3n shared-pole ratio:",
        f"- {sp.sstr(orbital_three_n_ratio())}",
        "Degree-1 interpolation residual at a third sample:",
        f"- {sp.sstr(obstruction.residual)}",
        "",
        "Audit rows:",
    ]
    for row in audit_rows():
        lines.append(
            f"- {row.status}: {row.observable}; verdict={row.verdict}; "
            f"collapse={row.collapse_boundary}"
        )
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
