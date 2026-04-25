from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from resonant_comparator_audit import ResonanceBudget, line_shape_budget


@dataclass(frozen=True)
class InterpolationObstruction:
    order: int
    nodes: tuple[sp.Symbol, ...]
    target: sp.Symbol
    residual: sp.Expr


@dataclass(frozen=True)
class NonlinearSecondOrderRow:
    target: str
    status: str
    dynamic_law: str
    static_comparator: str
    verdict: str
    collapse_boundary: str


def symbols() -> dict[str, sp.Symbol]:
    mu_chi, gamma_chi, omega_chi_sq = sp.symbols(
        "mu_chi gamma_chi omega_chi_sq",
        positive=True,
    )
    alpha, beta_f2 = sp.symbols("alpha beta_F2")
    c_y, c_chi = sp.symbols("c_Y c_chi")
    lambda_f_chi, lambda_chi2 = sp.symbols("lambda_Fchi lambda_chi2")
    z, u, v = sp.symbols("z u v")
    return {
        "mu_chi": mu_chi,
        "gamma_chi": gamma_chi,
        "omega_chi_sq": omega_chi_sq,
        "alpha": alpha,
        "beta_F2": beta_f2,
        "c_Y": c_y,
        "c_chi": c_chi,
        "lambda_Fchi": lambda_f_chi,
        "lambda_chi2": lambda_chi2,
        "z": z,
        "u": u,
        "v": v,
    }


def denominator(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return (
        syms["mu_chi"] * z**2
        + syms["gamma_chi"] * z
        + syms["omega_chi_sq"]
    )


def internal_linear_transfer(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return syms["alpha"] * syms["omega_chi_sq"] / denominator(z)


def observable_linear_transfer(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return syms["c_Y"] + syms["c_chi"] * internal_linear_transfer(z)


def nonlinear_drive_sideband_transfer(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return syms["c_chi"] * syms["beta_F2"] / denominator(z)


def nonlinear_readout_sideband_transfer(
    u: sp.Expr | None = None,
    v: sp.Expr | None = None,
) -> sp.Expr:
    syms = symbols()
    if u is None:
        u = syms["u"]
    if v is None:
        v = syms["v"]
    h_u = internal_linear_transfer(u)
    h_v = internal_linear_transfer(v)
    return sp.factor(
        syms["lambda_Fchi"] * (h_u + h_v)
        + syms["lambda_chi2"] * h_u * h_v
    )


def shared_denominator_residual() -> sp.Expr:
    syms = symbols()
    z0, z1 = sp.symbols("z0 z1")
    linear0 = observable_linear_transfer(z0)
    linear1 = observable_linear_transfer(z1)
    side0 = nonlinear_drive_sideband_transfer(z0)
    side1 = nonlinear_drive_sideband_transfer(z1)
    residual = (linear0 - syms["c_Y"]) * side1 - (linear1 - syms["c_Y"]) * side0
    return sp.factor(residual)


def drive_sideband_interpolation_obstruction(order: int = 1) -> InterpolationObstruction:
    if order < 0:
        raise ValueError("order must be nonnegative")
    syms = symbols()
    target = sp.Symbol("z_star")
    nodes = sp.symbols(f"z0:{order + 1}")
    z = syms["z"]
    target_function = nonlinear_drive_sideband_transfer(z)
    polynomial = 0
    for node in nodes:
        basis = sp.prod(
            (z - other) / (node - other)
            for other in nodes
            if other != node
        )
        polynomial += target_function.subs(z, node) * basis
    residual = sp.factor(
        target_function.subs(z, target) - polynomial.subs(z, target)
    )
    return InterpolationObstruction(
        order=order,
        nodes=tuple(nodes),
        target=target,
        residual=residual,
    )


def nonlinear_sideband_budget(
    polynomial_order: int,
    projection_nuisance: int = 0,
) -> ResonanceBudget:
    return line_shape_budget(polynomial_order, projection_nuisance)


def audit_rows() -> tuple[NonlinearSecondOrderRow, ...]:
    return (
        NonlinearSecondOrderRow(
            target="nonlinear-drive generated sideband",
            status="Proven",
            dynamic_law="C_beta/(mu_chi z^2+gamma_chi z+omega_chi_sq)",
            static_comparator="one static nonlinear coefficient at one generated line",
            verdict="sideband existence is not unique",
            collapse_boundary="one generated complex sample is exactly mimicked by one free static nonlinear coefficient",
        ),
        NonlinearSecondOrderRow(
            target="shared second-order denominator",
            status="Counterexample candidate",
            dynamic_law="G_L(z)=c_Y+A/D_2(z), S_beta(z)=B/D_2(z)",
            static_comparator="finite shared polynomial line and sideband comparators",
            verdict="candidate only after both finite budgets are exceeded or a calibrated cross-ratio is used",
            collapse_boundary="arbitrary per-frequency nonlinear or projection nuisance kills the shared-denominator test",
        ),
        NonlinearSecondOrderRow(
            target="nonlinear readout sideband",
            status="Counterexample candidate",
            dynamic_law="lambda_Fchi[H_2(u)+H_2(v)]+lambda_chi2 H_2(u)H_2(v)",
            static_comparator="finite local products F d^aF d^bF with shared coefficients",
            verdict="stronger than line creation only if lambda nuisances are counted and the shared H_2 law is overbudget",
            collapse_boundary="free readout coefficients at every sideband absorb finite underbudget data",
        ),
        NonlinearSecondOrderRow(
            target="resonance-assisted sideband sweep",
            status="Counterexample candidate",
            dynamic_law="generated lines inherit the same quadratic denominator or products of it",
            static_comparator="static nonlinear polynomial of generated frequency plus projection nuisance",
            verdict="requires M+2+K complex generated samples and resonance bracketing",
            collapse_boundary="underbudget samples, no resonance bracket, beta_F2=0, c_chi=0, or mu_chi=0",
        ),
    )


def audit_report() -> str:
    obstruction = drive_sideband_interpolation_obstruction(order=1)
    budget = nonlinear_sideband_budget(polynomial_order=1, projection_nuisance=1)
    lines = [
        "Nonlinear second-order internal mode",
        "",
        "Linear observable transfer:",
        f"- {sp.sstr(observable_linear_transfer())}",
        "Nonlinear-drive sideband transfer:",
        f"- {sp.sstr(nonlinear_drive_sideband_transfer())}",
        "Nonlinear-readout sideband transfer:",
        f"- {sp.sstr(nonlinear_readout_sideband_transfer())}",
        "Shared-denominator residual:",
        f"- {sp.sstr(shared_denominator_residual())}",
        "Degree-1 static sideband interpolation residual:",
        f"- {sp.sstr(obstruction.residual)}",
        "Example generated-sideband budget:",
        f"- M={budget.polynomial_order}, K={budget.projection_nuisance}, "
        f"minimum_samples={budget.minimum_frequency_samples}",
        "",
        "Audit rows:",
    ]
    for row in audit_rows():
        lines.append(
            f"- {row.status}: {row.target}; verdict={row.verdict}; "
            f"collapse={row.collapse_boundary}"
        )
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
