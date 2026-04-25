from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


I = sp.I


@dataclass(frozen=True)
class RatioAuditRow:
    target: str
    status: str
    dynamic_statement: str
    static_comparator: str
    verdict: str
    collapse_boundary: str


def symbols() -> dict[str, sp.Symbol]:
    tau_chi, u, v = sp.symbols("tau_chi u v", nonzero=True)
    c_y, beta, c_side = sp.symbols("c_Y beta C_side")
    return {
        "tau_chi": tau_chi,
        "u": u,
        "v": v,
        "c_Y": c_y,
        "beta": beta,
        "C_side": c_side,
    }


def h_pole(z: sp.Expr) -> sp.Expr:
    syms = symbols()
    return 1 / (1 + syms["tau_chi"] * z)


def g_linear(z: sp.Expr) -> sp.Expr:
    syms = symbols()
    return syms["c_Y"] + syms["beta"] * h_pole(z)


def shared_tau_ratio(u: sp.Expr | None = None, v: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if u is None:
        u = syms["u"]
    if v is None:
        v = syms["v"]
    return syms["C_side"] * h_pole(u + v) / (g_linear(u) * g_linear(v))


def shared_tau_diagnostic_residual(
    ratio: sp.Expr | None = None,
    u: sp.Expr | None = None,
    v: sp.Expr | None = None,
) -> sp.Expr:
    syms = symbols()
    if u is None:
        u = syms["u"]
    if v is None:
        v = syms["v"]
    if ratio is None:
        ratio = shared_tau_ratio(u, v)
    residual = ratio * g_linear(u) * g_linear(v) * (1 + syms["tau_chi"] * (u + v))
    return sp.cancel(residual - syms["C_side"])


def tau_zero_limit() -> sp.Expr:
    syms = symbols()
    return sp.simplify(shared_tau_ratio().subs(syms["tau_chi"], 0))


def no_sideband_limit() -> sp.Expr:
    syms = symbols()
    return sp.simplify(shared_tau_ratio().subs(syms["C_side"], 0))


def static_linear_polynomial(order: int, z: sp.Expr) -> sp.Expr:
    if order < 0:
        raise ValueError("order must be nonnegative")
    coeffs = sp.symbols(f"a0:{order + 1}")
    return sum(coeff * z**power for power, coeff in enumerate(coeffs))


def bivariate_monomial_count(total_degree: int) -> int:
    if total_degree < 0:
        raise ValueError("total_degree must be nonnegative")
    return (total_degree + 1) * (total_degree + 2) // 2


def static_sideband_polynomial(total_degree: int, u: sp.Expr | None = None, v: sp.Expr | None = None) -> sp.Expr:
    if total_degree < 0:
        raise ValueError("total_degree must be nonnegative")
    syms = symbols()
    if u is None:
        u = syms["u"]
    if v is None:
        v = syms["v"]
    terms: list[sp.Expr] = []
    coeff_index = 0
    for degree_u in range(total_degree + 1):
        for degree_v in range(total_degree + 1 - degree_u):
            coeff = sp.Symbol(f"b{coeff_index}")
            terms.append(coeff * u**degree_u * v**degree_v)
            coeff_index += 1
    return sum(terms)


def static_ratio_template(linear_order: int, sideband_degree: int) -> sp.Expr:
    syms = symbols()
    u = syms["u"]
    v = syms["v"]
    p_u = static_linear_polynomial(linear_order, u)
    p_v = static_linear_polynomial(linear_order, v)
    q_uv = static_sideband_polynomial(sideband_degree, u, v)
    return q_uv / (p_u * p_v)


def interpolation_budget(linear_order: int, sideband_degree: int) -> dict[str, int]:
    return {
        "linear_complex_samples": linear_order + 1,
        "linear_obstruction_sample": linear_order + 2,
        "sideband_pair_samples": bivariate_monomial_count(sideband_degree),
        "sideband_obstruction_pair": bivariate_monomial_count(sideband_degree) + 1,
    }


def ratio_audit_rows() -> tuple[RatioAuditRow, ...]:
    return (
        RatioAuditRow(
            target="shared-tau diagnostic residual",
            status="Proven",
            dynamic_statement="R(u,v)G(u)G(v)[1+tau_chi(u+v)]-C_side = 0",
            static_comparator="finite P_N(u), P_N(v), and Q_M(u,v)",
            verdict="the diagnostic is exact for the dynamic model",
            collapse_boundary="C_side=0 or tau_chi=0 removes the dynamic sideband pole content",
        ),
        RatioAuditRow(
            target="finite static nonlinear interpolation",
            status="Proven",
            dynamic_statement="R(u,v)=C_side H(u+v)/(G(u)G(v))",
            static_comparator="Q_M(u,v)/(P_N(u)P_N(v))",
            verdict="finite samples can be absorbed up to the interpolation budget",
            collapse_boundary="linear samples <= N+1 and sideband pairs <= (M+1)(M+2)/2",
        ),
        RatioAuditRow(
            target="overbudget shared-law test",
            status="Counterexample candidate",
            dynamic_statement="one tau_chi must fit linear poles and generated-line pole together",
            static_comparator="finite shared local derivative coefficients",
            verdict="distinguishable only after the linear or sideband interpolation budget is exceeded",
            collapse_boundary="arbitrary per-frequency or per-sideband nuisance destroys the test",
        ),
    )


def audit_report() -> str:
    budget = interpolation_budget(linear_order=1, sideband_degree=2)
    lines = [
        "Shared-tau ratio audit",
        "",
        "Dynamic ratio:",
        f"- {sp.sstr(shared_tau_ratio())}",
        "Diagnostic residual:",
        f"- {sp.sstr(shared_tau_diagnostic_residual())}",
        "tau_chi=0 limit:",
        f"- {sp.sstr(tau_zero_limit())}",
        "C_side=0 limit:",
        f"- {sp.sstr(no_sideband_limit())}",
        "Example static interpolation budget, N=1 and M=2:",
        f"- {budget}",
        "",
        "Audit rows:",
    ]
    for row in ratio_audit_rows():
        lines.append(
            f"- {row.status}: {row.target}; verdict={row.verdict}; "
            f"collapse={row.collapse_boundary}"
        )
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
