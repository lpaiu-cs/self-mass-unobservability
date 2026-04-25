from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from second_order_mode_response import symbols as mode_symbols
from second_order_mode_response import transfer_function


@dataclass(frozen=True)
class ProjectionAuditRow:
    channel: str
    status: str
    observed_transfer: str
    diagnostic: str
    verdict: str
    collapse_boundary: str


def symbols() -> dict[str, sp.Symbol]:
    syms = mode_symbols()
    gamma, kappa_sq = sp.symbols("Gamma kappa_sq")
    syms.update({"Gamma": gamma, "kappa_sq": kappa_sq})
    return syms


def acceleration_observed_transfer(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return syms["Gamma"] * transfer_function(z)


def range_observed_transfer(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return syms["Gamma"] * transfer_function(z) / (syms["kappa_sq"] + z**2)


def internal_denominator(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return syms["mu_chi"] * z**2 + syms["gamma_chi"] * z + syms["omega_chi_sq"]


def acceleration_diagnostic_residual(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    residual = (
        acceleration_observed_transfer(z) * internal_denominator(z)
        - syms["Gamma"] * syms["alpha"] * syms["omega_chi_sq"]
    )
    return sp.cancel(residual)


def range_diagnostic_residual(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    residual = (
        range_observed_transfer(z)
        * (syms["kappa_sq"] + z**2)
        * internal_denominator(z)
        - syms["Gamma"] * syms["alpha"] * syms["omega_chi_sq"]
    )
    return sp.cancel(residual)


def deproject_acceleration(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return sp.cancel(acceleration_observed_transfer(z) / syms["Gamma"])


def deproject_range(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    return sp.cancel(range_observed_transfer(z) * (syms["kappa_sq"] + z**2) / syms["Gamma"])


def projection_nuisance_budget(channel: str) -> int:
    if channel == "acceleration":
        return 1
    if channel == "range":
        return 2
    raise ValueError(f"unknown channel: {channel}")


def projection_audit_rows() -> tuple[ProjectionAuditRow, ...]:
    return (
        ProjectionAuditRow(
            channel="acceleration",
            status="Proven",
            observed_transfer="Gamma H2(z)",
            diagnostic="O_a(z)D_chi(z)-Gamma alpha omega_chi_sq = 0",
            verdict="internal resonance survives finite shared Gamma",
            collapse_boundary="Gamma=0 or arbitrary per-frequency acceleration nuisance",
        ),
        ProjectionAuditRow(
            channel="range",
            status="Proven",
            observed_transfer="Gamma H2(z)/(kappa_sq+z^2)",
            diagnostic="O_R(z)(kappa_sq+z^2)D_chi(z)-Gamma alpha omega_chi_sq = 0",
            verdict="internal resonance survives finite shared Gamma,kappa_sq after deprojection",
            collapse_boundary="Gamma=0, projection pole sampling, pole-cancelling projection, or arbitrary per-frequency range nuisance",
        ),
    )


def audit_report() -> str:
    lines = [
        "Second-order projection audit",
        "",
        "Acceleration channel:",
        f"- observed={sp.sstr(acceleration_observed_transfer())}",
        f"- diagnostic residual={sp.sstr(acceleration_diagnostic_residual())}",
        f"- deprojected={sp.sstr(deproject_acceleration())}",
        "Range channel:",
        f"- observed={sp.sstr(range_observed_transfer())}",
        f"- diagnostic residual={sp.sstr(range_diagnostic_residual())}",
        f"- deprojected={sp.sstr(deproject_range())}",
        "Projection nuisance budgets:",
        f"- acceleration={projection_nuisance_budget('acceleration')}",
        f"- range={projection_nuisance_budget('range')}",
        "",
        "Audit rows:",
    ]
    for row in projection_audit_rows():
        lines.append(
            f"- {row.status}: {row.channel}; verdict={row.verdict}; "
            f"collapse={row.collapse_boundary}"
        )
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
