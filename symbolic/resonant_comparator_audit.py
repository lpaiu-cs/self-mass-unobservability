from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from second_order_mode_response import resonance_frequency_squared
from second_order_mode_response import symbols as mode_symbols
from second_order_mode_response import transfer_function


@dataclass(frozen=True)
class ResonanceBudget:
    polynomial_order: int
    projection_nuisance: int
    static_parameter_budget: int
    minimum_frequency_samples: int
    resonance_bracket_samples: int


@dataclass(frozen=True)
class ResonanceAuditRow:
    target: str
    status: str
    dynamic_law: str
    static_comparator: str
    verdict: str
    collapse_boundary: str


def _check_nonnegative(name: str, value: int) -> None:
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")


def line_shape_budget(polynomial_order: int, projection_nuisance: int = 0) -> ResonanceBudget:
    _check_nonnegative("polynomial_order", polynomial_order)
    _check_nonnegative("projection_nuisance", projection_nuisance)
    budget = polynomial_order + 1 + projection_nuisance
    minimum_samples = budget + 1
    return ResonanceBudget(
        polynomial_order=polynomial_order,
        projection_nuisance=projection_nuisance,
        static_parameter_budget=budget,
        minimum_frequency_samples=minimum_samples,
        resonance_bracket_samples=max(3, minimum_samples),
    )


def exact_polynomial_collapse_residual() -> sp.Expr:
    syms = mode_symbols()
    z = syms["z"]
    p0, p1, p2 = sp.symbols("p0 p1 p2")
    polynomial = p0 + p1 * z + p2 * z**2
    residual = sp.expand(
        polynomial
        * (
            syms["mu_chi"] * z**2
            + syms["gamma_chi"] * z
            + syms["omega_chi_sq"]
        )
        - syms["alpha"] * syms["omega_chi_sq"]
    )
    return residual


def dynamic_peak_condition() -> sp.Expr:
    syms = mode_symbols()
    return sp.simplify(2 * syms["mu_chi"] * syms["omega_chi_sq"] - syms["gamma_chi"] ** 2)


def phase_wrap_pole() -> sp.Expr:
    syms = mode_symbols()
    return sp.simplify(syms["omega_chi_sq"] / syms["mu_chi"])


def sample_design_rows() -> tuple[tuple[str, ResonanceBudget], ...]:
    return (
        ("linear comparator, no projection nuisance", line_shape_budget(1, 0)),
        ("linear comparator, acceleration projection nuisance", line_shape_budget(1, 1)),
        ("linear comparator, range projection nuisance", line_shape_budget(1, 2)),
        ("quadratic comparator, range projection nuisance", line_shape_budget(2, 2)),
    )


def resonance_audit_rows() -> tuple[ResonanceAuditRow, ...]:
    return (
        ResonanceAuditRow(
            target="resonant line shape",
            status="Counterexample candidate",
            dynamic_law="alpha omega_chi_sq/(mu_chi z^2+gamma_chi z+omega_chi_sq)",
            static_comparator="finite shared polynomial P_N(z) plus K projection nuisances",
            verdict="requires more than N+1+K complex frequency samples",
            collapse_boundary="mu_chi=0, adiabatic band, or arbitrary per-frequency nuisance",
        ),
        ResonanceAuditRow(
            target="phase wrapping",
            status="Counterexample candidate",
            dynamic_law="B/A=gamma_chi Omega/(omega_chi_sq-mu_chi Omega^2)",
            static_comparator="finite polynomial phase can mimic finite underbudget samples",
            verdict="requires samples bracketing omega_chi_sq/mu_chi and exceeding the budget",
            collapse_boundary="no resonance bracket, overdamped/no peak regime, or underbudget design",
        ),
    )


def audit_report() -> str:
    lines = [
        "Resonant comparator audit",
        "",
        "Second-order transfer:",
        f"- {sp.sstr(transfer_function())}",
        "Exact polynomial collapse residual for quadratic trial P_2:",
        f"- {sp.sstr(exact_polynomial_collapse_residual())}",
        "Resonance peak condition expression:",
        f"- {sp.sstr(dynamic_peak_condition())} > 0",
        "Phase-wrap pole location:",
        f"- Omega^2={sp.sstr(phase_wrap_pole())}",
        "",
        "Sample designs:",
    ]
    for label, budget in sample_design_rows():
        lines.append(
            f"- {label}: budget={budget.static_parameter_budget}, "
            f"minimum_samples={budget.minimum_frequency_samples}, "
            f"bracket_samples={budget.resonance_bracket_samples}"
        )
    lines.append("")
    lines.append("Audit rows:")
    for row in resonance_audit_rows():
        lines.append(
            f"- {row.status}: {row.target}; verdict={row.verdict}; "
            f"collapse={row.collapse_boundary}"
        )
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
