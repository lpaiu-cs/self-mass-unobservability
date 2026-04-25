from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


I = sp.I


@dataclass(frozen=True)
class SecondOrderSummary:
    transfer: sp.Expr
    low_frequency_series: sp.Expr
    first_order_limit: sp.Expr
    resonance_frequency_squared: sp.Expr


def symbols() -> dict[str, sp.Symbol]:
    mu_chi, gamma_chi, omega_chi_sq = sp.symbols(
        "mu_chi gamma_chi omega_chi_sq",
        positive=True,
    )
    alpha, z, omega = sp.symbols("alpha z Omega")
    return {
        "mu_chi": mu_chi,
        "gamma_chi": gamma_chi,
        "omega_chi_sq": omega_chi_sq,
        "alpha": alpha,
        "z": z,
        "Omega": omega,
    }


def transfer_function(z: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if z is None:
        z = syms["z"]
    denominator = (
        syms["mu_chi"] * z**2
        + syms["gamma_chi"] * z
        + syms["omega_chi_sq"]
    )
    return syms["alpha"] * syms["omega_chi_sq"] / denominator


def frequency_transfer(omega: sp.Expr | None = None) -> sp.Expr:
    syms = symbols()
    if omega is None:
        omega = syms["Omega"]
    return transfer_function(I * omega)


def cos_sin_components(omega: sp.Expr | None = None) -> tuple[sp.Expr, sp.Expr]:
    syms = symbols()
    if omega is None:
        omega = syms["Omega"]
    real_den = syms["omega_chi_sq"] - syms["mu_chi"] * omega**2
    imag_den = syms["gamma_chi"] * omega
    norm = real_den**2 + imag_den**2
    cos_component = syms["alpha"] * syms["omega_chi_sq"] * real_den / norm
    sin_component = syms["alpha"] * syms["omega_chi_sq"] * imag_den / norm
    return sp.simplify(cos_component), sp.simplify(sin_component)


def phase_tangent(omega: sp.Expr | None = None) -> sp.Expr:
    if omega is None:
        omega = symbols()["Omega"]
    cos_component, sin_component = cos_sin_components(omega)
    return sp.simplify(sin_component / cos_component)


def low_frequency_series(order: int = 4) -> sp.Expr:
    syms = symbols()
    return sp.series(transfer_function(), syms["z"], 0, order).removeO()


def first_order_limit() -> sp.Expr:
    syms = symbols()
    return sp.simplify(transfer_function().subs(syms["mu_chi"], 0))


def resonance_frequency_squared() -> sp.Expr:
    syms = symbols()
    mu_chi = syms["mu_chi"]
    gamma_chi = syms["gamma_chi"]
    omega_chi_sq = syms["omega_chi_sq"]
    return sp.simplify(omega_chi_sq / mu_chi - gamma_chi**2 / (2 * mu_chi**2))


def second_order_summary() -> SecondOrderSummary:
    return SecondOrderSummary(
        transfer=transfer_function(),
        low_frequency_series=low_frequency_series(order=4),
        first_order_limit=first_order_limit(),
        resonance_frequency_squared=resonance_frequency_squared(),
    )


def audit_report() -> str:
    syms = symbols()
    cos_component, sin_component = cos_sin_components()
    summary = second_order_summary()
    lines = [
        "Second-order internal mode response",
        "",
        "Transfer:",
        f"- {sp.sstr(summary.transfer)}",
        "Frequency response:",
        f"- {sp.sstr(frequency_transfer())}",
        "Cos/sin components:",
        f"- cos={sp.sstr(cos_component)}",
        f"- sin={sp.sstr(sin_component)}",
        "Phase tangent:",
        f"- {sp.sstr(phase_tangent())}",
        "Low-frequency derivative expansion:",
        f"- {sp.sstr(summary.low_frequency_series)}",
        "First-order relaxation limit mu_chi=0:",
        f"- {sp.sstr(summary.first_order_limit)}",
        "Resonance condition:",
        f"- Omega_peak^2={sp.sstr(summary.resonance_frequency_squared)}",
        f"- exists when gamma_chi^2 < {sp.sstr(2 * syms['mu_chi'] * syms['omega_chi_sq'])}",
    ]
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
