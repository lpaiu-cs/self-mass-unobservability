from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class TwoFrequencyResponse:
    t: sp.Symbol
    omega1: sp.Symbol
    omega2: sp.Symbol
    tau: sp.Symbol
    alpha: sp.Symbol
    f1: sp.Symbol
    f2: sp.Symbol
    drive: sp.Expr
    chi: sp.Expr
    components: dict[sp.Symbol, tuple[sp.Expr, sp.Expr]]


def _component(alpha: sp.Expr, f: sp.Expr, omega: sp.Expr, tau: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    x = omega * tau
    in_phase = alpha * f / (1 + x**2)
    quadrature = alpha * f * x / (1 + x**2)
    return sp.simplify(in_phase), sp.simplify(quadrature)


def two_frequency_response() -> TwoFrequencyResponse:
    t = sp.Symbol("t", real=True)
    omega1 = sp.Symbol("Omega1", positive=True, real=True)
    omega2 = sp.Symbol("Omega2", positive=True, real=True)
    tau = sp.Symbol("tau_chi", positive=True, real=True)
    alpha = sp.Symbol("alpha", positive=True, real=True)
    f1 = sp.Symbol("F1", positive=True, real=True)
    f2 = sp.Symbol("F2", positive=True, real=True)

    a1, b1 = _component(alpha, f1, omega1, tau)
    a2, b2 = _component(alpha, f2, omega2, tau)
    drive = f1 * sp.cos(omega1 * t) + f2 * sp.cos(omega2 * t)
    chi = (
        a1 * sp.cos(omega1 * t)
        + b1 * sp.sin(omega1 * t)
        + a2 * sp.cos(omega2 * t)
        + b2 * sp.sin(omega2 * t)
    )

    return TwoFrequencyResponse(
        t=t,
        omega1=omega1,
        omega2=omega2,
        tau=tau,
        alpha=alpha,
        f1=f1,
        f2=f2,
        drive=drive,
        chi=sp.simplify(chi),
        components={omega1: (a1, b1), omega2: (a2, b2)},
    )


def relaxation_residual(response: TwoFrequencyResponse | None = None) -> sp.Expr:
    response = response or two_frequency_response()
    lhs = response.tau * sp.diff(response.chi, response.t) + response.chi
    rhs = response.alpha * response.drive
    return sp.trigsimp(sp.simplify(lhs - rhs))


def linear_sideband_residuals() -> dict[str, int]:
    """Return sideband amplitudes for the linear model.

    Linearity means there are no Omega1 + Omega2 or |Omega1 - Omega2| terms.
    """
    return {"sum": 0, "difference": 0}


def main() -> None:
    response = two_frequency_response()
    print(f"drive: {sp.sstr(response.drive)}")
    print(f"chi: {sp.sstr(response.chi)}")
    print(f"sidebands: {linear_sideband_residuals()}")
    print(f"residual: {sp.sstr(relaxation_residual(response))}")


if __name__ == "__main__":
    main()
