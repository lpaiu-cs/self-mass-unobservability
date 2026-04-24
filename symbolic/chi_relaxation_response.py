from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class RelaxationResponse:
    t: sp.Symbol
    omega: sp.Symbol
    tau: sp.Symbol
    alpha: sp.Symbol
    f0: sp.Symbol
    drive: sp.Expr
    chi: sp.Expr
    in_phase: sp.Expr
    quadrature: sp.Expr
    amplitude: sp.Expr
    phase_lag: sp.Expr


def transfer_function(
    *,
    alpha: sp.Expr | None = None,
    omega: sp.Expr | None = None,
    tau: sp.Expr | None = None,
) -> sp.Expr:
    """Return H(omega) for tau * dot(chi) + chi = alpha * F."""
    alpha = alpha if alpha is not None else sp.Symbol("alpha", positive=True, real=True)
    omega = omega if omega is not None else sp.Symbol("Omega", positive=True, real=True)
    tau = tau if tau is not None else sp.Symbol("tau_chi", positive=True, real=True)
    return alpha / (1 + sp.I * omega * tau)


def steady_state_response() -> RelaxationResponse:
    t = sp.Symbol("t", real=True)
    omega = sp.Symbol("Omega", positive=True, real=True)
    tau = sp.Symbol("tau_chi", positive=True, real=True)
    alpha = sp.Symbol("alpha", positive=True, real=True)
    f0 = sp.Symbol("F0", positive=True, real=True)
    x = omega * tau

    drive = f0 * sp.cos(omega * t)
    in_phase = alpha * f0 / (1 + x**2)
    quadrature = alpha * f0 * x / (1 + x**2)
    chi = in_phase * sp.cos(omega * t) + quadrature * sp.sin(omega * t)
    amplitude = alpha * f0 / sp.sqrt(1 + x**2)
    phase_lag = sp.atan(x)

    return RelaxationResponse(
        t=t,
        omega=omega,
        tau=tau,
        alpha=alpha,
        f0=f0,
        drive=drive,
        chi=sp.simplify(chi),
        in_phase=sp.simplify(in_phase),
        quadrature=sp.simplify(quadrature),
        amplitude=sp.simplify(amplitude),
        phase_lag=phase_lag,
    )


def relaxation_residual(response: RelaxationResponse | None = None) -> sp.Expr:
    response = response or steady_state_response()
    lhs = response.tau * sp.diff(response.chi, response.t) + response.chi
    rhs = response.alpha * response.drive
    return sp.trigsimp(sp.simplify(lhs - rhs))


def small_x_series(order: int = 5) -> dict[str, sp.Expr]:
    x = sp.symbols("x", positive=True, real=True)
    alpha, f0 = sp.symbols("alpha F0", positive=True, real=True)
    in_phase = alpha * f0 / (1 + x**2)
    quadrature = alpha * f0 * x / (1 + x**2)
    return {
        "in_phase": sp.series(in_phase, x, 0, order).removeO(),
        "quadrature": sp.series(quadrature, x, 0, order).removeO(),
    }


def response_summary() -> dict[str, sp.Expr]:
    response = steady_state_response()
    return {
        "drive": response.drive,
        "chi": response.chi,
        "in_phase": response.in_phase,
        "quadrature": response.quadrature,
        "amplitude": response.amplitude,
        "phase_lag": response.phase_lag,
        "transfer": transfer_function(
            alpha=response.alpha,
            omega=response.omega,
            tau=response.tau,
        ),
    }


def main() -> None:
    for key, value in response_summary().items():
        print(f"{key}: {sp.sstr(value)}")
    print(f"residual: {sp.sstr(relaxation_residual())}")


if __name__ == "__main__":
    main()
