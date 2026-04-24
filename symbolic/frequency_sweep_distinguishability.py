from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = REPO_ROOT / "outputs" / "tsv" / "dynamic_chi_multifrequency_audit.tsv"
JSON_PATH = REPO_ROOT / "outputs" / "json" / "dynamic_chi_multifrequency_audit.json"

TABLE_COLUMNS = (
    "term",
    "claim_status",
    "derivative_order",
    "N_frequencies",
    "comparator_basis",
    "verdict",
    "surviving_target",
    "notes",
)

ALLOWED_CLAIM_STATUSES = {
    "Proven",
    "Imported from prior work",
    "Conjectural",
    "Counterexample candidate",
}


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


@dataclass(frozen=True)
class MultiFrequencyAuditRow:
    term: str
    claim_status: str
    derivative_order: str
    N_frequencies: str
    comparator_basis: str
    verdict: str
    surviving_target: str
    notes: str

    def __post_init__(self) -> None:
        if self.claim_status not in ALLOWED_CLAIM_STATUSES:
            allowed = ", ".join(sorted(ALLOWED_CLAIM_STATUSES))
            raise ValueError(
                f"claim_status must be one of {allowed}; got {self.claim_status!r}"
            )

    def to_row(self) -> dict[str, str]:
        return {
            "term": self.term,
            "claim_status": self.claim_status,
            "derivative_order": self.derivative_order,
            "N_frequencies": self.N_frequencies,
            "comparator_basis": self.comparator_basis,
            "verdict": self.verdict,
            "surviving_target": self.surviving_target,
            "notes": self.notes,
        }


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


def audit_rows() -> list[MultiFrequencyAuditRow]:
    return [
        MultiFrequencyAuditRow(
            term="single_frequency_static_only",
            claim_status="Counterexample candidate",
            derivative_order="not applicable",
            N_frequencies="1",
            comparator_basis="instantaneous static F coefficient",
            verdict="distinguishable from static-only basis",
            surviving_target="finite-Omega tau_chi quadrature",
            notes=(
                "This is the weakest target. It does not survive once a "
                "frequency-local dot F coefficient is allowed."
            ),
        ),
        MultiFrequencyAuditRow(
            term="single_frequency_dotF",
            claim_status="Proven",
            derivative_order="1",
            N_frequencies="1",
            comparator_basis="frequency-local {F, dot F}",
            verdict="degenerate",
            surviving_target="none",
            notes=(
                "One sine quadrature can be fit exactly by one unconstrained "
                "dot F Wilson coefficient at the known frequency."
            ),
        ),
        MultiFrequencyAuditRow(
            term="real_shared_coefficients_interpolation_zone",
            claim_status="Proven",
            derivative_order="N",
            N_frequencies="K <= floor((N+1)/2)",
            comparator_basis="real shared-coefficient degree-N derivative EFT",
            verdict="degenerate by even/odd-channel interpolation",
            surviving_target="none at finite sampled set",
            notes=(
                "For positive frequencies, real coefficients split into an even "
                "channel and an odd channel in u=Omega^2; the odd channel is "
                "the limiting interpolation count."
            ),
        ),
        MultiFrequencyAuditRow(
            term="real_shared_coefficients_first_obstruction",
            claim_status="Counterexample candidate",
            derivative_order="N",
            N_frequencies="floor((N+1)/2)+1",
            comparator_basis="real shared-coefficient degree-N derivative EFT",
            verdict="distinguishable for distinct positive frequencies",
            surviving_target="shared tau_chi odd-channel pole residual",
            notes=(
                "The first exact physical obstruction is the interpolation "
                "residual of -beta tau_chi/(1+tau_chi^2 Omega^2)."
            ),
        ),
        MultiFrequencyAuditRow(
            term="complex_shared_coefficients_interpolation_zone",
            claim_status="Proven",
            derivative_order="N",
            N_frequencies="K <= N+1",
            comparator_basis="complex shared-coefficient degree-N polynomial",
            verdict="degenerate by polynomial interpolation",
            surviving_target="none at finite sampled set",
            notes=(
                "This is an intentionally overpowered comparator. It gives a "
                "conservative no-distinguishability region."
            ),
        ),
        MultiFrequencyAuditRow(
            term="complex_shared_coefficients_first_obstruction",
            claim_status="Counterexample candidate",
            derivative_order="N",
            N_frequencies="N+2",
            comparator_basis="complex shared-coefficient degree-N polynomial",
            verdict="distinguishable for distinct off-pole samples",
            surviving_target="full rational transfer pole residual",
            notes=(
                "The residual is beta(-tau_chi)^(N+1) prod(z*-zj) divided by "
                "the sampled pole denominators."
            ),
        ),
        MultiFrequencyAuditRow(
            term="low_frequency_taylor_band",
            claim_status="Proven",
            derivative_order="N",
            N_frequencies="any within |Omega tau_chi| <= rho",
            comparator_basis="Taylor derivative EFT through order N",
            verdict="operationally degenerate below truncation tolerance",
            surviving_target="none if tolerance >= |beta| rho^(N+1)",
            notes=(
                "The exact Taylor residual is beta(-tau_chi z)^(N+1)/"
                "(1+tau_chi z)."
            ),
        ),
        MultiFrequencyAuditRow(
            term="linear_two_frequency_sideband_check",
            claim_status="Proven",
            derivative_order="not applicable",
            N_frequencies="2",
            comparator_basis="sideband search basis",
            verdict="null",
            surviving_target="none",
            notes=(
                "Linear two-frequency forcing produces only the input "
                "frequencies. Multi-frequency novelty is transfer-law "
                "distinguishability, not sideband creation, in this MVP."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in audit_rows()]


def validate_rows(rows: Iterable[dict[str, str]]) -> None:
    expected = set(TABLE_COLUMNS)
    for index, row in enumerate(rows, start=1):
        actual = set(row)
        if actual != expected:
            raise ValueError(
                f"row {index} has columns {sorted(actual)}, "
                f"expected {sorted(expected)}"
            )
        if row["claim_status"] not in ALLOWED_CLAIM_STATUSES:
            raise ValueError(f"row {index} has invalid claim_status")
        if not all(str(row[column]).strip() for column in TABLE_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


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


def payload(order: int = 3) -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "model": "tau_chi dot chi + chi = alpha F(t)",
            "readout_transfer": "G(z) = c_Y + beta/(1 + tau_chi z)",
            "beta": "alpha c_chi",
            "z": "i Omega",
            "scope": "symbolic multi-frequency audit; no empirical runtime",
        },
        "sample_order_used_for_formulas": order,
        "formula_summary": summary(order),
        "collapse_conditions": [
            "beta = alpha c_chi = 0",
            "tau_chi = 0",
            "sampled frequencies not distinct",
            "complex-sample pole condition 1 + tau_chi z_j = 0",
            "operational low-frequency tolerance above Taylor residual bound",
        ],
        "rows": rows,
    }


def write_outputs(data: dict[str, object] | None = None) -> None:
    data = data or payload()
    rows = data["rows"]
    if not isinstance(rows, list):
        raise TypeError("payload rows must be a list")

    TSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with TSV_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TABLE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    data = payload()
    write_outputs(data)
    for key, value in data["formula_summary"].items():
        print(f"{key}: {value}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
