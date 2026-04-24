from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp

from chi_relaxation_response import transfer_function


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = REPO_ROOT / "outputs" / "tsv" / "dynamic_chi_basis_audit.tsv"
JSON_PATH = REPO_ROOT / "outputs" / "json" / "dynamic_chi_basis_audit.json"

ALLOWED_OBSERVABLE_STATUSES = {"observable", "absorbed", "degenerate", "null"}
ALLOWED_CLAIM_STATUSES = {
    "Proven",
    "Imported from prior work",
    "Conjectural",
    "Counterexample candidate",
}

AUDIT_COLUMNS = (
    "term",
    "claim_status",
    "drive_scope",
    "comparator_basis",
    "projection_result",
    "observable_status",
    "absorption_target",
    "notes",
)


@dataclass(frozen=True)
class BasisAuditRow:
    term: str
    claim_status: str
    drive_scope: str
    comparator_basis: str
    projection_result: str
    observable_status: str
    absorption_target: str
    notes: str

    def __post_init__(self) -> None:
        if self.claim_status not in ALLOWED_CLAIM_STATUSES:
            allowed = ", ".join(sorted(ALLOWED_CLAIM_STATUSES))
            raise ValueError(
                f"claim_status must be one of {allowed}; got {self.claim_status!r}"
            )
        if self.observable_status not in ALLOWED_OBSERVABLE_STATUSES:
            allowed = ", ".join(sorted(ALLOWED_OBSERVABLE_STATUSES))
            raise ValueError(
                "observable_status must be one of "
                f"{allowed}; got {self.observable_status!r}"
            )

    def to_row(self) -> dict[str, str]:
        return asdict(self)


def symbols() -> dict[str, sp.Symbol]:
    names = "Omega tau_chi alpha c_Y c_chi F0"
    return {
        name: sp.Symbol(name, positive=True, real=True)
        for name in names.split()
    }


def readout_transfer() -> sp.Expr:
    sym = symbols()
    return sym["c_Y"] + sym["c_chi"] * transfer_function(
        alpha=sym["alpha"],
        omega=sym["Omega"],
        tau=sym["tau_chi"],
    )


def monochromatic_readout_components() -> dict[str, sp.Expr]:
    sym = symbols()
    omega = sym["Omega"]
    tau = sym["tau_chi"]
    alpha = sym["alpha"]
    c_y = sym["c_Y"]
    c_chi = sym["c_chi"]
    f0 = sym["F0"]
    x = omega * tau
    return {
        "in_phase": sp.simplify(f0 * (c_y + c_chi * alpha / (1 + x**2))),
        "quadrature": sp.simplify(f0 * c_chi * alpha * x / (1 + x**2)),
        "quadrature_over_in_phase_chi_only": x,
    }


def static_only_residual() -> dict[str, sp.Expr]:
    components = monochromatic_readout_components()
    return {
        "matched_static_amplitude": components["in_phase"],
        "unmatched_quadrature": components["quadrature"],
    }


def single_frequency_derivative_fit() -> dict[str, sp.Expr]:
    """Return coefficients in a0 F + a1 dF/dt that match one drive frequency."""

    sym = symbols()
    omega = sym["Omega"]
    f0 = sym["F0"]
    components = monochromatic_readout_components()
    return {
        "a0": sp.simplify(components["in_phase"] / f0),
        "a1": sp.simplify(-components["quadrature"] / (omega * f0)),
    }


def derivative_series_transfer(order: int) -> sp.Expr:
    if order < 0:
        raise ValueError("order must be non-negative")
    sym = symbols()
    omega = sym["Omega"]
    tau = sym["tau_chi"]
    alpha = sym["alpha"]
    c_y = sym["c_Y"]
    c_chi = sym["c_chi"]
    chi_series = alpha * sum((-sp.I * omega * tau) ** n for n in range(order + 1))
    return sp.expand(c_y + c_chi * chi_series)


def derivative_series_residual(order: int) -> sp.Expr:
    return sp.simplify(readout_transfer() - derivative_series_transfer(order))


def polynomial_obstruction(order: int) -> dict[str, str]:
    """Show why a finite polynomial derivative transfer cannot equal the pole.

    The comparator is P_N(i Omega). Exact equality to alpha/(1+i Omega tau)
    would require (1+i Omega tau) P_N - alpha = 0 as a polynomial. The highest
    coefficient forces the highest derivative coefficient to vanish; recursion
    then forces all coefficients to vanish and leaves -alpha = 0.
    """

    if order < 0:
        raise ValueError("order must be non-negative")

    omega, tau, alpha = sp.symbols("Omega tau_chi alpha", nonzero=True)
    coeffs = sp.symbols(f"d0:{order + 1}")
    polynomial = sum(coeffs[n] * (sp.I * omega) ** n for n in range(order + 1))
    obstruction = sp.Poly(
        sp.expand((1 + sp.I * omega * tau) * polynomial - alpha),
        omega,
    )
    coefficients = obstruction.all_coeffs()
    return {
        "polynomial_identity": sp.sstr(obstruction.as_expr()),
        "highest_power_coefficient": sp.sstr(coefficients[0]),
        "constant_coefficient_after_recursion": "-alpha",
        "conclusion": (
            "no finite polynomial in i Omega matches alpha/(1+i Omega tau_chi) "
            "when alpha and tau_chi are nonzero"
        ),
    }


def audit_rows() -> list[BasisAuditRow]:
    return [
        BasisAuditRow(
            term="static_in_phase_amplitude",
            claim_status="Proven",
            drive_scope="single monochromatic drive",
            comparator_basis="instantaneous static F coefficient",
            projection_result="in-phase cos(Omega t) component is fit by one amplitude",
            observable_status="degenerate",
            absorption_target="static sensitivity coefficient c_Y_eff",
            notes=(
                "The cos component alone is not novel; it can be absorbed into "
                "a frequency-local static amplitude."
            ),
        ),
        BasisAuditRow(
            term="static_quadrature_phase_lag",
            claim_status="Counterexample candidate",
            drive_scope="single monochromatic drive with finite Omega tau_chi",
            comparator_basis="instantaneous static F coefficient only",
            projection_result="sin(Omega t) quadrature remains outside static basis",
            observable_status="observable",
            absorption_target="none in instantaneous static basis",
            notes=(
                "The unmatched quadrature is proportional to "
                "alpha c_chi F0 Omega tau_chi/(1+Omega^2 tau_chi^2)."
            ),
        ),
        BasisAuditRow(
            term="single_frequency_dotF_fit",
            claim_status="Proven",
            drive_scope="one known drive frequency",
            comparator_basis="local derivative basis {F, dot F} with free coefficients",
            projection_result="one derivative coefficient fits the quadrature exactly",
            observable_status="degenerate",
            absorption_target="frequency-local dot F Wilson coefficient",
            notes=(
                "At one frequency, a0 F + a1 dot F fits both quadratures. "
                "This does not prove a dynamical internal state."
            ),
        ),
        BasisAuditRow(
            term="adiabatic_derivative_series",
            claim_status="Proven",
            drive_scope="Omega tau_chi << 1",
            comparator_basis="finite local derivative EFT through chosen order N",
            projection_result="absorbed order-by-order up to the truncation error",
            observable_status="absorbed",
            absorption_target="derivative Wilson coefficients through order N",
            notes=(
                "The expansion of (1+tau_chi d/dt)^-1 reproduces the response "
                "to any fixed order after choosing N."
            ),
        ),
        BasisAuditRow(
            term="finite_tau_rational_transfer",
            claim_status="Counterexample candidate",
            drive_scope="frequency-dependent response or varying drive frequency",
            comparator_basis="finite static or finite-order polynomial derivative basis",
            projection_result="rational transfer with pole is outside any finite polynomial basis",
            observable_status="observable",
            absorption_target="none for finite shared-coefficient polynomial basis",
            notes=(
                "The surviving object is not one quadrature point but the shared "
                "transfer relation alpha c_chi/(1+i Omega tau_chi)."
            ),
        ),
        BasisAuditRow(
            term="tau_zero_collapse",
            claim_status="Proven",
            drive_scope="tau_chi = 0",
            comparator_basis="instantaneous static F coefficient",
            projection_result="quadrature vanishes and chi = alpha F",
            observable_status="null",
            absorption_target="static sensitivity redefinition",
            notes="The readout reduces to c_Y -> c_Y + alpha c_chi.",
        ),
        BasisAuditRow(
            term="omega_zero_collapse",
            claim_status="Proven",
            drive_scope="Omega = 0",
            comparator_basis="instantaneous static or secular basis",
            projection_result="phase lag vanishes",
            observable_status="null",
            absorption_target="static or secular offset",
            notes="Zero-frequency forcing has no delayed quadrature.",
        ),
        BasisAuditRow(
            term="alpha_cchi_zero_collapse",
            claim_status="Proven",
            drive_scope="alpha c_chi = 0",
            comparator_basis="any timing basis",
            projection_result="internal state is either undriven or unread",
            observable_status="null",
            absorption_target="not applicable",
            notes="There is no dynamic contribution in the measured readout.",
        ),
        BasisAuditRow(
            term="linear_two_frequency_sidebands",
            claim_status="Proven",
            drive_scope="linear two-frequency drive",
            comparator_basis="sideband search basis",
            projection_result="no Omega1+Omega2 or |Omega1-Omega2| rows are generated",
            observable_status="null",
            absorption_target="not applicable",
            notes="Linearity gives superposition only; sidebands need nonlinear input or readout.",
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in audit_rows()]


def validate_rows(rows: Iterable[dict[str, str]]) -> None:
    expected = set(AUDIT_COLUMNS)
    for index, row in enumerate(rows, start=1):
        actual = set(row)
        if actual != expected:
            raise ValueError(
                f"row {index} has columns {sorted(actual)}, "
                f"expected {sorted(expected)}"
            )
        if row["claim_status"] not in ALLOWED_CLAIM_STATUSES:
            raise ValueError(f"row {index} has invalid claim_status")
        if row["observable_status"] not in ALLOWED_OBSERVABLE_STATUSES:
            raise ValueError(f"row {index} has invalid observable_status")
        if not all(str(row[column]).strip() for column in AUDIT_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(AUDIT_COLUMNS),
        "assumptions": {
            "model": "tau_chi dot chi + chi = alpha F(t)",
            "readout": "delta(m/m0) = c_Y F + c_chi chi",
            "main_drive": "F(t) = F0 cos(Omega t)",
            "scope": "symbolic basis audit; no empirical runtime",
        },
        "readout_transfer": sp.sstr(readout_transfer()),
        "monochromatic_components": {
            key: sp.sstr(value)
            for key, value in monochromatic_readout_components().items()
        },
        "static_only_residual": {
            key: sp.sstr(value)
            for key, value in static_only_residual().items()
        },
        "single_frequency_derivative_fit": {
            key: sp.sstr(value)
            for key, value in single_frequency_derivative_fit().items()
        },
        "derivative_series_residual_order_3": sp.sstr(derivative_series_residual(3)),
        "polynomial_obstruction_order_3": polynomial_obstruction(3),
        "rows": rows,
    }


def write_outputs(data: dict[str, object] | None = None) -> None:
    data = data or payload()
    rows = data["rows"]
    if not isinstance(rows, list):
        raise TypeError("payload rows must be a list")

    TSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with TSV_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_COLUMNS, delimiter="\t")
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
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
