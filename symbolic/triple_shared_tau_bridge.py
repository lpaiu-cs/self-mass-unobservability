from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp

from frequency_sweep_distinguishability import (
    ALLOWED_CLAIM_STATUSES,
    polynomial_sample_limit,
    real_coefficient_sample_limit,
    real_odd_channel_residual,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = REPO_ROOT / "outputs" / "tsv" / "dynamic_chi_triple_shared_tau_bridge.tsv"
JSON_PATH = REPO_ROOT / "outputs" / "json" / "dynamic_chi_triple_shared_tau_bridge.json"

BRIDGE_COLUMNS = (
    "term",
    "claim_status",
    "carrier_inventory",
    "N_frequencies",
    "comparator_basis",
    "verdict",
    "surviving_target",
    "notes",
)


@dataclass(frozen=True)
class TripleBridgeRow:
    term: str
    claim_status: str
    carrier_inventory: str
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
        return asdict(self)


def symbols() -> dict[str, sp.Symbol]:
    names = "Omega_in Omega_out tau_chi beta c_Y F_in F_out phi_in phi_out t"
    return {
        name: sp.Symbol(name, positive=True, real=True)
        for name in names.split()
    }


def forcing_model() -> dict[str, sp.Expr]:
    sym = symbols()
    omega_in = sym["Omega_in"]
    omega_out = sym["Omega_out"]
    f_in = sym["F_in"]
    f_out = sym["F_out"]
    phi_in = sym["phi_in"]
    phi_out = sym["phi_out"]
    t = sym["t"]
    return {
        "inner_carrier": f_in * sp.cos(omega_in * t + phi_in),
        "outer_carrier": f_out * sp.cos(omega_out * t + phi_out),
        "total": (
            f_in * sp.cos(omega_in * t + phi_in)
            + f_out * sp.cos(omega_out * t + phi_out)
        ),
    }


def transfer_at_carriers() -> dict[str, sp.Expr]:
    sym = symbols()
    beta = sym["beta"]
    tau = sym["tau_chi"]
    c_y = sym["c_Y"]
    omega_in = sym["Omega_in"]
    omega_out = sym["Omega_out"]
    return {
        "G_in": c_y + beta / (1 + sp.I * omega_in * tau),
        "G_out": c_y + beta / (1 + sp.I * omega_out * tau),
        "difference": sp.simplify(
            beta / (1 + sp.I * omega_in * tau)
            - beta / (1 + sp.I * omega_out * tau)
        ),
    }


def real_two_frequency_verdict(order: int) -> dict[str, str | int]:
    boundary = real_coefficient_sample_limit(order)
    frequency_count = 2
    obstruction_count = boundary["first_positive_frequency_obstruction_count"]
    if frequency_count >= obstruction_count:
        verdict = "distinguishable"
        target = "shared tau_chi odd-channel pole residual"
    else:
        verdict = "degenerate at two carriers"
        target = "none without more distinct carrier frequencies"
    return {
        "degree": order,
        "N_frequencies": frequency_count,
        "first_obstruction_count": obstruction_count,
        "verdict": verdict,
        "surviving_target": target,
    }


def complex_two_frequency_verdict(order: int) -> dict[str, str | int]:
    boundary = polynomial_sample_limit(order)
    frequency_count = 2
    obstruction_count = boundary["first_exact_obstruction_count"]
    if frequency_count >= obstruction_count:
        verdict = "distinguishable"
        target = "full rational transfer pole residual"
    else:
        verdict = "degenerate at two carriers"
        target = "none without more distinct carrier frequencies"
    return {
        "degree": order,
        "N_frequencies": frequency_count,
        "first_obstruction_count": obstruction_count,
        "verdict": verdict,
        "surviving_target": target,
    }


def two_carrier_real_odd_residual_order1() -> sp.Expr:
    return real_odd_channel_residual(1)


def bridge_rows() -> list[TripleBridgeRow]:
    real_n1 = real_two_frequency_verdict(1)
    real_n2 = real_two_frequency_verdict(2)
    real_n3 = real_two_frequency_verdict(3)
    complex_n0 = complex_two_frequency_verdict(0)
    complex_n1 = complex_two_frequency_verdict(1)
    return [
        TripleBridgeRow(
            term="triple_inner_outer_inventory",
            claim_status="Proven",
            carrier_inventory="Omega_in, Omega_out",
            N_frequencies="2 if Omega_in != Omega_out and both amplitudes are nonzero",
            comparator_basis="frequency inventory only",
            verdict="two existing GR carrier samples available",
            surviving_target="shared tau_chi transfer-law test input",
            notes=(
                "The carriers need not be new tensor harmonics. Dynamic-chi uses "
                "existing distinct carrier frequencies as samples of the same "
                "relaxation transfer."
            ),
        ),
        TripleBridgeRow(
            term="single_carrier_local_derivative",
            claim_status="Proven",
            carrier_inventory="either Omega_in or Omega_out alone",
            N_frequencies="1",
            comparator_basis="frequency-local {F, dot F}",
            verdict="degenerate",
            surviving_target="none",
            notes=(
                "Each carrier by itself has the Request 10.1 one-frequency "
                "degeneracy."
            ),
        ),
        TripleBridgeRow(
            term="two_carrier_real_degree1",
            claim_status="Counterexample candidate",
            carrier_inventory="Omega_in, Omega_out",
            N_frequencies="2",
            comparator_basis="real shared-coefficient degree-1 derivative EFT",
            verdict=str(real_n1["verdict"]),
            surviving_target=str(real_n1["surviving_target"]),
            notes=(
                "Two distinct positive frequencies are the first obstruction "
                "for a real degree-1 comparator."
            ),
        ),
        TripleBridgeRow(
            term="two_carrier_real_degree2",
            claim_status="Counterexample candidate",
            carrier_inventory="Omega_in, Omega_out",
            N_frequencies="2",
            comparator_basis="real shared-coefficient degree-2 derivative EFT",
            verdict=str(real_n2["verdict"]),
            surviving_target=str(real_n2["surviving_target"]),
            notes=(
                "The odd channel still has only one coefficient through degree 2, "
                "so the second positive frequency exposes the pole relation."
            ),
        ),
        TripleBridgeRow(
            term="two_carrier_real_degree3_or_higher",
            claim_status="Proven",
            carrier_inventory="Omega_in, Omega_out",
            N_frequencies="2",
            comparator_basis="real shared-coefficient degree-N derivative EFT, N >= 3",
            verdict=str(real_n3["verdict"]),
            surviving_target=str(real_n3["surviving_target"]),
            notes=(
                "Two carriers are insufficient once the real derivative EFT has "
                "enough odd-channel coefficients. More distinct carrier "
                "frequencies are required."
            ),
        ),
        TripleBridgeRow(
            term="two_carrier_complex_degree0",
            claim_status="Counterexample candidate",
            carrier_inventory="Omega_in, Omega_out",
            N_frequencies="2",
            comparator_basis="complex shared-coefficient degree-0 polynomial",
            verdict=str(complex_n0["verdict"]),
            surviving_target=str(complex_n0["surviving_target"]),
            notes=(
                "This overpowered complex static comparator is still broken by "
                "two distinct samples."
            ),
        ),
        TripleBridgeRow(
            term="two_carrier_complex_degree1_or_higher",
            claim_status="Proven",
            carrier_inventory="Omega_in, Omega_out",
            N_frequencies="2",
            comparator_basis="complex shared-coefficient degree-N polynomial, N >= 1",
            verdict=str(complex_n1["verdict"]),
            surviving_target=str(complex_n1["surviving_target"]),
            notes=(
                "An overpowered complex degree-1 polynomial can interpolate two "
                "samples, so the conservative complex-comparator gate needs at "
                "least three carriers."
            ),
        ),
        TripleBridgeRow(
            term="commensurate_or_missing_carrier_collapse",
            claim_status="Proven",
            carrier_inventory="Omega_in = Omega_out or F_in F_out = 0",
            N_frequencies="less than 2",
            comparator_basis="any multi-frequency bridge comparator",
            verdict="collapses to one-frequency degeneracy",
            surviving_target="none",
            notes=(
                "The bridge requires two distinct positive carrier frequencies "
                "with nonzero readout amplitudes."
            ),
        ),
        TripleBridgeRow(
            term="linear_sideband_check",
            claim_status="Proven",
            carrier_inventory="Omega_in, Omega_out",
            N_frequencies="2",
            comparator_basis="sideband search basis",
            verdict="null",
            surviving_target="none",
            notes=(
                "The linear one-state model still produces no Omega_in +/- "
                "Omega_out sidebands. The bridge is transfer-law based."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in bridge_rows()]


def validate_rows(rows: Iterable[dict[str, str]]) -> None:
    expected = set(BRIDGE_COLUMNS)
    for index, row in enumerate(rows, start=1):
        actual = set(row)
        if actual != expected:
            raise ValueError(
                f"row {index} has columns {sorted(actual)}, "
                f"expected {sorted(expected)}"
            )
        if row["claim_status"] not in ALLOWED_CLAIM_STATUSES:
            raise ValueError(f"row {index} has invalid claim_status")
        if not all(str(row[column]).strip() for column in BRIDGE_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(BRIDGE_COLUMNS),
        "assumptions": {
            "model": "tau_chi dot chi + chi = alpha F(t)",
            "forcing": (
                "F_in cos(Omega_in t + phi_in) + "
                "F_out cos(Omega_out t + phi_out)"
            ),
            "carrier_source": "existing hierarchical-triple inner/outer GR clock carriers",
            "scope": "symbolic carrier bridge; no triple timing runtime",
        },
        "forcing_model": {
            key: sp.sstr(value)
            for key, value in forcing_model().items()
        },
        "transfer_at_carriers": {
            key: sp.sstr(value)
            for key, value in transfer_at_carriers().items()
        },
        "real_degree_boundaries": {
            str(order): real_two_frequency_verdict(order)
            for order in range(0, 5)
        },
        "complex_degree_boundaries": {
            str(order): complex_two_frequency_verdict(order)
            for order in range(0, 4)
        },
        "two_carrier_real_odd_residual_order1": sp.sstr(
            two_carrier_real_odd_residual_order1()
        ),
        "collapse_conditions": [
            "Omega_in = Omega_out",
            "F_in = 0 or F_out = 0",
            "beta = alpha c_chi = 0",
            "tau_chi = 0",
            "low-frequency band below derivative truncation tolerance",
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
        writer = csv.DictWriter(handle, fieldnames=BRIDGE_COLUMNS, delimiter="\t")
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
    print(f"forcing: {data['forcing_model']['total']}")
    print(f"G_in: {data['transfer_at_carriers']['G_in']}")
    print(f"G_out: {data['transfer_at_carriers']['G_out']}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
