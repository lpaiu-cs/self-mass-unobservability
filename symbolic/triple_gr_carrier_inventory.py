from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from frequency_sweep_distinguishability import (
    ALLOWED_CLAIM_STATUSES,
    polynomial_sample_limit,
    real_coefficient_sample_limit,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = REPO_ROOT / "outputs" / "tsv" / "dynamic_chi_triple_gr_carrier_inventory.tsv"
JSON_PATH = REPO_ROOT / "outputs" / "json" / "dynamic_chi_triple_gr_carrier_inventory.json"

INVENTORY_COLUMNS = (
    "term",
    "claim_status",
    "carrier_set",
    "N_frequencies",
    "comparator_basis",
    "verdict",
    "surviving_target",
    "notes",
)


@dataclass(frozen=True)
class CarrierInventoryRow:
    term: str
    claim_status: str
    carrier_set: str
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


def carrier_inventory() -> dict[str, str]:
    return {
        "inner_monopole": "Omega_in",
        "outer_monopole": "Omega_out",
        "outer_dipole_combination": "abs(Omega_in - Omega_out)",
        "eccentric_outer_dipole_sidebands": (
            "abs(2 Omega_in - Omega_out), abs(Omega_in - 2 Omega_out)"
        ),
    }


def three_carrier_conditions() -> list[str]:
    return [
        "Omega_in > 0",
        "Omega_out > 0",
        "Omega_in != Omega_out",
        "Omega_in != 2 Omega_out",
        "Omega_out != 2 Omega_in",
        "inner, outer, and outer-dipole projection amplitudes are nonzero",
        "projection nuisance is finite-dimensional, not arbitrary per carrier",
    ]


def real_k_frequency_verdict(order: int, frequency_count: int) -> dict[str, str | int]:
    boundary = real_coefficient_sample_limit(order)
    obstruction_count = boundary["first_positive_frequency_obstruction_count"]
    if frequency_count >= obstruction_count:
        verdict = "distinguishable"
        target = "shared tau_chi odd-channel pole residual"
    else:
        verdict = f"degenerate at {frequency_count} carriers"
        target = "none without more distinct carrier frequencies or order prior"
    return {
        "degree": order,
        "N_frequencies": frequency_count,
        "first_obstruction_count": obstruction_count,
        "verdict": verdict,
        "surviving_target": target,
    }


def complex_k_frequency_verdict(order: int, frequency_count: int) -> dict[str, str | int]:
    boundary = polynomial_sample_limit(order)
    obstruction_count = boundary["first_exact_obstruction_count"]
    if frequency_count >= obstruction_count:
        verdict = "distinguishable"
        target = "full rational transfer pole residual"
    else:
        verdict = f"degenerate at {frequency_count} carriers"
        target = "none without more distinct carrier frequencies or order prior"
    return {
        "degree": order,
        "N_frequencies": frequency_count,
        "first_obstruction_count": obstruction_count,
        "verdict": verdict,
        "surviving_target": target,
    }


def inventory_rows() -> list[CarrierInventoryRow]:
    real_n4 = real_k_frequency_verdict(4, 3)
    real_n5 = real_k_frequency_verdict(5, 3)
    complex_n1 = complex_k_frequency_verdict(1, 3)
    complex_n2 = complex_k_frequency_verdict(2, 3)
    return [
        CarrierInventoryRow(
            term="inner_outer_combination_inventory",
            claim_status="Counterexample candidate",
            carrier_set="Omega_in, Omega_out, abs(Omega_in - Omega_out)",
            N_frequencies="3 under generic nonresonant hierarchy",
            comparator_basis="frequency inventory only",
            verdict="third GR carrier available conditionally",
            surviving_target="three-sample shared tau_chi bridge input",
            notes=(
                "The third carrier is the existing GR outer-dipole combination "
                "frequency. It need not be a new clock-sector shape."
            ),
        ),
        CarrierInventoryRow(
            term="inner_outer_pair_baseline",
            claim_status="Proven",
            carrier_set="Omega_in, Omega_out",
            N_frequencies="2",
            comparator_basis="real shared-coefficient derivative EFT",
            verdict="distinguishes degree 1 and 2 only",
            surviving_target="low-order shared tau_chi pole residual",
            notes=(
                "This restates Request 10.2b: two carriers are useful but not "
                "enough for higher derivative order."
            ),
        ),
        CarrierInventoryRow(
            term="three_carrier_real_degree4",
            claim_status="Counterexample candidate",
            carrier_set="Omega_in, Omega_out, abs(Omega_in - Omega_out)",
            N_frequencies="3",
            comparator_basis="real shared-coefficient degree-N derivative EFT, N <= 4",
            verdict=str(real_n4["verdict"]),
            surviving_target=str(real_n4["surviving_target"]),
            notes=(
                "Three distinct positive carriers reach the first obstruction "
                "count for real degree 3 and 4 comparators."
            ),
        ),
        CarrierInventoryRow(
            term="three_carrier_real_degree5_or_higher",
            claim_status="Proven",
            carrier_set="Omega_in, Omega_out, abs(Omega_in - Omega_out)",
            N_frequencies="3",
            comparator_basis="real shared-coefficient degree-N derivative EFT, N >= 5",
            verdict=str(real_n5["verdict"]),
            surviving_target=str(real_n5["surviving_target"]),
            notes=(
                "Degree 5 raises the real positive-frequency obstruction count "
                "to four samples, so the three-carrier bridge is insufficient."
            ),
        ),
        CarrierInventoryRow(
            term="three_carrier_complex_degree1",
            claim_status="Counterexample candidate",
            carrier_set="Omega_in, Omega_out, abs(Omega_in - Omega_out)",
            N_frequencies="3",
            comparator_basis="complex shared-coefficient degree-N polynomial, N <= 1",
            verdict=str(complex_n1["verdict"]),
            surviving_target=str(complex_n1["surviving_target"]),
            notes=(
                "Three samples are enough to beat the deliberately generous "
                "complex degree-1 polynomial comparator."
            ),
        ),
        CarrierInventoryRow(
            term="three_carrier_complex_degree2_or_higher",
            claim_status="Proven",
            carrier_set="Omega_in, Omega_out, abs(Omega_in - Omega_out)",
            N_frequencies="3",
            comparator_basis="complex shared-coefficient degree-N polynomial, N >= 2",
            verdict=str(complex_n2["verdict"]),
            surviving_target=str(complex_n2["surviving_target"]),
            notes=(
                "Complex degree 2 needs four distinct samples. The conservative "
                "complex gate remains open without another carrier or an order prior."
            ),
        ),
        CarrierInventoryRow(
            term="eccentric_outer_dipole_sideband_extension",
            claim_status="Counterexample candidate",
            carrier_set="abs(2 Omega_in - Omega_out), abs(Omega_in - 2 Omega_out)",
            N_frequencies="optional fourth/fifth samples if eccentric sidebands are usable",
            comparator_basis="same shared-coefficient derivative gates",
            verdict="potential extension",
            surviving_target="higher-order comparator pressure",
            notes=(
                "The Request 8 triple rank audit showed these sidebands in the "
                "GR carrier vector. They are not new scalar-EFT shapes, but can "
                "serve as additional frequency samples if their projections are "
                "not arbitrary nuisance amplitudes."
            ),
        ),
        CarrierInventoryRow(
            term="resonance_or_projection_collapse",
            claim_status="Proven",
            carrier_set="commensurate or unprojected carriers",
            N_frequencies="less than 3",
            comparator_basis="any three-carrier bridge comparator",
            verdict="collapses to two-carrier or one-carrier bridge",
            surviving_target="none beyond lower inventory result",
            notes=(
                "The third sample is lost if Omega_in=Omega_out, "
                "Omega_in=2 Omega_out, Omega_out=2 Omega_in, or the "
                "outer-dipole projection amplitude is zero/arbitrary."
            ),
        ),
        CarrierInventoryRow(
            term="arbitrary_projection_failure",
            claim_status="Proven",
            carrier_set="any finite carrier set",
            N_frequencies="any",
            comparator_basis="unconstrained complex projection nuisance per carrier",
            verdict="collapses",
            surviving_target="none",
            notes=(
                "If each carrier has an arbitrary complex projection amplitude, "
                "the projection absorbs the shared transfer law point by point."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in inventory_rows()]


def validate_rows(rows: Iterable[dict[str, str]]) -> None:
    expected = set(INVENTORY_COLUMNS)
    for index, row in enumerate(rows, start=1):
        actual = set(row)
        if actual != expected:
            raise ValueError(
                f"row {index} has columns {sorted(actual)}, "
                f"expected {sorted(expected)}"
            )
        if row["claim_status"] not in ALLOWED_CLAIM_STATUSES:
            raise ValueError(f"row {index} has invalid claim_status")
        if not all(str(row[column]).strip() for column in INVENTORY_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(INVENTORY_COLUMNS),
        "assumptions": {
            "carrier_source": "existing GR hierarchical-triple clock carriers",
            "base_carriers": "inner monopole, outer monopole, outer-dipole combination",
            "scope": "symbolic inventory and comparator-count audit; no runtime",
        },
        "carrier_inventory": carrier_inventory(),
        "three_carrier_conditions": three_carrier_conditions(),
        "real_degree_boundaries_for_three_carriers": {
            str(order): real_k_frequency_verdict(order, 3)
            for order in range(0, 7)
        },
        "complex_degree_boundaries_for_three_carriers": {
            str(order): complex_k_frequency_verdict(order, 3)
            for order in range(0, 5)
        },
        "rows": rows,
    }


def write_outputs(data: dict[str, object] | None = None) -> None:
    data = data or payload()
    rows = data["rows"]
    if not isinstance(rows, list):
        raise TypeError("payload rows must be a list")

    TSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with TSV_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=INVENTORY_COLUMNS, delimiter="\t")
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
    print(f"carrier_inventory: {data['carrier_inventory']}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
