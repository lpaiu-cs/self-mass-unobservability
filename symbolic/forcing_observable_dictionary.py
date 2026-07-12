from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = REPO_ROOT / "outputs" / "tsv" / "forcing_observable_dictionary.tsv"
JSON_PATH = REPO_ROOT / "outputs" / "json" / "forcing_observable_dictionary.json"

TABLE_COLUMNS = (
    "term",
    "claim_status",
    "forcing_class",
    "frequency_source",
    "observable_projection",
    "verdict",
    "missing_assumption",
    "notes",
)

ALLOWED_CLAIM_STATUSES = {
    "Proven",
    "Imported from prior work",
    "Conjectural",
    "Counterexample candidate",
}


@dataclass(frozen=True)
class DictionaryRow:
    term: str
    claim_status: str
    forcing_class: str
    frequency_source: str
    observable_projection: str
    verdict: str
    missing_assumption: str
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
            "forcing_class": self.forcing_class,
            "frequency_source": self.frequency_source,
            "observable_projection": self.observable_projection,
            "verdict": self.verdict,
            "missing_assumption": self.missing_assumption,
            "notes": self.notes,
        }


def symbols() -> dict[str, sp.Symbol]:
    names = "p e n t Omega tau_chi beta c_Y F_k Lambda"
    return {
        name: sp.Symbol(name, real=True)
        for name in names.split()
    }


def power_law_orbital_harmonics_to_e2() -> dict[str, sp.Expr]:
    sym = symbols()
    p = sym["p"]
    e = sym["e"]
    n = sym["n"]
    t = sym["t"]
    phase = n * t
    return {
        "constant": 1 + p * (p - 1) * e**2 / 4,
        "first_harmonic": p * e * sp.cos(phase),
        "second_harmonic": p * (p + 3) * e**2 * sp.cos(2 * phase) / 4,
        "full_to_e2": (
            1
            + p * e * sp.cos(phase)
            + p * (p - 1) * e**2 / 4
            + p * (p + 3) * e**2 * sp.cos(2 * phase) / 4
        ),
    }


def power_law_series_residual_to_e2() -> sp.Expr:
    sym = symbols()
    p = sym["p"]
    e = sym["e"]
    n = sym["n"]
    t = sym["t"]
    phase = n * t

    r_over_a_to_e2 = 1 - e * sp.cos(phase) + e**2 * sp.sin(phase) ** 2
    derived = sp.series(r_over_a_to_e2 ** (-p), e, 0, 3).removeO()
    formula = power_law_orbital_harmonics_to_e2()["full_to_e2"]
    return sp.trigsimp(sp.expand(derived - formula))


def observable_harmonic_components() -> dict[str, sp.Expr]:
    sym = symbols()
    omega = sym["Omega"]
    tau = sym["tau_chi"]
    beta = sym["beta"]
    c_y = sym["c_Y"]
    f_k = sym["F_k"]
    projection = sym["Lambda"]
    x = omega * tau
    return {
        "cos_coefficient": sp.simplify(
            projection * f_k * (c_y + beta / (1 + x**2))
        ),
        "sin_coefficient": sp.simplify(
            projection * f_k * beta * x / (1 + x**2)
        ),
        "calibrated_cos": sp.simplify(c_y + beta / (1 + x**2)),
        "calibrated_sin": sp.simplify(beta * x / (1 + x**2)),
    }


def orbital_two_harmonic_dotf_obstruction() -> sp.Expr:
    sym = symbols()
    n = sym["n"]
    tau = sym["tau_chi"]
    beta = sym["beta"]
    u0 = n**2
    u_star = (2 * n) ** 2
    return sp.factor(
        beta
        * tau**3
        * (u_star - u0)
        / ((1 + tau**2 * u0) * (1 + tau**2 * u_star))
    )


def dictionary_rows() -> list[DictionaryRow]:
    return [
        DictionaryRow(
            term="orbital_power_law_harmonics",
            claim_status="Counterexample candidate",
            forcing_class="normalized invariant Y(r)/Y_* = (a/r)^p",
            frequency_source="single eccentric orbit gives n and 2n through O(e^2)",
            observable_projection="linear channel O_hat = Lambda G F_hat",
            verdict="concrete single-system sweep candidate",
            missing_assumption="nonzero e, p, p+3, beta, tau_chi, and Lambda",
            notes=(
                "This supplies the first concrete F(Y(t)) template without "
                "using data or clock-sector assumptions."
            ),
        ),
        DictionaryRow(
            term="calibrated_linear_projection",
            claim_status="Counterexample candidate",
            forcing_class="any known harmonic drive F_k cos(Omega_k t)",
            frequency_source="one or more known drive frequencies",
            observable_projection="known nonzero Lambda or finite-dimensional Lambda(Omega)",
            verdict="preserves the pole relation after calibration",
            missing_assumption="projection must not be arbitrary per frequency",
            notes=(
                "If Lambda is common, O_hat/(Lambda F_hat) recovers G. If "
                "Lambda is a finite polynomial nuisance, raise the comparator "
                "order accordingly."
            ),
        ),
        DictionaryRow(
            term="multi_system_sweep",
            claim_status="Conjectural",
            forcing_class="same drive construction across several systems",
            frequency_source="different orbital frequencies across systems",
            observable_projection="system-specific but modeled projection",
            verdict="usable only with shared tau_chi relation",
            missing_assumption="common tau_chi or independently specified body-dependence",
            notes=(
                "Without shared parameters, each system can be fit by its own "
                "effective derivative coefficients."
            ),
        ),
        DictionaryRow(
            term="independent_two_tone_forcing",
            claim_status="Counterexample candidate",
            forcing_class="F = F1 cos(Omega1 t) + F2 cos(Omega2 t)",
            frequency_source="two known distinct drives in one system",
            observable_projection="linear channel O_hat = Lambda G F_hat",
            verdict="two transfer samples but no sidebands",
            missing_assumption="known distinct drives and nonzero projection",
            notes=(
                "Enough to beat real degree-1 or degree-2 derivative EFTs, "
                "not arbitrary higher finite order."
            ),
        ),
        DictionaryRow(
            term="arbitrary_projection_failure",
            claim_status="Proven",
            forcing_class="any",
            frequency_source="any finite sampled set",
            observable_projection="unconstrained complex Lambda(Omega_k)",
            verdict="dictionary collapses",
            missing_assumption="finite or calibrated observable projection",
            notes=(
                "An arbitrary frequency-local projection nuisance can absorb "
                "the pole relation."
            ),
        ),
        DictionaryRow(
            term="circular_or_flat_drive_failure",
            claim_status="Proven",
            forcing_class="power-law orbital harmonic template",
            frequency_source="e=0, p=0, or p=-3 at O(e^2)",
            observable_projection="any",
            verdict="insufficient harmonic sweep",
            missing_assumption="nonzero harmonic content in selected invariant",
            notes=(
                "No eccentric modulation, no first harmonic, or no second "
                "harmonic at this order removes the two-point test."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in dictionary_rows()]


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


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "drive_template": "F(Y(t)) = (a/r(t))^p through O(e^2)",
            "observable_template": "O_hat(Omega_k) = Lambda(Omega_k) G(i Omega_k) F_hat(Omega_k)",
            "transfer": "G(z) = c_Y + beta/(1 + tau_chi z)",
            "scope": "symbolic forcing/observable dictionary; no data fitting",
        },
        "power_law_orbital_harmonics_to_e2": {
            key: sp.sstr(value)
            for key, value in power_law_orbital_harmonics_to_e2().items()
        },
        "power_law_series_residual_to_e2": sp.sstr(power_law_series_residual_to_e2()),
        "observable_harmonic_components": {
            key: sp.sstr(value)
            for key, value in observable_harmonic_components().items()
        },
        "orbital_two_harmonic_dotf_obstruction": sp.sstr(
            orbital_two_harmonic_dotf_obstruction()
        ),
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
    print(f"harmonics: {data['power_law_orbital_harmonics_to_e2']}")
    print(f"projection: {data['observable_harmonic_components']}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
