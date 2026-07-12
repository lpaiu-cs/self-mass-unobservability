from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = REPO_ROOT / "outputs" / "tsv" / "nonlinear_sideband_test.tsv"
JSON_PATH = REPO_ROOT / "outputs" / "json" / "nonlinear_sideband_test.json"

TABLE_COLUMNS = (
    "term",
    "claim_status",
    "model",
    "generated_frequencies",
    "verdict",
    "collapse_condition",
    "projection_note",
    "notes",
)

ALLOWED_CLAIM_STATUSES = {
    "Proven",
    "Imported from prior work",
    "Conjectural",
    "Counterexample candidate",
}


@dataclass(frozen=True)
class SidebandRow:
    term: str
    claim_status: str
    model: str
    generated_frequencies: str
    verdict: str
    collapse_condition: str
    projection_note: str
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
            "model": self.model,
            "generated_frequencies": self.generated_frequencies,
            "verdict": self.verdict,
            "collapse_condition": self.collapse_condition,
            "projection_note": self.projection_note,
            "notes": self.notes,
        }


def symbols() -> dict[str, sp.Symbol]:
    names = (
        "Omega1 Omega2 tau_chi alpha F1 F2 beta_F2 c_chi "
        "lambda_Fchi lambda_chi2 p e n Lambda_side"
    )
    return {name: sp.Symbol(name, real=True) for name in names.split()}


def relaxation_component(frequency: sp.Expr, amplitude: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    sym = symbols()
    x = frequency * sym["tau_chi"]
    return (
        sp.simplify(amplitude / (1 + x**2)),
        sp.simplify(amplitude * x / (1 + x**2)),
    )


def linear_chi_components() -> dict[str, sp.Expr]:
    sym = symbols()
    omega1 = sym["Omega1"]
    omega2 = sym["Omega2"]
    tau = sym["tau_chi"]
    alpha = sym["alpha"]
    f1 = sym["F1"]
    f2 = sym["F2"]
    x1 = omega1 * tau
    x2 = omega2 * tau
    return {
        "A1": sp.simplify(alpha * f1 / (1 + x1**2)),
        "B1": sp.simplify(alpha * f1 * x1 / (1 + x1**2)),
        "A2": sp.simplify(alpha * f2 / (1 + x2**2)),
        "B2": sp.simplify(alpha * f2 * x2 / (1 + x2**2)),
    }


def nonlinear_drive_rhs_sidebands() -> dict[str, sp.Expr]:
    sym = symbols()
    return {
        "dc": sp.simplify(sym["beta_F2"] * (sym["F1"] ** 2 + sym["F2"] ** 2) / 2),
        "2Omega1": sp.simplify(sym["beta_F2"] * sym["F1"] ** 2 / 2),
        "2Omega2": sp.simplify(sym["beta_F2"] * sym["F2"] ** 2 / 2),
        "Omega1_plus_Omega2": sp.simplify(sym["beta_F2"] * sym["F1"] * sym["F2"]),
        "Omega1_minus_Omega2": sp.simplify(sym["beta_F2"] * sym["F1"] * sym["F2"]),
    }


def nonlinear_drive_chi_sidebands() -> dict[str, tuple[sp.Expr, sp.Expr]]:
    sym = symbols()
    rhs = nonlinear_drive_rhs_sidebands()
    return {
        "2Omega1": relaxation_component(2 * sym["Omega1"], rhs["2Omega1"]),
        "2Omega2": relaxation_component(2 * sym["Omega2"], rhs["2Omega2"]),
        "Omega1_plus_Omega2": relaxation_component(
            sym["Omega1"] + sym["Omega2"],
            rhs["Omega1_plus_Omega2"],
        ),
        "Omega1_minus_Omega2": relaxation_component(
            sym["Omega1"] - sym["Omega2"],
            rhs["Omega1_minus_Omega2"],
        ),
    }


def nonlinear_readout_sum_sideband() -> dict[str, sp.Expr]:
    sym = symbols()
    comp = linear_chi_components()
    a1, b1, a2, b2 = comp["A1"], comp["B1"], comp["A2"], comp["B2"]
    f1, f2 = sym["F1"], sym["F2"]
    lam_fchi = sym["lambda_Fchi"]
    lam_chi2 = sym["lambda_chi2"]
    return {
        "cos": sp.simplify(
            lam_fchi * (f1 * a2 + f2 * a1) / 2
            + lam_chi2 * (a1 * a2 - b1 * b2)
        ),
        "sin": sp.simplify(
            lam_fchi * (f1 * b2 + f2 * b1) / 2
            + lam_chi2 * (a1 * b2 + b1 * a2)
        ),
    }


def nonlinear_readout_difference_sideband() -> dict[str, sp.Expr]:
    sym = symbols()
    comp = linear_chi_components()
    a1, b1, a2, b2 = comp["A1"], comp["B1"], comp["A2"], comp["B2"]
    f1, f2 = sym["F1"], sym["F2"]
    lam_fchi = sym["lambda_Fchi"]
    lam_chi2 = sym["lambda_chi2"]
    return {
        "cos": sp.simplify(
            lam_fchi * (f1 * a2 + f2 * a1) / 2
            + lam_chi2 * (a1 * a2 + b1 * b2)
        ),
        "sin": sp.simplify(
            lam_fchi * (-f1 * b2 + f2 * b1) / 2
            + lam_chi2 * (-a1 * b2 + b1 * a2)
        ),
    }


def orbital_n_2n_input_amplitudes() -> dict[str, sp.Expr]:
    sym = symbols()
    p, e = sym["p"], sym["e"]
    return {
        "n": sp.simplify(p * e),
        "2n": sp.simplify(p * (p + 3) * e**2 / 4),
    }


def orbital_3n_nonlinear_drive_amplitude() -> sp.Expr:
    sym = symbols()
    amplitudes = orbital_n_2n_input_amplitudes()
    return sp.simplify(sym["beta_F2"] * amplitudes["n"] * amplitudes["2n"])


def linear_projection_sideband_output(input_sideband: sp.Expr | None = None) -> sp.Expr:
    sym = symbols()
    input_sideband = input_sideband if input_sideband is not None else sp.Integer(0)
    return sp.simplify(sym["Lambda_side"] * input_sideband)


def sideband_rows() -> list[SidebandRow]:
    return [
        SidebandRow(
            term="linear_baseline",
            claim_status="Proven",
            model="tau dot chi + chi = alpha F; linear readout",
            generated_frequencies="Omega1, Omega2 only",
            verdict="no sideband",
            collapse_condition="not applicable",
            projection_note="linear projection cannot create absent frequencies",
            notes="This is the baseline failure that motivates M5.",
        ),
        SidebandRow(
            term="nonlinear_drive_F2",
            claim_status="Proven",
            model="tau dot chi + chi = alpha F + beta_F2 F^2",
            generated_frequencies="2Omega1, 2Omega2, Omega1+Omega2, |Omega1-Omega2|",
            verdict="sidebands generated",
            collapse_condition="beta_F2 F1 F2 = 0 for mixed sidebands",
            projection_note="sidebands survive any nonzero linear projection at that frequency",
            notes="The relaxation transfer adds a phase lag at each generated frequency.",
        ),
        SidebandRow(
            term="nonlinear_readout_Fchi_chi2",
            claim_status="Proven",
            model="linear chi with lambda_Fchi F chi + lambda_chi2 chi^2 readout",
            generated_frequencies="Omega1+Omega2 and |Omega1-Omega2|",
            verdict="sidebands generated",
            collapse_condition="lambda_Fchi = lambda_chi2 = 0 or one harmonic amplitude vanishes",
            projection_note="linear projection can hide but not synthesize these frequencies",
            notes="The coefficients inherit A_j and B_j, so sideband phase carries tau_chi.",
        ),
        SidebandRow(
            term="orbital_3n_sideband",
            claim_status="Counterexample candidate",
            model="M3 orbital F=(a/r)^p plus beta_F2 F^2",
            generated_frequencies="3n from n+2n",
            verdict="new harmonic candidate",
            collapse_condition="beta_F2 = 0, e = 0, p = 0, or p = -3",
            projection_note="3n is absent from the O(e^2) linear drive",
            notes="This is the cleanest bridge from M3 forcing to M5 sideband.",
        ),
        SidebandRow(
            term="projection_robustness",
            claim_status="Proven",
            model="linear time-invariant projection",
            generated_frequencies="none by projection alone",
            verdict="sideband cannot be faked by linear projection of absent input",
            collapse_condition="projection is nonlinear or time-varying",
            projection_note="Lambda(i nu) only multiplies existing sideband amplitude",
            notes="This makes sidebands stronger than one-frequency quadrature anomalies.",
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in sideband_rows()]


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
    drive_chi_sidebands = nonlinear_drive_chi_sidebands()
    return {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "input": "F(t)=F1 cos(Omega1 t)+F2 cos(Omega2 t)",
            "nonlinear_drive": "tau_chi dot chi + chi = alpha F + beta_F2 F^2",
            "nonlinear_readout": "m/m0 includes lambda_Fchi F chi + lambda_chi2 chi^2",
            "scope": "symbolic minimal nonlinear sideband audit; no data fitting",
        },
        "nonlinear_drive_rhs_sidebands": {
            key: sp.sstr(value)
            for key, value in nonlinear_drive_rhs_sidebands().items()
        },
        "nonlinear_drive_chi_sidebands": {
            key: {"cos": sp.sstr(value[0]), "sin": sp.sstr(value[1])}
            for key, value in drive_chi_sidebands.items()
        },
        "nonlinear_readout_sum_sideband": {
            key: sp.sstr(value)
            for key, value in nonlinear_readout_sum_sideband().items()
        },
        "nonlinear_readout_difference_sideband": {
            key: sp.sstr(value)
            for key, value in nonlinear_readout_difference_sideband().items()
        },
        "orbital_n_2n_input_amplitudes": {
            key: sp.sstr(value)
            for key, value in orbital_n_2n_input_amplitudes().items()
        },
        "orbital_3n_nonlinear_drive_amplitude": sp.sstr(
            orbital_3n_nonlinear_drive_amplitude()
        ),
        "linear_projection_of_absent_sideband": sp.sstr(
            linear_projection_sideband_output()
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
    print(f"nonlinear_drive_rhs_sidebands: {data['nonlinear_drive_rhs_sidebands']}")
    print(f"orbital_3n_nonlinear_drive_amplitude: {data['orbital_3n_nonlinear_drive_amplitude']}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
