from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = REPO_ROOT / "outputs" / "tsv" / "projection_channel_audit.tsv"
JSON_PATH = REPO_ROOT / "outputs" / "json" / "projection_channel_audit.json"

TABLE_COLUMNS = (
    "term",
    "claim_status",
    "projection_channel",
    "lambda_form",
    "verdict",
    "collapse_condition",
    "next_requirement",
    "notes",
)

ALLOWED_CLAIM_STATUSES = {
    "Proven",
    "Imported from prior work",
    "Conjectural",
    "Counterexample candidate",
}


@dataclass(frozen=True)
class ProjectionRow:
    term: str
    claim_status: str
    projection_channel: str
    lambda_form: str
    verdict: str
    collapse_condition: str
    next_requirement: str
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
            "projection_channel": self.projection_channel,
            "lambda_form": self.lambda_form,
            "verdict": self.verdict,
            "collapse_condition": self.collapse_condition,
            "next_requirement": self.next_requirement,
            "notes": self.notes,
        }


def symbols() -> dict[str, sp.Symbol]:
    names = "z Omega tau_chi beta c_Y Gamma kappa F_hat O_hat"
    return {name: sp.Symbol(name, real=True) for name in names.split()}


def relaxation_transfer() -> sp.Expr:
    sym = symbols()
    return sym["c_Y"] + sym["beta"] / (1 + sym["tau_chi"] * sym["z"])


def acceleration_projection() -> sp.Expr:
    sym = symbols()
    return sym["Gamma"]


def range_projection() -> sp.Expr:
    sym = symbols()
    return sym["Gamma"] / (sym["kappa"] ** 2 + sym["z"] ** 2)


def observed_acceleration_transfer() -> sp.Expr:
    return sp.factor(acceleration_projection() * relaxation_transfer())


def observed_range_transfer() -> sp.Expr:
    return sp.factor(range_projection() * relaxation_transfer())


def range_deprojection_residual() -> sp.Expr:
    sym = symbols()
    observed = observed_range_transfer()
    deprojected = sp.simplify(
        observed * (sym["kappa"] ** 2 + sym["z"] ** 2) / sym["Gamma"]
    )
    return sp.factor(sp.simplify(deprojected - relaxation_transfer()))


def range_equation_residual() -> sp.Expr:
    sym = symbols()
    r_hat = observed_range_transfer() * sym["F_hat"]
    lhs = (sym["kappa"] ** 2 + sym["z"] ** 2) * r_hat
    rhs = sym["Gamma"] * relaxation_transfer() * sym["F_hat"]
    return sp.factor(sp.simplify(lhs - rhs))


def range_relaxation_pole_residue() -> sp.Expr:
    sym = symbols()
    pole = -1 / sym["tau_chi"]
    transfer = observed_range_transfer()
    residue = sp.limit((sym["z"] - pole) * transfer, sym["z"], pole)
    return sp.factor(sp.simplify(residue))


def projection_pole_condition() -> sp.Expr:
    sym = symbols()
    return sym["kappa"] ** 2 + sym["z"] ** 2


def arbitrary_projection_solution() -> sp.Expr:
    sym = symbols()
    return sp.simplify(sym["O_hat"] / (relaxation_transfer() * sym["F_hat"]))


def projection_rows() -> list[ProjectionRow]:
    return [
        ProjectionRow(
            term="acceleration_constant_projection",
            claim_status="Counterexample candidate",
            projection_channel="delta a_hat = Gamma q_hat",
            lambda_form="Gamma",
            verdict="pole relation survives up to one scale",
            collapse_condition="Gamma = 0 or beta tau_chi = 0",
            next_requirement="calibrate Gamma or use ratios insensitive to scale",
            notes=(
                "A constant projection cannot create arbitrary frequency "
                "dependence. It only rescales c_Y and beta."
            ),
        ),
        ProjectionRow(
            term="range_known_projection",
            claim_status="Proven",
            projection_channel="ddot R + kappa^2 R = Gamma q_A",
            lambda_form="Gamma/(kappa^2 + z^2)",
            verdict="exactly deprojects to G(z) off projection poles",
            collapse_condition="Gamma = 0 or kappa^2 + z_k^2 = 0",
            next_requirement="known Gamma and kappa",
            notes=(
                "Multiplying by (kappa^2+z^2)/Gamma recovers the same "
                "relaxation transfer."
            ),
        ),
        ProjectionRow(
            term="range_shared_finite_nuisance",
            claim_status="Counterexample candidate",
            projection_channel="shared range response",
            lambda_form="Gamma/(kappa^2 + z^2) with shared Gamma,kappa",
            verdict="finite nuisance, not arbitrary per-frequency nuisance",
            collapse_condition="too few samples to fit both pole and projection nuisance",
            next_requirement="independent calibration or higher-frequency count",
            notes=(
                "The relaxation pole survives, but exact distinguishability "
                "must include the shared projection parameters."
            ),
        ),
        ProjectionRow(
            term="range_relaxation_pole_survival",
            claim_status="Proven",
            projection_channel="range response",
            lambda_form="Gamma/(kappa^2 + z^2)",
            verdict="relaxation pole remains in observed transfer",
            collapse_condition="Gamma beta = 0 or projection zero at relaxation pole",
            next_requirement="avoid pole-cancelling projection",
            notes=(
                "At z=-1/tau_chi the numerator of the range transfer is "
                "Gamma beta."
            ),
        ),
        ProjectionRow(
            term="arbitrary_frequency_projection",
            claim_status="Proven",
            projection_channel="O_hat(z_k)=Lambda_k G(z_k) F_hat(z_k)",
            lambda_form="unconstrained complex Lambda_k",
            verdict="projection degeneracy",
            collapse_condition="Lambda_k free at every sampled frequency",
            next_requirement="finite-dimensional or calibrated projection model",
            notes=(
                "Choosing Lambda_k=O_hat_k/(G_k F_hat_k) absorbs the pole "
                "point by point."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in projection_rows()]


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
            "body_response": "q_hat(z) = G(z) F_hat(z)",
            "transfer": "G(z) = c_Y + beta/(1 + tau_chi z)",
            "acceleration_channel": "delta a_hat = Gamma q_hat",
            "range_channel": "ddot R + kappa^2 R = Gamma q_A",
            "scope": "symbolic projection audit; no empirical data channel",
        },
        "relaxation_transfer": sp.sstr(relaxation_transfer()),
        "acceleration_projection": sp.sstr(acceleration_projection()),
        "range_projection": sp.sstr(range_projection()),
        "observed_acceleration_transfer": sp.sstr(observed_acceleration_transfer()),
        "observed_range_transfer": sp.sstr(observed_range_transfer()),
        "range_deprojection_residual": sp.sstr(range_deprojection_residual()),
        "range_equation_residual": sp.sstr(range_equation_residual()),
        "range_relaxation_pole_residue": sp.sstr(range_relaxation_pole_residue()),
        "projection_pole_condition": sp.sstr(projection_pole_condition()),
        "arbitrary_projection_solution": sp.sstr(arbitrary_projection_solution()),
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
    print(f"range_deprojection_residual: {data['range_deprojection_residual']}")
    print(f"range_relaxation_pole_residue: {data['range_relaxation_pole_residue']}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
