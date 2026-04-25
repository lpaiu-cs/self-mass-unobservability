from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp

from frequency_sweep_distinguishability import ALLOWED_CLAIM_STATUSES
from triple_projection_nuisance_gate import (
    ALLOWED_BRIDGE_VERDICTS,
    ALLOWED_RUNTIME_RELEVANCE,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = (
    REPO_ROOT
    / "outputs"
    / "tsv"
    / "dynamic_chi_triple_projection_manifold_gate.tsv"
)
JSON_PATH = (
    REPO_ROOT
    / "outputs"
    / "json"
    / "dynamic_chi_triple_projection_manifold_gate.json"
)

TABLE_COLUMNS = (
    "term",
    "claim_status",
    "carrier_set",
    "projection_manifold",
    "parameterization",
    "real_parameter_count",
    "generic_jacobian_rank",
    "nuisance_class",
    "bridge_verdict",
    "runtime_worthiness",
    "notes",
)


@dataclass(frozen=True)
class ProjectionManifoldRow:
    term: str
    claim_status: str
    carrier_set: str
    projection_manifold: str
    parameterization: str
    real_parameter_count: str
    generic_jacobian_rank: str
    nuisance_class: str
    bridge_verdict: str
    runtime_worthiness: str
    notes: str

    def __post_init__(self) -> None:
        if self.claim_status not in ALLOWED_CLAIM_STATUSES:
            allowed = ", ".join(sorted(ALLOWED_CLAIM_STATUSES))
            raise ValueError(
                f"claim_status must be one of {allowed}; got {self.claim_status!r}"
            )
        if self.bridge_verdict not in ALLOWED_BRIDGE_VERDICTS:
            allowed = ", ".join(sorted(ALLOWED_BRIDGE_VERDICTS))
            raise ValueError(
                f"bridge_verdict must be one of {allowed}; "
                f"got {self.bridge_verdict!r}"
            )
        if self.runtime_worthiness not in ALLOWED_RUNTIME_RELEVANCE:
            allowed = ", ".join(sorted(ALLOWED_RUNTIME_RELEVANCE))
            raise ValueError(
                f"runtime_worthiness must be one of {allowed}; "
                f"got {self.runtime_worthiness!r}"
            )

    def to_row(self) -> dict[str, str]:
        return asdict(self)


def symbols() -> dict[str, sp.Symbol]:
    names = "A_in A_out A_c H phi_in phi_out sigma"
    sym = {name: sp.Symbol(name, real=True) for name in names.split()}
    for name in ("x_in", "y_in", "x_out", "y_out", "x_c", "y_c"):
        sym[name] = sp.Symbol(name, real=True)
    return sym


def carrier_set() -> str:
    return "Omega_in, Omega_out, abs(Omega_in-Omega_out)"


def physical_phase_locked_vector(linked_amplitude: bool = False) -> tuple[list[sp.Expr], list[sp.Symbol]]:
    sym = symbols()
    a_in = sym["A_in"]
    a_out = sym["A_out"]
    a_c = sym["H"] * a_in * a_out if linked_amplitude else sym["A_c"]
    phi_in = sym["phi_in"]
    phi_out = sym["phi_out"]
    sigma = sym["sigma"]
    phi_c = phi_in - phi_out + sigma
    outputs = [
        a_in * sp.cos(phi_in),
        a_in * sp.sin(phi_in),
        a_out * sp.cos(phi_out),
        a_out * sp.sin(phi_out),
        a_c * sp.cos(phi_c),
        a_c * sp.sin(phi_c),
    ]
    if linked_amplitude:
        params = [sym["A_in"], sym["A_out"], sym["H"], sym["phi_in"], sym["phi_out"], sym["sigma"]]
    else:
        params = [
            sym["A_in"],
            sym["A_out"],
            sym["A_c"],
            sym["phi_in"],
            sym["phi_out"],
            sym["sigma"],
        ]
    return outputs, params


def calibrated_phase_locked_vector() -> tuple[list[sp.Expr], list[sp.Symbol]]:
    outputs, params = physical_phase_locked_vector(linked_amplitude=True)
    sym = symbols()
    return outputs, [param for param in params if param != sym["H"]]


def phase_locked_without_extra_orientation_vector() -> tuple[list[sp.Expr], list[sp.Symbol]]:
    outputs, params = physical_phase_locked_vector(linked_amplitude=False)
    sym = symbols()
    substitutions = {sym["sigma"]: 0}
    reduced_outputs = [sp.simplify(output.subs(substitutions)) for output in outputs]
    reduced_params = [param for param in params if param != sym["sigma"]]
    return reduced_outputs, reduced_params


def arbitrary_complex_vector() -> tuple[list[sp.Expr], list[sp.Symbol]]:
    sym = symbols()
    outputs = [
        sym["x_in"],
        sym["y_in"],
        sym["x_out"],
        sym["y_out"],
        sym["x_c"],
        sym["y_c"],
    ]
    params = outputs.copy()
    return outputs, params


def generic_substitutions() -> dict[sp.Symbol, sp.Expr]:
    sym = symbols()
    return {
        sym["A_in"]: sp.Integer(2),
        sym["A_out"]: sp.Integer(3),
        sym["A_c"]: sp.Integer(5),
        sym["H"]: sp.Integer(7),
        sym["phi_in"]: 0,
        sym["phi_out"]: sp.pi / 2,
        sym["sigma"]: sp.pi / 3,
        sym["x_in"]: sp.Integer(2),
        sym["y_in"]: sp.Integer(3),
        sym["x_out"]: sp.Integer(5),
        sym["y_out"]: sp.Integer(7),
        sym["x_c"]: sp.Integer(11),
        sym["y_c"]: sp.Integer(13),
    }


def jacobian_rank(outputs: list[sp.Expr], params: list[sp.Symbol]) -> int:
    matrix = sp.Matrix(outputs).jacobian(params)
    evaluated = sp.simplify(matrix.subs(generic_substitutions()))
    return int(evaluated.rank())


def phase_lock_residual() -> sp.Expr:
    sym = symbols()
    phi_c = sym["phi_in"] - sym["phi_out"]
    return sp.simplify(sp.sin(phi_c - sym["phi_in"] + sym["phi_out"]))


def rank_summary() -> dict[str, int]:
    phase_locked_outputs, phase_locked_params = phase_locked_without_extra_orientation_vector()
    linked_outputs, linked_params = calibrated_phase_locked_vector()
    arbitrary_outputs, arbitrary_params = arbitrary_complex_vector()
    return {
        "observable_vector_real_dimension": 6,
        "phase_locked_rank": jacobian_rank(phase_locked_outputs, phase_locked_params),
        "linked_amplitude_rank": jacobian_rank(linked_outputs, linked_params),
        "arbitrary_complex_rank": jacobian_rank(arbitrary_outputs, arbitrary_params),
    }


def physical_projection_rows() -> list[ProjectionManifoldRow]:
    ranks = rank_summary()
    carriers = carrier_set()
    return [
        ProjectionManifoldRow(
            term="calibrated_common_projection",
            claim_status="Counterexample candidate",
            carrier_set=carriers,
            projection_manifold="calibrated/common real projection",
            parameterization="Lambda_k=Gamma or known Lambda_k(theta)",
            real_parameter_count="0 fitted projection nuisance after calibration",
            generic_jacobian_rank="0 projection-rank contribution",
            nuisance_class="finite shared calibrated",
            bridge_verdict="distinguishable",
            runtime_worthiness="runtime-motivated",
            notes=(
                "This is the clean 10.4 positive gate: carrier ratios or "
                "deprojection recover the shared G(i Omega_k)."
            ),
        ),
        ProjectionManifoldRow(
            term="phase_locked_outer_dipole_projection",
            claim_status="Counterexample candidate",
            carrier_set=carriers,
            projection_manifold="finite shared real manifold",
            parameterization=(
                "Lambda_in=A_in e^{i phi_in}, Lambda_out=A_out e^{i phi_out}, "
                "Lambda_c=A_c e^{i(phi_in-phi_out)}"
            ),
            real_parameter_count="5 real parameters",
            generic_jacobian_rank=str(ranks["phase_locked_rank"]),
            nuisance_class="finite shared, not arbitrary per carrier",
            bridge_verdict="conditional",
            runtime_worthiness="runtime-motivated",
            notes=(
                "The combination-carrier phase is locked to the two orbital "
                "phases. The manifold has one constraint inside the six-real "
                "dimensional complex carrier vector."
            ),
        ),
        ProjectionManifoldRow(
            term="amplitude_linked_outer_dipole_projection",
            claim_status="Counterexample candidate",
            carrier_set=carriers,
            projection_manifold="finite shared real manifold with amplitude link",
            parameterization=(
                "Lambda_c=H A_in A_out e^{i(phi_in-phi_out+sigma)} with calibrated H"
            ),
            real_parameter_count="5 fitted real parameters when sigma floats; 4 if sigma is calibrated",
            generic_jacobian_rank=str(ranks["linked_amplitude_rank"]),
            nuisance_class="finite shared, geometry-linked",
            bridge_verdict="conditional",
            runtime_worthiness="runtime-motivated",
            notes=(
                "A product-like outer-dipole amplitude link is still finite and "
                "shared. It strengthens the case that the projection is not "
                "three independent complex amplitudes."
            ),
        ),
        ProjectionManifoldRow(
            term="finite_shared_complex_projection",
            claim_status="Counterexample candidate",
            carrier_set=carriers,
            projection_manifold="finite shared complex nuisance",
            parameterization="Lambda_k=ell0+ell1 z_k or equivalent shared complex model",
            real_parameter_count="4 real parameters for degree-1 complex template",
            generic_jacobian_rank="at most 4 before carrier-response parameters",
            nuisance_class="finite shared but overpowered",
            bridge_verdict="conditional",
            runtime_worthiness="conditional",
            notes=(
                "This class can consume much of the three-carrier pressure. It "
                "needs either a projection prior, more carriers, or a stricter "
                "physical derivation."
            ),
        ),
        ProjectionManifoldRow(
            term="effective_per_carrier_complex_projection",
            claim_status="Proven",
            carrier_set=carriers,
            projection_manifold="arbitrary per-carrier complex nuisance",
            parameterization=(
                "Lambda_in, Lambda_out, Lambda_c are independent complex numbers"
            ),
            real_parameter_count="6 real parameters",
            generic_jacobian_rank=str(ranks["arbitrary_complex_rank"]),
            nuisance_class="effectively arbitrary per carrier",
            bridge_verdict="collapse",
            runtime_worthiness="not-runtime-motivated",
            notes=(
                "This rank fills the full carrier-amplitude space, so the "
                "projection can absorb the shared transfer relation pointwise."
            ),
        ),
        ProjectionManifoldRow(
            term="phase_unlocked_combination_nuisance",
            claim_status="Proven",
            carrier_set=carriers,
            projection_manifold="carrier phases independently floated",
            parameterization="A_k e^{i phi_k} for k in {in,out,c}",
            real_parameter_count="6 real parameters",
            generic_jacobian_rank="6",
            nuisance_class="effectively arbitrary per carrier",
            bridge_verdict="collapse",
            runtime_worthiness="not-runtime-motivated",
            notes=(
                "If the combination phase is allowed to float independently of "
                "the inner and outer phases, the physical phase-lock is lost and "
                "the gate reduces to arbitrary complex projection."
            ),
        ),
        ProjectionManifoldRow(
            term="runtime_worthiness_gate",
            claim_status="Conjectural",
            carrier_set=carriers,
            projection_manifold="actual hierarchical-triple timing projection",
            parameterization="must be derived for the named timing model",
            real_parameter_count="model-dependent",
            generic_jacobian_rank="model-dependent",
            nuisance_class="finite shared if phase locks and geometry links are enforced",
            bridge_verdict="conditional",
            runtime_worthiness="conditional",
            notes=(
                "Runtime work becomes scientifically justified only if the "
                "actual projection is closer to the finite phase-locked manifold "
                "than to independent complex carrier amplitudes."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in physical_projection_rows()]


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
        if row["bridge_verdict"] not in ALLOWED_BRIDGE_VERDICTS:
            raise ValueError(f"row {index} has invalid bridge_verdict")
        if row["runtime_worthiness"] not in ALLOWED_RUNTIME_RELEVANCE:
            raise ValueError(f"row {index} has invalid runtime_worthiness")
        if not all(str(row[column]).strip() for column in TABLE_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "carrier_set": carrier_set(),
            "observable_vector": "three complex carrier amplitudes = six real coordinates",
            "scope": "symbolic physical projection-manifold gate; no runtime",
        },
        "rank_summary": rank_summary(),
        "phase_lock_residual": sp.sstr(phase_lock_residual()),
        "phase_lock_constraint": (
            "arg Lambda_c - arg Lambda_in + arg Lambda_out = 0 mod pi "
            "for the minimal outer-dipole phase-locked model"
        ),
        "runtime_worthiness_rule": {
            "runtime-motivated": (
                "projection is calibrated or finite shared with phase/geometric "
                "links, so runtime can test a specific shared-tau target"
            ),
            "conditional": "projection rank or priors must be closed before runtime",
            "not-runtime-motivated": (
                "projection spans independent complex carrier amplitudes and "
                "collapses the bridge"
            ),
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
    print(f"rank_summary: {data['rank_summary']}")
    print(f"phase_lock_residual: {data['phase_lock_residual']}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
