from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from frequency_sweep_distinguishability import ALLOWED_CLAIM_STATUSES


REPO_ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = (
    REPO_ROOT
    / "outputs"
    / "tsv"
    / "dynamic_chi_external_nutimo_handoff_packet.tsv"
)
JSON_PATH = (
    REPO_ROOT
    / "outputs"
    / "json"
    / "dynamic_chi_external_nutimo_handoff_packet.json"
)

TABLE_COLUMNS = (
    "packet_item",
    "claim_status",
    "pilot_stage",
    "external_owner_action",
    "required_input",
    "requested_artifact",
    "validation_rule",
    "pass_interpretation",
    "fail_interpretation",
    "notes",
)


@dataclass(frozen=True)
class HandoffPacketRow:
    packet_item: str
    claim_status: str
    pilot_stage: str
    external_owner_action: str
    required_input: str
    requested_artifact: str
    validation_rule: str
    pass_interpretation: str
    fail_interpretation: str
    notes: str

    def __post_init__(self) -> None:
        if self.claim_status not in ALLOWED_CLAIM_STATUSES:
            allowed = ", ".join(sorted(ALLOWED_CLAIM_STATUSES))
            raise ValueError(
                f"claim_status must be one of {allowed}; got {self.claim_status!r}"
            )

    def to_row(self) -> dict[str, str]:
        return asdict(self)


def target_carriers() -> list[str]:
    return ["Omega_in", "Omega_out", "abs(Omega_in-Omega_out)"]


def expected_return_files() -> dict[str, str]:
    return {
        "configuration_manifest.json": (
            "active Nutimo configuration, specialcase flag, delay flags, "
            "fitted parameter names, and target carrier extraction settings"
        ),
        "baseline_fit_summary.json": (
            "baseline parameter vector, timing span, residual RMS or chi2 "
            "summary, and code/data release identifiers"
        ),
        "finite_jacobian.npy_or_tsv": (
            "ntoa by nfit residual derivative matrix for the standard fitted "
            "parameter set"
        ),
        "carrier_projection_rank.json": (
            "six-real-coordinate three-carrier projection rank, singular values, "
            "and rank threshold"
        ),
        "dynamic_chi_test_column.tsv": (
            "synthetic dynamic-chi residual column and carrier-level complex "
            "amplitudes for the shared-tau transfer"
        ),
        "synthetic_injection_recovery.json": (
            "refit result for the minimal shared-tau injection, including "
            "pre/post carrier amplitudes and residual projection onto the "
            "finite nuisance span"
        ),
        "decision_summary.md": (
            "one-page verdict using the Request 10.7 pass/fail language"
        ),
    }


def handoff_rows() -> list[HandoffPacketRow]:
    carriers = ", ".join(target_carriers())
    return [
        HandoffPacketRow(
            packet_item="environment_boundary",
            claim_status="Conjectural",
            pilot_stage="preflight",
            external_owner_action=(
                "Declare the pre-provisioned Nutimo/Tempo2/Boost/Minuit "
                "environment; do not report dependency repair as science."
            ),
            required_input="public Nutimo source/data release and runnable environment",
            requested_artifact="baseline_fit_summary.json",
            validation_rule=(
                "must identify code release, data release, executable path, "
                "and whether residuals plus finite gradients can be computed"
            ),
            pass_interpretation="runtime gate can start",
            fail_interpretation="handoff stops; no implementation-class evidence",
            notes="This keeps the pilot from becoming an infrastructure project.",
        ),
        HandoffPacketRow(
            packet_item="configuration_manifest",
            claim_status="Counterexample candidate",
            pilot_stage="configuration_closure",
            external_owner_action=(
                "Export active fit parameters, delay flags, specialcase flag, "
                "target carrier list, and any harmonic nuisance blocks."
            ),
            required_input="baseline par/tim and active Nutimo run configuration",
            requested_artifact="configuration_manifest.json",
            validation_rule=(
                "RN_PL-like or equivalent harmonic soak-up must be explicitly "
                "listed as enabled or disabled for the target carriers"
            ),
            pass_interpretation="standard finite state/geometry core is testable",
            fail_interpretation="collapse if target-carrier harmonic soak-up is enabled",
            notes="This is the first stop rule from Request 10.7.",
        ),
        HandoffPacketRow(
            packet_item="target_carrier_inventory",
            claim_status="Counterexample candidate",
            pilot_stage="configuration_closure",
            external_owner_action=(
                f"Extract the three target carrier frequencies: {carriers}."
            ),
            required_input="inner/outer periods and outer-dipole combination component",
            requested_artifact="configuration_manifest.json",
            validation_rule=(
                "frequencies must be positive, distinct, nonresonant at the "
                "2:1 collapse cases, and nonzero in projection"
            ),
            pass_interpretation="three-carrier shared-tau bridge remains available",
            fail_interpretation="falls back to weaker two-carrier or one-carrier gate",
            notes="This is not a new harmonic search; it uses existing GR carriers.",
        ),
        HandoffPacketRow(
            packet_item="finite_parameter_jacobian",
            claim_status="Counterexample candidate",
            pilot_stage="named_jacobian_rank_gate",
            external_owner_action=(
                "Compute residual derivatives with respect to the standard "
                "finite fitted parameter set near the baseline solution."
            ),
            required_input="runnable residual function and fitted parameter scales",
            requested_artifact="finite_jacobian.npy_or_tsv",
            validation_rule=(
                "columns must correspond to named fitted parameters, not "
                "carrier-local Fourier amplitudes"
            ),
            pass_interpretation="finite nuisance span is available for rank test",
            fail_interpretation="cannot classify named runtime projection",
            notes="This is the minimal useful runtime calculation before posterior work.",
        ),
        HandoffPacketRow(
            packet_item="carrier_projection_rank",
            claim_status="Counterexample candidate",
            pilot_stage="named_jacobian_rank_gate",
            external_owner_action=(
                "Project the finite nuisance Jacobian onto the six real "
                "coordinates of the three carrier complex amplitudes."
            ),
            required_input="finite_jacobian.npy_or_tsv and carrier extraction operator",
            requested_artifact="carrier_projection_rank.json",
            validation_rule=(
                "rank <= 5 supports finite shared phase/geometric manifold; "
                "rank 6/6 is the per-carrier complex collapse class"
            ),
            pass_interpretation="dynamic-chi remains runtime-motivated if test column also survives",
            fail_interpretation="named implementation behaves like carrier-local soak-up",
            notes="Report singular values and rank threshold, not just an integer.",
        ),
        HandoffPacketRow(
            packet_item="dynamic_chi_test_column",
            claim_status="Counterexample candidate",
            pilot_stage="named_jacobian_rank_gate",
            external_owner_action=(
                "Construct the synthetic shared-tau residual column using "
                "G(z)=c_Y+beta/(1+tau_chi z) on the target carriers."
            ),
            required_input="chosen tau_chi, beta, c_Y, carrier amplitudes, and projection convention",
            requested_artifact="dynamic_chi_test_column.tsv",
            validation_rule=(
                "column must be tested for membership in the finite nuisance "
                "Jacobian span"
            ),
            pass_interpretation="proceed to minimal synthetic injection",
            fail_interpretation="dynamic-chi collapses into named finite nuisance span",
            notes="Use several finite tau_chi values away from zero-frequency collapse.",
        ),
        HandoffPacketRow(
            packet_item="minimal_synthetic_injection",
            claim_status="Counterexample candidate",
            pilot_stage="synthetic_injection",
            external_owner_action=(
                "Inject the shared-tau carrier perturbation and refit only the "
                "standard finite nuisance parameters."
            ),
            required_input="baseline residuals, finite Jacobian or fit routine, dynamic_chi_test_column",
            requested_artifact="synthetic_injection_recovery.json",
            validation_rule=(
                "post-refit carrier relation must be tested against one shared "
                "tau_chi and against the admitted finite derivative comparator"
            ),
            pass_interpretation="promote to external runtime experiment",
            fail_interpretation="record named implementation conditional no-go",
            notes="This is a sanity pilot, not a detection or final posterior.",
        ),
        HandoffPacketRow(
            packet_item="decision_summary",
            claim_status="Conjectural",
            pilot_stage="decision",
            external_owner_action=(
                "Return a one-page verdict using the exact labels: "
                "runtime-motivated, conditional, collapse."
            ),
            required_input="all prior handoff artifacts",
            requested_artifact="decision_summary.md",
            validation_rule=(
                "must state which stop rule triggered, or why external runtime "
                "work is promoted"
            ),
            pass_interpretation="dynamic-chi is a named runtime-worthy positive branch",
            fail_interpretation="dynamic-chi closes for this named implementation class",
            notes="Do not convert a failed pilot into an open-ended source chase.",
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in handoff_rows()]


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
        "scope": (
            "external handoff packet for Request 10.7; this repo does not "
            "perform the Nutimo runtime run"
        ),
        "target_carriers": target_carriers(),
        "expected_return_files": expected_return_files(),
        "hard_stop_rules": [
            "dependency repair replaces the science gate",
            "RN_PL-like harmonic soak-up is enabled on target carriers",
            "target carrier inventory collapses below the required sample count",
            "carrier projection rank is 6/6",
            "dynamic-chi test column lies in the finite nuisance span",
        ],
        "promotion_rule": (
            "promote only if configuration closure passes, projection rank is "
            "<=5, dynamic-chi column survives the finite nuisance span, and "
            "minimal synthetic shared-tau injection is recovered after refit"
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
    print(f"target_carriers: {', '.join(data['target_carriers'])}")
    print(f"expected_return_files: {len(data['expected_return_files'])}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
