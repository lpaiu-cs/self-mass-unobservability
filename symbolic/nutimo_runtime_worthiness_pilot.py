from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

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
    / "dynamic_chi_nutimo_runtime_worthiness_pilot.tsv"
)
JSON_PATH = (
    REPO_ROOT
    / "outputs"
    / "json"
    / "dynamic_chi_nutimo_runtime_worthiness_pilot.json"
)

TABLE_COLUMNS = (
    "gate",
    "claim_status",
    "pilot_stage",
    "required_inputs",
    "pass_condition",
    "fail_condition",
    "projection_class_if_pass",
    "projection_class_if_fail",
    "bridge_verdict_if_pass",
    "bridge_verdict_if_fail",
    "runtime_action",
    "notes",
)

ALLOWED_GATE_ACTIONS = {
    "prepare_external_environment_only_if_needed",
    "proceed_to_jacobian_gate",
    "stop_not_runtime_motivated",
    "proceed_to_synthetic_injection",
    "promote_to_external_runtime_experiment",
    "record_conditional_no_go",
}


@dataclass(frozen=True)
class RuntimePilotGate:
    gate: str
    claim_status: str
    pilot_stage: str
    required_inputs: str
    pass_condition: str
    fail_condition: str
    projection_class_if_pass: str
    projection_class_if_fail: str
    bridge_verdict_if_pass: str
    bridge_verdict_if_fail: str
    runtime_action: str
    notes: str

    def __post_init__(self) -> None:
        if self.claim_status not in ALLOWED_CLAIM_STATUSES:
            allowed = ", ".join(sorted(ALLOWED_CLAIM_STATUSES))
            raise ValueError(
                f"claim_status must be one of {allowed}; got {self.claim_status!r}"
            )
        for field_name in ("bridge_verdict_if_pass", "bridge_verdict_if_fail"):
            value = getattr(self, field_name)
            if value not in ALLOWED_BRIDGE_VERDICTS:
                allowed = ", ".join(sorted(ALLOWED_BRIDGE_VERDICTS))
                raise ValueError(f"{field_name} must be one of {allowed}; got {value!r}")
        if self.runtime_action not in ALLOWED_GATE_ACTIONS:
            allowed = ", ".join(sorted(ALLOWED_GATE_ACTIONS))
            raise ValueError(
                f"runtime_action must be one of {allowed}; got {self.runtime_action!r}"
            )

    def to_row(self) -> dict[str, str]:
        return asdict(self)


def carrier_set() -> str:
    return "Omega_in, Omega_out, abs(Omega_in-Omega_out)"


def pilot_stage_order() -> list[str]:
    return [
        "source_and_environment_boundary",
        "configuration_closure",
        "named_jacobian_rank_gate",
        "minimal_synthetic_shared_tau_injection",
        "runtime_promotion_decision",
    ]


def jacobian_rank_gate_verdict(
    projection_rank: int,
    dynamic_chi_column_in_span: bool,
) -> dict[str, str]:
    if projection_rank >= 6 or dynamic_chi_column_in_span:
        return {
            "bridge_verdict": "collapse",
            "runtime_worthiness": "not-runtime-motivated",
            "classification": (
                "effective per-carrier complex freedom or dynamic-chi column "
                "absorbed by the named nuisance span"
            ),
        }
    return {
        "bridge_verdict": "conditional",
        "runtime_worthiness": "runtime-motivated",
        "classification": (
            "finite shared projection rank with an unabsorbed dynamic-chi "
            "test column"
        ),
    }


def minimal_external_artifacts() -> dict[str, str]:
    return {
        "configuration_manifest": (
            "active parameter blocks, specialcase flag, delay flags, target "
            "carrier list, and whether RN_PL-like harmonic soak-up is disabled"
        ),
        "finite_jacobian": (
            "residual derivative matrix for the standard fitted parameter set "
            "near the baseline solution"
        ),
        "carrier_projection_rank_summary": (
            "effective rank of the three-carrier projection nuisance manifold; "
            "rank <= 5 keeps the 10.5 phase/geometric lock gate alive"
        ),
        "dynamic_chi_column": (
            "synthetic shared-tau carrier response column evaluated at "
            "Omega_in, Omega_out, and abs(Omega_in-Omega_out)"
        ),
        "injection_recovery_table": (
            "carrier amplitudes/phases before and after refit, with the shared "
            "G(z)=c_Y+beta/(1+tau_chi z) relation tested across all carriers"
        ),
    }


def pilot_rows() -> list[RuntimePilotGate]:
    carriers = carrier_set()
    return [
        RuntimePilotGate(
            gate="external_environment_boundary",
            claim_status="Conjectural",
            pilot_stage="source_and_environment_boundary",
            required_inputs=(
                "pre-provisioned Nutimo/Tempo2/Boost/Minuit environment or "
                "external collaborator runtime"
            ),
            pass_condition=(
                "environment can compute residuals and finite-parameter "
                "gradients without changing the scientific comparator class"
            ),
            fail_condition=(
                "local dependency repair becomes open-ended or changes the "
                "timing model configuration"
            ),
            projection_class_if_pass="not yet classified",
            projection_class_if_fail="not applicable",
            bridge_verdict_if_pass="conditional",
            bridge_verdict_if_fail="collapse",
            runtime_action="prepare_external_environment_only_if_needed",
            notes=(
                "This gate authorizes runtime only as a science test. It is not "
                "a request to reopen local build infrastructure work."
            ),
        ),
        RuntimePilotGate(
            gate="configuration_closure_standard_core",
            claim_status="Counterexample candidate",
            pilot_stage="configuration_closure",
            required_inputs=(
                "active par/tim configuration, delay flags, fitted parameter "
                "map, specialcase flag, and target carrier list"
            ),
            pass_condition=(
                "standard integrated three-body core is active and explicit "
                "harmonic soak-up terms are disabled for target carriers"
            ),
            fail_condition=(
                "RN_PL-like per-harmonic amplitude/phase nuisance or "
                "equivalent carrier-local soak-up is enabled"
            ),
            projection_class_if_pass="finite shared state/geometry manifold",
            projection_class_if_fail="effectively per-harmonic complex nuisance",
            bridge_verdict_if_pass="conditional",
            bridge_verdict_if_fail="collapse",
            runtime_action="proceed_to_jacobian_gate",
            notes=(
                "This is the runtime version of the 10.6 source-inspection "
                "boundary."
            ),
        ),
        RuntimePilotGate(
            gate="carrier_inventory_closure",
            claim_status="Counterexample candidate",
            pilot_stage="configuration_closure",
            required_inputs=(
                "nonzero fitted projections for "
                f"{carriers}; inner and outer periods from the named model"
            ),
            pass_condition=(
                "three distinct positive carrier samples are present and not "
                "lost to 2:1 resonance, zero projection, or missing component"
            ),
            fail_condition=(
                "carrier inventory collapses to one or two usable samples"
            ),
            projection_class_if_pass="three-sample finite carrier inventory",
            projection_class_if_fail="insufficient carrier inventory",
            bridge_verdict_if_pass="conditional",
            bridge_verdict_if_fail="collapse",
            runtime_action="proceed_to_jacobian_gate",
            notes=(
                "Two carriers only test low-order real derivative comparators; "
                "the named pilot targets the stronger three-carrier bridge."
            ),
        ),
        RuntimePilotGate(
            gate="named_jacobian_projection_rank",
            claim_status="Counterexample candidate",
            pilot_stage="named_jacobian_rank_gate",
            required_inputs=(
                "finite residual Jacobian for standard fitted parameters and "
                "a carrier-projection extraction around the baseline fit"
            ),
            pass_condition=(
                "effective projection nuisance rank is <=5 in the six-real "
                "three-carrier amplitude space and the dynamic-chi test column "
                "is not in the fitted nuisance span"
            ),
            fail_condition=(
                "effective projection rank reaches 6/6 or the dynamic-chi "
                "column is absorbed by the fitted parameter span"
            ),
            projection_class_if_pass="finite shared phase/geometric manifold",
            projection_class_if_fail="effectively per-carrier complex projection",
            bridge_verdict_if_pass="conditional",
            bridge_verdict_if_fail="collapse",
            runtime_action="proceed_to_synthetic_injection",
            notes=(
                "This is the decisive runtime-worthiness gate. It should be run "
                "before any posterior or long MCMC experiment."
            ),
        ),
        RuntimePilotGate(
            gate="minimal_synthetic_shared_tau_injection",
            claim_status="Counterexample candidate",
            pilot_stage="minimal_synthetic_shared_tau_injection",
            required_inputs=(
                "synthetic carrier perturbation generated by "
                "G(z)=c_Y+beta/(1+tau_chi z) on the three target carriers"
            ),
            pass_condition=(
                "after refitting standard finite parameters, the remaining "
                "carrier relation is still best described by one shared "
                "tau_chi rather than by the admitted finite derivative "
                "comparator"
            ),
            fail_condition=(
                "standard finite parameters or admitted harmonic nuisance "
                "remove the shared-tau residual across the target carriers"
            ),
            projection_class_if_pass="runtime-visible shared-tau relation",
            projection_class_if_fail="named-model nuisance collapse",
            bridge_verdict_if_pass="distinguishable",
            bridge_verdict_if_fail="collapse",
            runtime_action="promote_to_external_runtime_experiment",
            notes=(
                "This is a synthetic sanity pilot only. It is not a final "
                "J0337 constraint or detection claim."
            ),
        ),
        RuntimePilotGate(
            gate="runtime_promotion_decision",
            claim_status="Conjectural",
            pilot_stage="runtime_promotion_decision",
            required_inputs=(
                "configuration manifest, rank-gate result, and synthetic "
                "injection recovery table"
            ),
            pass_condition=(
                "standard core is finite shared, rank gate is <=5, and shared "
                "tau_chi injection survives the finite nuisance refit"
            ),
            fail_condition=(
                "any earlier gate collapses or requires arbitrary per-carrier "
                "complex nuisance to describe the named fit"
            ),
            projection_class_if_pass="conditional positive named runtime branch",
            projection_class_if_fail="named implementation no-go or conditional no-go",
            bridge_verdict_if_pass="distinguishable",
            bridge_verdict_if_fail="collapse",
            runtime_action="record_conditional_no_go",
            notes=(
                "A pass promotes external runtime work to a real observable "
                "experiment; a fail closes the current dynamic-chi branch for "
                "this named implementation class."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in pilot_rows()]


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
        if row["bridge_verdict_if_pass"] not in ALLOWED_BRIDGE_VERDICTS:
            raise ValueError(f"row {index} has invalid bridge_verdict_if_pass")
        if row["bridge_verdict_if_fail"] not in ALLOWED_BRIDGE_VERDICTS:
            raise ValueError(f"row {index} has invalid bridge_verdict_if_fail")
        if row["runtime_action"] not in ALLOWED_GATE_ACTIONS:
            raise ValueError(f"row {index} has invalid runtime_action")
        if not all(str(row[column]).strip() for column in TABLE_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "scope": (
                "external runtime-worthiness pilot contract; no local Nutimo "
                "build or TOA fit is performed by this script"
            ),
            "carrier_set": carrier_set(),
            "dynamic_transfer": "G(z)=c_Y+beta/(1+tau_chi z)",
            "named_model_boundary": (
                "standard Nutimo core is finite state/geometry unless the "
                "fit enables explicit harmonic soak-up on target carriers"
            ),
        },
        "stage_order": pilot_stage_order(),
        "minimal_external_artifacts": minimal_external_artifacts(),
        "rank_gate_examples": {
            "rank5_unabsorbed": jacobian_rank_gate_verdict(
                projection_rank=5,
                dynamic_chi_column_in_span=False,
            ),
            "rank6_collapse": jacobian_rank_gate_verdict(
                projection_rank=6,
                dynamic_chi_column_in_span=False,
            ),
            "rank5_but_column_absorbed": jacobian_rank_gate_verdict(
                projection_rank=5,
                dynamic_chi_column_in_span=True,
            ),
        },
        "pilot_verdict": (
            "conditional runtime pilot is scientifically motivated, but only "
            "after configuration closure and before any long posterior run"
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
    print(f"pilot_verdict: {data['pilot_verdict']}")
    print(f"rank_gate_examples: {data['rank_gate_examples']}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
