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
    / "dynamic_chi_named_timing_model_projection_audit.tsv"
)
JSON_PATH = (
    REPO_ROOT
    / "outputs"
    / "json"
    / "dynamic_chi_named_timing_model_projection_audit.json"
)

TABLE_COLUMNS = (
    "term",
    "claim_status",
    "named_model_component",
    "source_anchor",
    "implementation_evidence",
    "projection_nuisance_class",
    "inferred_rank_class",
    "bridge_verdict",
    "runtime_worthiness",
    "notes",
)


@dataclass(frozen=True)
class NamedTimingModelRow:
    term: str
    claim_status: str
    named_model_component: str
    source_anchor: str
    implementation_evidence: str
    projection_nuisance_class: str
    inferred_rank_class: str
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


def source_release_ledger() -> dict[str, str]:
    return {
        "2025_zenodo": "https://zenodo.org/records/13899771",
        "2025_paper": (
            "https://www.aanda.org/component/article?access=doi&doi="
            "10.1051%2F0004-6361%2F202452100"
        ),
        "2020_zenodo": "https://zenodo.org/records/3778978",
        "local_source_scope": (
            "Nutimo 2025 source inspected from nutimo.tar.bz2; no build, "
            "runtime, TOA fit, or dependency provisioning was attempted."
        ),
    }


def carrier_set() -> str:
    return "Omega_in, Omega_out, abs(Omega_in-Omega_out)"


def named_model_rows() -> list[NamedTimingModelRow]:
    carriers = carrier_set()
    return [
        NamedTimingModelRow(
            term="public_named_model_source_path",
            claim_status="Proven",
            named_model_component="Nutimo/J0337 public release path",
            source_anchor="Request 5 notes; Zenodo 2025 record 13899771",
            implementation_evidence=(
                "The public release path contains Nutimo source and J0337 data "
                "products; this audit inspects the source path only."
            ),
            projection_nuisance_class="source-inspection boundary",
            inferred_rank_class="not a projection model",
            bridge_verdict="conditional",
            runtime_worthiness="conditional",
            notes=(
                "Public source availability is imported as a setup fact. The "
                "runtime stack remains intentionally out of scope here."
            ),
        ),
        NamedTimingModelRow(
            term="finite_fittable_parameter_map",
            claim_status="Proven",
            named_model_component="Parametres fittable parameter list",
            source_anchor="Parameters.h:63-122",
            implementation_evidence=(
                "The core parameter map is a finite list of spin, inner-orbit, "
                "outer-orbit, mass, orientation, SEP, sky, and extra-body "
                "parameters rather than Fourier amplitudes per carrier."
            ),
            projection_nuisance_class="finite shared physical manifold",
            inferred_rank_class="finite parameter rank, not per-carrier complex",
            bridge_verdict="conditional",
            runtime_worthiness="runtime-motivated",
            notes=(
                f"For the {carriers} bridge this supports the finite-manifold "
                "side of the 10.5 gate, pending actual fit configuration."
            ),
        ),
        NamedTimingModelRow(
            term="motion_parameter_recomputes_state_vectors",
            claim_status="Proven",
            named_model_component="parameter perturbation and state recomputation",
            source_anchor="Parameters.cpp:2216-2278, 2418-2455",
            implementation_evidence=(
                "Changes to orbital, SEP, and extra-body parameters set "
                "motion_changed and route through Compute_state_vectors."
            ),
            projection_nuisance_class="finite state-vector projection",
            inferred_rank_class="shared nonlinear geometry rank",
            bridge_verdict="conditional",
            runtime_worthiness="runtime-motivated",
            notes=(
                "Carrier amplitudes are generated by recomputed worldline and "
                "delay maps, not by independent complex coefficients at each "
                "frequency."
            ),
        ),
        NamedTimingModelRow(
            term="delay_projection_from_shared_states",
            claim_status="Proven",
            named_model_component="Fittriple delay computation",
            source_anchor="Fittriple-compute.cpp:167-260; Delay_brut.cpp:512-735",
            implementation_evidence=(
                "Einstein, Shapiro, aberration, and geometric delay calls take "
                "the shared pulsar/companion state vectors sp, si, so and "
                "line-of-sight geometry as inputs."
            ),
            projection_nuisance_class="finite shared physical manifold",
            inferred_rank_class="shared state and geometry rank",
            bridge_verdict="conditional",
            runtime_worthiness="runtime-motivated",
            notes=(
                "This is the named-code analog of the phase/geometric-lock "
                "assumption used in the symbolic 10.5 gate."
            ),
        ),
        NamedTimingModelRow(
            term="geometric_delay_line_of_sight_projection",
            claim_status="Proven",
            named_model_component="geometric/Roemer projection",
            source_anchor="Delay_brut.cpp:715-809",
            implementation_evidence=(
                "The geometric delay projects the pulsar state against the "
                "SSB-to-pulsar-system line of sight and optional Kopeikin/"
                "Shklovskii terms, with component delay details available."
            ),
            projection_nuisance_class="finite shared geometry projection",
            inferred_rank_class="geometry-locked finite rank",
            bridge_verdict="conditional",
            runtime_worthiness="runtime-motivated",
            notes=(
                "The projection is parameterized by physical geometry; it is "
                "not written as carrier-local complex amplitude freedom."
            ),
        ),
        NamedTimingModelRow(
            term="finite_parameter_gradient_program",
            claim_status="Proven",
            named_model_component="Compute_gradients design matrix path",
            source_anchor="Compute_gradients.cpp:54-199",
            implementation_evidence=(
                "The gradient utility differentiates each residual with respect "
                "to the finite fitted parameter list and writes an ntoa by npar "
                "Jacobian."
            ),
            projection_nuisance_class="finite fitted-parameter Jacobian",
            inferred_rank_class="finite model derivative rank",
            bridge_verdict="conditional",
            runtime_worthiness="runtime-motivated",
            notes=(
                "A future runtime hand-off can test the exact rank class by "
                "adding a dynamic-chi column to this finite Jacobian."
            ),
        ),
        NamedTimingModelRow(
            term="spin_phase_and_mean_offset_nuisance",
            claim_status="Proven",
            named_model_component="phase residual and mean removal",
            source_anchor="Fittriple-compute.cpp:431-453",
            implementation_evidence=(
                "Residuals are computed from spin phase, spin derivatives, "
                "dphase0, turn counts, and an optional weighted mean removal."
            ),
            projection_nuisance_class="finite spin and offset nuisance",
            inferred_rank_class="absorbs secular or constant pieces only",
            bridge_verdict="conditional",
            runtime_worthiness="conditional",
            notes=(
                "This nuisance is expected and does not by itself supply "
                "independent complex freedom at the three carrier frequencies."
            ),
        ),
        NamedTimingModelRow(
            term="rn_pl_specialcase_harmonic_nuisance",
            claim_status="Proven",
            named_model_component="RN_PL special case",
            source_anchor="Fittriple-compute.cpp:343-383",
            implementation_evidence=(
                "The RN_PL branch introduces explicit sinusoid amplitudes and "
                "phases ds_Ak, ds_phik for multiple harmonics."
            ),
            projection_nuisance_class="effectively per-harmonic complex nuisance",
            inferred_rank_class="can reach per-carrier complex rank if enabled",
            bridge_verdict="collapse",
            runtime_worthiness="not-runtime-motivated",
            notes=(
                "Dynamic-chi tests must disable or explicitly isolate this "
                "special case; otherwise it can absorb carrier-level phase and "
                "amplitude relations by construction."
            ),
        ),
        NamedTimingModelRow(
            term="named_model_projection_class_verdict",
            claim_status="Counterexample candidate",
            named_model_component="standard Nutimo triple timing core",
            source_anchor=(
                "Combined Parameters.h, Parameters.cpp, Fittriple-compute.cpp, "
                "Delay_brut.cpp, Compute_gradients.cpp inspection"
            ),
            implementation_evidence=(
                "The standard core is finite parameter/state/geometry based; "
                "the explicit harmonic special case is a separable collapse "
                "risk rather than the default projection form."
            ),
            projection_nuisance_class="hybrid but core is finite shared physical manifold",
            inferred_rank_class=(
                "conditional finite shared rank; collapses if arbitrary "
                "harmonic nuisance is admitted on target carriers"
            ),
            bridge_verdict="conditional",
            runtime_worthiness="runtime-motivated",
            notes=(
                "The named source inspection is positive enough to justify an "
                "external runtime hand-off only for a configuration that keeps "
                "the finite geometric model and excludes per-harmonic carrier "
                "soak-up terms."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in named_model_rows()]


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


def verdict_summary(rows: list[dict[str, str]] | None = None) -> dict[str, str]:
    rows = rows or rows_as_dicts()
    by_term = {row["term"]: row for row in rows}
    core = by_term["named_model_projection_class_verdict"]
    collapse = by_term["rn_pl_specialcase_harmonic_nuisance"]
    return {
        "standard_core_bridge_verdict": core["bridge_verdict"],
        "standard_core_runtime_worthiness": core["runtime_worthiness"],
        "specialcase_bridge_verdict": collapse["bridge_verdict"],
        "specialcase_runtime_worthiness": collapse["runtime_worthiness"],
        "implementation_class": (
            "conditional positive: standard core is finite shared; "
            "arbitrary harmonic nuisance collapses the bridge if enabled"
        ),
    }


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "scope": "named source-inspection audit only; no runtime",
            "carrier_set": carrier_set(),
            "dynamic_chi_target": "shared finite-tau_chi transfer law G(z)",
            "gate_imported_from_request_10_5": (
                "finite shared phase/geometric projection is runtime-motivated; "
                "independent per-carrier complex projection collapses"
            ),
        },
        "source_releases": source_release_ledger(),
        "verdict_summary": verdict_summary(rows),
        "runtime_handoff_rule": {
            "runtime-motivated": (
                "Use the standard integrated three-body finite-parameter core "
                "and test a dynamic-chi column against its Jacobian."
            ),
            "not-runtime-motivated": (
                "Do not claim a dynamic-chi observable if arbitrary harmonic "
                "amplitude/phase nuisance is enabled on the target carriers."
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
    print(f"verdict_summary: {data['verdict_summary']}")
    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
