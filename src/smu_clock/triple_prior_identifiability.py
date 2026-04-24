"""Conditional identifiability audit for the leading triple clock kernel.

The rank audit showed that the decoupled outer-dipole clock kernel is a scalar
multiple of the GR carrier. This module asks what kind of external prior would
be needed to break that amplitude degeneracy in a normalized Fisher model.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import inf, sqrt
from typing import Iterable


PRIOR_COLUMNS = (
    "prior_set",
    "carrier_amplitude_degeneracy",
    "clock_amplitude_identifiability",
    "carrier_prior_precision",
    "fisher_rank",
    "marginal_clock_information",
    "condition_number",
    "verdict",
    "notes",
)


@dataclass(frozen=True)
class PriorScenario:
    """One external-prior scenario for the normalized amplitude model."""

    name: str
    carrier_prior_precision: float
    notes: str


@dataclass(frozen=True)
class PriorAuditRow:
    """One row in the prior-breaking identifiability table."""

    prior_set: str
    carrier_amplitude_degeneracy: str
    clock_amplitude_identifiability: str
    carrier_prior_precision: str
    fisher_rank: str
    marginal_clock_information: str
    condition_number: str
    verdict: str
    notes: str

    def to_row(self) -> dict[str, str]:
        return asdict(self)


def scenarios() -> tuple[PriorScenario, ...]:
    """Return the bounded prior ladder used by the audit."""

    return (
        PriorScenario(
            name="none",
            carrier_prior_precision=0.0,
            notes="No external mass, geometry, or carrier-amplitude prior.",
        ),
        PriorScenario(
            name="mass_only",
            carrier_prior_precision=0.0,
            notes=(
                "Mass information alone does not fix the projected outer-dipole "
                "carrier amplitude when geometry remains free."
            ),
        ),
        PriorScenario(
            name="mass_plus_geometry_weak",
            carrier_prior_precision=0.1,
            notes=(
                "Weak external geometry information gives a finite but small "
                "carrier-amplitude prior in normalized data units."
            ),
        ),
        PriorScenario(
            name="outer_amplitude_moderate",
            carrier_prior_precision=1.0,
            notes=(
                "A direct or effectively equivalent outer-carrier amplitude "
                "prior comparable to the timing information."
            ),
        ),
        PriorScenario(
            name="idealized_strong_prior",
            carrier_prior_precision=100.0,
            notes=(
                "Idealized strong external prior on the GR outer-dipole carrier "
                "amplitude."
            ),
        ),
    )


def normalized_fisher_matrix(prior_precision: float) -> tuple[tuple[float, float], ...]:
    """Return Fisher matrix for [A_GR_dipole, delta_alpha_clock].

    The normalized data columns are identical, so the data Fisher block is
    [[1, 1], [1, 1]]. The external prior constrains only A_GR_dipole.
    """

    if prior_precision < 0:
        raise ValueError("prior precision must be non-negative")
    return ((1.0 + prior_precision, 1.0), (1.0, 1.0))


def fisher_rank(prior_precision: float) -> int:
    return 1 if prior_precision == 0 else 2


def marginal_clock_information(prior_precision: float) -> float:
    """Return Schur-complement information for the clock amplitude."""

    if prior_precision == 0:
        return 0.0
    return prior_precision / (1.0 + prior_precision)


def condition_number(prior_precision: float) -> float:
    """Return the 2x2 Fisher condition number."""

    if prior_precision == 0:
        return inf
    trace = 2.0 + prior_precision
    determinant = prior_precision
    discriminant = sqrt(trace * trace - 4.0 * determinant)
    lambda_max = (trace + discriminant) / 2.0
    lambda_min = (trace - discriminant) / 2.0
    return lambda_max / lambda_min


def classify(prior_precision: float) -> tuple[str, str, str]:
    """Classify degeneracy, identifiability, and verdict."""

    info = marginal_clock_information(prior_precision)
    if prior_precision == 0:
        return ("exact", "none", "no-go")
    if info < 0.25:
        return ("near-exact", "weak", "conditional-no-go")
    if info < 0.75:
        return ("partial", "conditional", "conditional")
    return ("broken by prior", "practical", "motivates-runtime")


def format_float(value: float) -> str:
    if value == inf:
        return "inf"
    return f"{value:.6g}"


def audit_rows() -> list[PriorAuditRow]:
    rows: list[PriorAuditRow] = []
    for scenario in scenarios():
        degeneracy, identifiability, verdict = classify(
            scenario.carrier_prior_precision
        )
        rows.append(
            PriorAuditRow(
                prior_set=scenario.name,
                carrier_amplitude_degeneracy=degeneracy,
                clock_amplitude_identifiability=identifiability,
                carrier_prior_precision=format_float(
                    scenario.carrier_prior_precision
                ),
                fisher_rank=str(fisher_rank(scenario.carrier_prior_precision)),
                marginal_clock_information=format_float(
                    marginal_clock_information(scenario.carrier_prior_precision)
                ),
                condition_number=format_float(
                    condition_number(scenario.carrier_prior_precision)
                ),
                verdict=verdict,
                notes=scenario.notes,
            )
        )
    return rows


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in audit_rows()]


def validate_rows(rows: Iterable[dict[str, str]]) -> None:
    expected = set(PRIOR_COLUMNS)
    for index, row in enumerate(rows, start=1):
        actual = set(row)
        if actual != expected:
            raise ValueError(
                f"row {index} has columns {sorted(actual)}, "
                f"expected {sorted(expected)}"
            )
        if not all(str(row[column]).strip() for column in PRIOR_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(PRIOR_COLUMNS),
        "assumptions": {
            "model": "normalized two-amplitude Fisher audit",
            "parameters": "[A_GR_dipole, delta_alpha_clock]",
            "data_fisher": "[[1, 1], [1, 1]]",
            "prior": "diagonal precision on A_GR_dipole only",
            "free_fall_sector": "GR",
            "propagation_sector": "GR",
            "clock_sector": "decoupled EFT only",
            "clock_order": "leading O(c^-2)",
            "hierarchy_order": "leading O(epsilon=a_in/a_out)",
            "status": "conditional prior-breaking audit; not a runtime fit",
        },
        "interpretation": {
            "marginal_clock_information": (
                "Schur complement p/(1+p), where p is carrier prior precision"
            ),
            "core_degeneracy": (
                "Without a carrier-amplitude prior, delta_alpha_clock is exactly "
                "degenerate with A_GR_dipole."
            ),
            "runtime_gate": (
                "Full triple runtime work is only motivated if realistic external "
                "priors reach the conditional/practical rows."
            ),
        },
        "rows": rows,
    }
