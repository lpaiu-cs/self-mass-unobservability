"""Symbolic basis-rank audit for the leading triple clock kernel.

This is a dictionary-level rank check, not a covariance fit and not a timing
runtime. It asks whether the decoupled outer-dipole clock kernels add a timing
basis outside a standard nested-Kepler triple nuisance span.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable


AUDIT_COLUMNS = (
    "term",
    "harmonic_basis",
    "minimal_span_status",
    "expanded_span_status",
    "observable_status",
    "absorption_target",
    "notes",
)

ALLOWED_STATUSES = {"observable", "absorbed", "degenerate", "null"}


@dataclass(frozen=True)
class BasisAuditRow:
    """One row in the triple basis-rank audit table."""

    term: str
    harmonic_basis: str
    minimal_span_status: str
    expanded_span_status: str
    observable_status: str
    absorption_target: str
    notes: str

    def __post_init__(self) -> None:
        if self.observable_status not in ALLOWED_STATUSES:
            allowed = ", ".join(sorted(ALLOWED_STATUSES))
            raise ValueError(
                f"observable_status must be one of {allowed}; "
                f"got {self.observable_status!r}"
            )

    def to_row(self) -> dict[str, str]:
        return asdict(self)


def leading_outer_dipole_harmonic_vector() -> dict[str, str]:
    """Return the leading small-e harmonic vector for the dipole kernel.

    This expands

        [(1 - e_in cos lambda_in)/(1 - e_out cos lambda_out)^2]
        cos(lambda_in - lambda_out)

    through first order in inner and outer eccentricity, using a coplanar
    phase model for the basis-rank audit.
    """

    return {
        "cos(lambda_in - lambda_out)": "1",
        "cos(2 lambda_in - lambda_out)": "-e_in/2",
        "cos(lambda_in - 2 lambda_out)": "e_out",
        "cos(lambda_in)": "e_out",
        "cos(lambda_out)": "-e_in/2",
    }


def nuisance_span_payload() -> dict[str, object]:
    """Return the nuisance spans used by the audit."""

    return {
        "minimal_span": [
            "secular spin phase/frequency",
            "inner-only Einstein-delay basis: sin/cos(lambda_in)",
            "outer-only redshift-delay basis: sin/cos(lambda_out)",
        ],
        "expanded_standard_triple_span": [
            "minimal_span",
            "GR outer-gradient clock carrier amplitude",
            "GR outer-gradient clock carrier phase quadratures",
            "inner/outer orbital phase and amplitude redefinitions",
        ],
        "rank_conclusion": (
            "The leading decoupled outer-dipole kernels add a combination vector "
            "outside the minimal span, but do not add rank once the GR "
            "outer-gradient carrier/geometry nuisance span is included."
        ),
    }


def audit_rows() -> list[BasisAuditRow]:
    """Return classified rank-audit rows."""

    vector = ", ".join(
        f"{basis}:{coefficient}"
        for basis, coefficient in leading_outer_dipole_harmonic_vector().items()
    )
    combination_note = (
        "Contains combination frequencies lambda_in-lambda_out, "
        "2lambda_in-lambda_out, and lambda_in-2lambda_out, plus inner/outer "
        "sidebands at O(e)."
    )
    return [
        BasisAuditRow(
            term="gr_outer_dipole_carrier",
            harmonic_basis=vector,
            minimal_span_status="outside minimal spin/inner-only/outer-only span",
            expanded_span_status="basis generator inside standard GR triple carrier span",
            observable_status="degenerate",
            absorption_target="GR triple carrier / geometry nuisance",
            notes=(
                f"{combination_note} This row defines the standard carrier "
                "basis used by the expanded nuisance span."
            ),
        ),
        BasisAuditRow(
            term="zeta1_outer_dipole_kernel",
            harmonic_basis=f"(zeta_1 s_A) * [{vector}]",
            minimal_span_status="adds one combination vector outside minimal span",
            expanded_span_status="rank unchanged after adding GR carrier vector",
            observable_status="degenerate",
            absorption_target="GR outer-gradient clock amplitude / geometry nuisance",
            notes=(
                "The linear clock-sector kernel is a scalar multiple of the GR "
                "outer-dipole carrier. It is not independent in timing-basis "
                "rank unless external priors fix the carrier amplitude tightly."
            ),
        ),
        BasisAuditRow(
            term="zeta2_outer_dipole_kernel",
            harmonic_basis=f"(zeta_2 s_A^2) * [{vector}]",
            minimal_span_status="adds one combination vector outside minimal span",
            expanded_span_status="rank unchanged after adding GR carrier vector",
            observable_status="degenerate",
            absorption_target="GR outer-gradient clock amplitude / geometry nuisance",
            notes=(
                "The quadratic clock-sector kernel is a scalar multiple of the "
                "same GR outer-dipole carrier. It is also mutually degenerate "
                "with the zeta1 kernel in a single-clock basis audit."
            ),
        ),
        BasisAuditRow(
            term="effective_clock_outer_dipole_amplitude",
            harmonic_basis=(
                "(zeta_1 s_A + zeta_2 s_A^2) * [GR outer-dipole vector]"
            ),
            minimal_span_status="one effective combination amplitude outside minimal span",
            expanded_span_status="one effective amplitude degenerate with carrier amplitude",
            observable_status="degenerate",
            absorption_target="single effective GR-carrier amplitude rescaling",
            notes=(
                "At this level the two EFT coefficients enter as one effective "
                "single-clock amplitude. This is a leading triple no-go "
                "candidate, not a final theorem."
            ),
        ),
    ]


def rows_as_dicts() -> list[dict[str, str]]:
    return [row.to_row() for row in audit_rows()]


def validate_rows(rows: Iterable[dict[str, str]]) -> None:
    expected = set(AUDIT_COLUMNS)
    for index, row in enumerate(rows, start=1):
        actual = set(row)
        if actual != expected:
            raise ValueError(
                f"row {index} has columns {sorted(actual)}, "
                f"expected {sorted(expected)}"
            )
        if row["observable_status"] not in ALLOWED_STATUSES:
            raise ValueError(f"row {index} has invalid observable_status")
        if not all(str(row[column]).strip() for column in AUDIT_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


def payload() -> dict[str, object]:
    rows = rows_as_dicts()
    validate_rows(rows)
    return {
        "schema": list(AUDIT_COLUMNS),
        "assumptions": {
            "free_fall_sector": "GR",
            "propagation_sector": "GR",
            "clock_sector": "decoupled EFT only",
            "carrier_orbit": "Newtonian nested-Kepler hierarchical triple",
            "clock_order": "leading O(c^-2)",
            "hierarchy_order": "leading O(epsilon=a_in/a_out)",
            "eccentricity_order_for_harmonic_audit": "first order in e_in and e_out",
            "status": "basis-rank audit; not a covariance fit or final theorem",
        },
        "nuisance_spans": nuisance_span_payload(),
        "leading_outer_dipole_vector": leading_outer_dipole_harmonic_vector(),
        "rows": rows,
    }
