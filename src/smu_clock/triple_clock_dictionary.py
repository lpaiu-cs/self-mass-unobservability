"""Leading hierarchical-triple outer-potential clock dictionary.

This module records symbolic dictionary rows for a Newtonian nested-Kepler
carrier at leading clock order O(c^-2) and leading nontrivial order in
epsilon = a_in / a_out. It is not a numerical timing model.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

from smu_clock.binary_pk_dictionary import ALLOWED_STATUSES, TABLE_COLUMNS


@dataclass(frozen=True)
class TripleDictionaryTerm:
    """One classified leading-order triple clock dictionary row."""

    term: str
    expansion_order: str
    harmonic: str
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


EXPANSION_ORDER = (
    "O(c^-2), Newtonian nested-Kepler carrier, "
    "leading nontrivial O(epsilon=a_in/a_out)"
)

COMPONENTS = (
    {
        "prefix": "baseline",
        "factor": "1",
        "kind": "GR carrier",
        "inner_target": "inner Einstein delay gamma",
        "outer_target": "outer Einstein/redshift delay parameter",
        "dipole_status": "degenerate",
        "dipole_target": "GR triple clock carrier / fitted triple parameters",
    },
    {
        "prefix": "zeta1_sA",
        "factor": "zeta_1 s_A",
        "kind": "linear decoupled clock term",
        "inner_target": "inner Einstein delay gamma",
        "outer_target": "outer Einstein/redshift delay parameter",
        "dipole_status": "observable",
        "dipole_target": "none at dictionary level; check covariance later",
    },
    {
        "prefix": "zeta2_sA2",
        "factor": "zeta_2 s_A^2",
        "kind": "quadratic decoupled clock term",
        "inner_target": "inner Einstein delay gamma",
        "outer_target": "outer Einstein/redshift delay parameter",
        "dipole_status": "observable",
        "dipole_target": "none at dictionary level; check covariance later",
    },
)


def notation_payload() -> dict[str, str]:
    """Return the notation used by the leading triple dictionary."""

    return {
        "inner_bodies": "A is the clock pulsar, B is the inner companion",
        "outer_body": "C orbits the inner-binary barycenter",
        "M_in": "m_A + m_B",
        "beta_A": "m_B / M_in, so rho_A = beta_A r_in",
        "epsilon": "a_in / a_out",
        "r_in": "a_in (1 - e_in cos E_in)",
        "R_out": "a_out (1 - e_out cos E_out)",
        "cos_Phi": "dot(rhat_in, Rhat_out), carrying relative orientation",
        "U_inner": "G m_B / r_in",
        "U_outer_expanded": (
            "G m_C/R_out + G m_C beta_A r_in cos_Phi/R_out^2 + O(epsilon^2)"
        ),
        "alpha_A": "1 + zeta_1 s_A + zeta_2 s_A^2",
        "clock_rate": (
            "d tau_A/dt = 1 - v_A^2/(2 c^2) "
            "- alpha_A (U_inner + U_outer)/c^2"
        ),
    }


def projection_payload() -> dict[str, str]:
    """Return symbolic projection structures used by the table rows."""

    return {
        "inner_monopole_projector": (
            "dt/dE_in = (1 - e_in cos E_in)/n_in; "
            "U_inner dt/dE_in integrates to E_in secular plus sin E_in"
        ),
        "outer_monopole_projector": (
            "dt/dE_out = (1 - e_out cos E_out)/n_out; "
            "outer monopole integrates to E_out secular plus sin E_out"
        ),
        "outer_dipole_projector": (
            "Delta_dipole = integral dt K_dipole "
            "[(1 - e_in cos E_in)/(1 - e_out cos E_out)^2] cos_Phi"
        ),
        "outer_dipole_kernel": (
            "K_dipole = -G m_C beta_A epsilon factor/(a_out c^2)"
        ),
        "harmonic_labels": "secular, inner-only, outer-only, combination",
    }


def _inner_rows(component: dict[str, str]) -> list[TripleDictionaryTerm]:
    prefix = component["prefix"]
    factor = component["factor"]
    kind = component["kind"]
    return [
        TripleDictionaryTerm(
            term=f"inner_{prefix}_secular",
            expansion_order=EXPANSION_ORDER,
            harmonic="secular",
            observable_status="absorbed",
            absorption_target="spin phase / spin frequency",
            notes=(
                f"{kind}; U_inner contribution with factor {factor}. "
                "Coordinate-time secular trend is absorbed by spin parameters."
            ),
        ),
        TripleDictionaryTerm(
            term=f"inner_{prefix}_einstein_sinE",
            expansion_order=EXPANSION_ORDER,
            harmonic="inner-only",
            observable_status="degenerate",
            absorption_target=component["inner_target"],
            notes=(
                f"{kind}; after secular absorption, residual is proportional to "
                f"-G e_in m_B ({factor}) sin E_in/(a_in c^2 n_in)."
            ),
        ),
    ]


def _outer_monopole_rows(component: dict[str, str]) -> list[TripleDictionaryTerm]:
    prefix = component["prefix"]
    factor = component["factor"]
    kind = component["kind"]
    return [
        TripleDictionaryTerm(
            term=f"outer_monopole_{prefix}_secular",
            expansion_order=EXPANSION_ORDER,
            harmonic="secular",
            observable_status="absorbed",
            absorption_target="spin phase / spin frequency",
            notes=(
                f"{kind}; outer monopole U_C(R_out) with factor {factor}. "
                "Coordinate-time secular trend is absorbed by spin parameters."
            ),
        ),
        TripleDictionaryTerm(
            term=f"outer_monopole_{prefix}_sinE",
            expansion_order=EXPANSION_ORDER,
            harmonic="outer-only",
            observable_status="degenerate",
            absorption_target=component["outer_target"],
            notes=(
                f"{kind}; after secular absorption, residual is proportional to "
                f"-G e_out m_C ({factor}) sin E_out/(a_out c^2 n_out)."
            ),
        ),
    ]


def _outer_dipole_row(component: dict[str, str]) -> TripleDictionaryTerm:
    prefix = component["prefix"]
    factor = component["factor"]
    kind = component["kind"]
    kernel = (
        "-G m_C beta_A epsilon "
        f"({factor}) [(1 - e_in cos E_in)/(1 - e_out cos E_out)^2] "
        "cos_Phi/(a_out c^2)"
    )
    return TripleDictionaryTerm(
        term=f"outer_dipole_{prefix}_combination",
        expansion_order=EXPANSION_ORDER,
        harmonic="combination",
        observable_status=component["dipole_status"],
        absorption_target=component["dipole_target"],
        notes=(
            f"{kind}; leading O(epsilon) outer-gradient clock term. "
            f"Residual structure is integral dt of {kernel}. "
            "This is not a pure spin, inner-only gamma, or outer-only gamma basis."
        ),
    )


def build_triple_dictionary_terms() -> list[TripleDictionaryTerm]:
    """Build classified rows for the leading triple clock dictionary."""

    rows: list[TripleDictionaryTerm] = []
    for component in COMPONENTS:
        rows.extend(_inner_rows(component))
        rows.extend(_outer_monopole_rows(component))
        rows.append(_outer_dipole_row(component))
    return rows


def triple_dictionary_rows() -> list[dict[str, str]]:
    """Return machine-readable table rows in stable column order."""

    return [term.to_row() for term in build_triple_dictionary_terms()]


def validate_triple_rows(rows: Iterable[dict[str, str]]) -> None:
    """Validate schema and classification completeness."""

    expected = set(TABLE_COLUMNS)
    allowed_harmonics = {"secular", "inner-only", "outer-only", "combination"}
    for index, row in enumerate(rows, start=1):
        actual = set(row)
        if actual != expected:
            raise ValueError(
                f"row {index} has columns {sorted(actual)}, "
                f"expected {sorted(expected)}"
            )
        if row["observable_status"] not in ALLOWED_STATUSES:
            raise ValueError(f"row {index} has invalid observable_status")
        if row["harmonic"] not in allowed_harmonics:
            raise ValueError(f"row {index} has invalid harmonic label")
        if not all(str(row[column]).strip() for column in TABLE_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")


def payload() -> dict[str, object]:
    """Return the complete JSON payload for the leading triple dictionary."""

    rows = triple_dictionary_rows()
    validate_triple_rows(rows)
    return {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "free_fall_sector": "GR",
            "propagation_sector": "GR",
            "clock_sector": "decoupled EFT only",
            "carrier_orbit": "Newtonian nested-Kepler hierarchical triple",
            "post_keplerian_order": "leading clock O(c^-2)",
            "epsilon_order": "through leading nontrivial O(a_in/a_out)",
            "status": "leading triple outer-potential dictionary; not a final theorem",
        },
        "notation": notation_payload(),
        "projection": projection_payload(),
        "rows": rows,
    }
