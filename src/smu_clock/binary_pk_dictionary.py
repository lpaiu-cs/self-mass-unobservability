"""Binary clock-sector dictionary skeleton.

This module expands only the decoupled clock-rate map on a Newtonian binary
carrier orbit. It deliberately avoids triple terms, estimator calls, and source
chasing.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import sympy as sp


TABLE_COLUMNS = (
    "term",
    "expansion_order",
    "harmonic",
    "observable_status",
    "absorption_target",
    "notes",
)

ALLOWED_STATUSES = {"observable", "absorbed", "degenerate", "null"}


@dataclass(frozen=True)
class DictionaryTerm:
    """One observable-classification row for the binary dictionary."""

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


@dataclass(frozen=True)
class BinaryClockExpansion:
    """Symbolic pieces of d tau_A / dt at the declared truncation."""

    eccentricity_order: int
    alpha_a: sp.Expr
    delta_alpha_a: sp.Expr
    inverse_radius: sp.Expr
    potential_u_ext: sp.Expr
    velocity_squared_a: sp.Expr
    kinetic_rate_term: sp.Expr
    gr_potential_rate_term: sp.Expr
    zeta1_rate_term: sp.Expr
    zeta2_rate_term: sp.Expr

    def to_jsonable(self) -> dict[str, str | int]:
        return {
            "eccentricity_order": self.eccentricity_order,
            "alpha_a": str(self.alpha_a),
            "delta_alpha_a": str(self.delta_alpha_a),
            "inverse_radius": str(self.inverse_radius),
            "potential_u_ext": str(self.potential_u_ext),
            "velocity_squared_a": str(self.velocity_squared_a),
            "kinetic_rate_term": str(self.kinetic_rate_term),
            "gr_potential_rate_term": str(self.gr_potential_rate_term),
            "zeta1_rate_term": str(self.zeta1_rate_term),
            "zeta2_rate_term": str(self.zeta2_rate_term),
        }


def symbols() -> dict[str, sp.Symbol]:
    """Return the symbols used by the binary clock skeleton."""

    names = "G c m_A m_B M a e E s_A zeta_1 zeta_2"
    return {name: sp.Symbol(name) for name in names.split()}


def truncate_eccentricity(expr: sp.Expr, e: sp.Symbol, order: int) -> sp.Expr:
    """Expand and truncate an expression through O(e**order)."""

    if order < 0:
        raise ValueError("eccentricity order must be non-negative")
    return sp.series(expr, e, 0, order + 1).removeO().expand()


def derive_binary_clock_expansion(eccentricity_order: int = 1) -> BinaryClockExpansion:
    """Derive the symbolic binary clock-rate pieces for body A.

    The carrier orbit is Newtonian. The clock expansion is kept at leading
    post-Keplerian clock order O(c^-2), with eccentricity truncated by
    ``eccentricity_order``.
    """

    sym = symbols()
    G = sym["G"]
    c = sym["c"]
    m_b = sym["m_B"]
    total_mass = sym["M"]
    a = sym["a"]
    e = sym["e"]
    eccentric_anomaly = sym["E"]
    s_a = sym["s_A"]
    zeta_1 = sym["zeta_1"]
    zeta_2 = sym["zeta_2"]

    cos_e = sp.cos(eccentric_anomaly)
    inverse_radius_exact = 1 / (a * (1 - e * cos_e))
    inverse_radius = truncate_eccentricity(
        inverse_radius_exact, e, eccentricity_order
    )

    alpha_a = 1 + zeta_1 * s_a + zeta_2 * s_a**2
    delta_alpha_a = alpha_a - 1
    potential_u_ext = G * m_b * inverse_radius

    barycentric_factor_a = (m_b / total_mass) ** 2
    relative_velocity_squared = G * total_mass * (2 * inverse_radius - 1 / a)
    velocity_squared_a = truncate_eccentricity(
        barycentric_factor_a * relative_velocity_squared, e, eccentricity_order
    )

    kinetic_rate_term = truncate_eccentricity(
        -velocity_squared_a / (2 * c**2), e, eccentricity_order
    )
    gr_potential_rate_term = truncate_eccentricity(
        -potential_u_ext / c**2, e, eccentricity_order
    )
    zeta1_rate_term = truncate_eccentricity(
        -(zeta_1 * s_a) * potential_u_ext / c**2, e, eccentricity_order
    )
    zeta2_rate_term = truncate_eccentricity(
        -(zeta_2 * s_a**2) * potential_u_ext / c**2, e, eccentricity_order
    )

    return BinaryClockExpansion(
        eccentricity_order=eccentricity_order,
        alpha_a=alpha_a,
        delta_alpha_a=delta_alpha_a,
        inverse_radius=inverse_radius,
        potential_u_ext=potential_u_ext,
        velocity_squared_a=velocity_squared_a,
        kinetic_rate_term=kinetic_rate_term,
        gr_potential_rate_term=gr_potential_rate_term,
        zeta1_rate_term=zeta1_rate_term,
        zeta2_rate_term=zeta2_rate_term,
    )


def build_dictionary_terms(eccentricity_order: int = 1) -> list[DictionaryTerm]:
    """Build the first binary dictionary table stub.

    Rows are classification scaffolds at O(c^-2), not final claims.
    """

    order = f"O(c^-2), Newtonian carrier orbit, through O(e^{eccentricity_order})"
    return [
        DictionaryTerm(
            term="gr_kinetic_clock_rate_constant",
            expansion_order=order,
            harmonic="constant in d_tau/dt; secular in pulse phase",
            observable_status="absorbed",
            absorption_target="spin phase / spin frequency",
            notes="GR baseline term from v_A^2; included only to anchor the projection.",
        ),
        DictionaryTerm(
            term="gr_kinetic_clock_rate_cosE",
            expansion_order=order,
            harmonic="cos E in d_tau/dt; Einstein-delay harmonic after integration",
            observable_status="degenerate",
            absorption_target="Einstein delay gamma",
            notes="GR baseline harmonic; not a decoupled clock-sector signal.",
        ),
        DictionaryTerm(
            term="gr_potential_clock_rate_constant",
            expansion_order=order,
            harmonic="constant in d_tau/dt; secular in pulse phase",
            observable_status="absorbed",
            absorption_target="spin phase / spin frequency",
            notes="GR external-potential baseline from U_ext.",
        ),
        DictionaryTerm(
            term="gr_potential_clock_rate_cosE",
            expansion_order=order,
            harmonic="cos E in d_tau/dt; Einstein-delay harmonic after integration",
            observable_status="degenerate",
            absorption_target="Einstein delay gamma",
            notes="Standard binary potential harmonic at this truncation.",
        ),
        DictionaryTerm(
            term="zeta1_sensitivity_potential_constant",
            expansion_order=order,
            harmonic="constant in d_tau/dt; secular in pulse phase",
            observable_status="absorbed",
            absorption_target="spin phase / spin frequency",
            notes="Linear sensitivity correction rescales the orbit-averaged potential term.",
        ),
        DictionaryTerm(
            term="zeta1_sensitivity_potential_cosE",
            expansion_order=order,
            harmonic="cos E in d_tau/dt; Einstein-delay harmonic after integration",
            observable_status="degenerate",
            absorption_target="Einstein delay gamma",
            notes="Linear sensitivity correction inherits the binary U_ext harmonic.",
        ),
        DictionaryTerm(
            term="zeta2_sensitivity_squared_potential_constant",
            expansion_order=order,
            harmonic="constant in d_tau/dt; secular in pulse phase",
            observable_status="absorbed",
            absorption_target="spin phase / spin frequency",
            notes="Quadratic sensitivity correction rescales the orbit-averaged potential term.",
        ),
        DictionaryTerm(
            term="zeta2_sensitivity_squared_potential_cosE",
            expansion_order=order,
            harmonic="cos E in d_tau/dt; Einstein-delay harmonic after integration",
            observable_status="degenerate",
            absorption_target="Einstein delay gamma",
            notes="Quadratic sensitivity correction inherits the binary U_ext harmonic.",
        ),
        DictionaryTerm(
            term="independent_clock_harmonic_beyond_cosE",
            expansion_order=order,
            harmonic="none at declared truncation",
            observable_status="null",
            absorption_target="none",
            notes="Placeholder row: no independent harmonic is generated by the multiplicative U_ext correction at this skeleton order; not a final theorem.",
        ),
    ]


def dictionary_rows(eccentricity_order: int = 1) -> list[dict[str, str]]:
    """Return the table rows as dictionaries in stable column order."""

    return [term.to_row() for term in build_dictionary_terms(eccentricity_order)]


def validate_rows(rows: Iterable[dict[str, str]]) -> None:
    """Validate schema and classification completeness."""

    expected = set(TABLE_COLUMNS)
    for index, row in enumerate(rows, start=1):
        actual = set(row)
        if actual != expected:
            raise ValueError(
                f"row {index} has columns {sorted(actual)}, "
                f"expected {sorted(expected)}"
            )
        if row["observable_status"] not in ALLOWED_STATUSES:
            raise ValueError(f"row {index} has invalid observable_status")
        if not all(str(row[column]).strip() for column in TABLE_COLUMNS):
            raise ValueError(f"row {index} contains an empty field")
