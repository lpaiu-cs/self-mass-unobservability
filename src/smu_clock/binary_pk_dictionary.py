"""Binary clock-sector dictionary at leading post-Keplerian order.

The O(e^1) engine expands the decoupled clock-rate map on a Newtonian binary
carrier orbit, splits each rate term into constant and cos(E) harmonics, and
maps those rate harmonics into timing-residual structures. The exact-in-e
engine rewrites the same clock EFT in d tau/dE form before integrating. Both
paths deliberately avoid triple systems, external estimators, and source
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
RATE_HARMONICS = ("constant", "cos(E)")


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
class HarmonicSplit:
    """Constant/cos(E) split of a rate expression at O(e^1)."""

    constant: sp.Expr
    cos_e: sp.Expr
    remainder: sp.Expr

    def coefficient(self, harmonic: str) -> sp.Expr:
        if harmonic == "constant":
            return self.constant
        if harmonic == "cos(E)":
            return self.cos_e
        raise ValueError(f"unsupported harmonic {harmonic!r}")

    def nonzero_harmonics(self) -> tuple[str, ...]:
        harmonics: list[str] = []
        if not is_zero(self.constant):
            harmonics.append("constant")
        if not is_zero(self.cos_e):
            harmonics.append("cos(E)")
        if not is_zero(self.remainder):
            harmonics.append("other")
        return tuple(harmonics)

    def to_jsonable(self) -> dict[str, str]:
        return {
            "constant": sstr(self.constant),
            "cos(E)": sstr(self.cos_e),
            "remainder": sstr(self.remainder),
        }


@dataclass(frozen=True)
class RateTerm:
    """One symbolic piece of d tau_A / dt."""

    name: str
    label: str
    expression: sp.Expr
    split: HarmonicSplit
    notes: str

    def to_jsonable(self) -> dict[str, object]:
        return {
            "name": self.name,
            "label": self.label,
            "expression": sstr(self.expression),
            "harmonic_split": self.split.to_jsonable(),
            "nonzero_harmonics": list(self.split.nonzero_harmonics()),
            "notes": self.notes,
        }


@dataclass(frozen=True)
class RateProjection:
    """Projection from a rate harmonic to a timing-residual structure."""

    source_term: str
    source_label: str
    rate_harmonic: str
    rate_coefficient: sp.Expr
    residual_structure: str
    residual_expression: sp.Expr
    observable_status: str
    absorption_target: str
    notes: str

    @property
    def row_term(self) -> str:
        suffix = "cosE" if self.rate_harmonic == "cos(E)" else self.rate_harmonic
        return f"{self.source_term}_{suffix}"

    def to_dictionary_term(self, expansion_order: str) -> DictionaryTerm:
        harmonic = (
            f"rate: {self.rate_harmonic}; residual: {self.residual_structure}"
        )
        notes = (
            f"{self.notes} rate_coefficient={sstr(self.rate_coefficient)}; "
            f"residual={sstr(self.residual_expression)}"
        )
        return DictionaryTerm(
            term=self.row_term,
            expansion_order=expansion_order,
            harmonic=harmonic,
            observable_status=self.observable_status,
            absorption_target=self.absorption_target,
            notes=notes,
        )

    def to_jsonable(self) -> dict[str, str]:
        return {
            "term": self.row_term,
            "source_term": self.source_term,
            "source_label": self.source_label,
            "rate_harmonic": self.rate_harmonic,
            "rate_coefficient": sstr(self.rate_coefficient),
            "residual_structure": self.residual_structure,
            "residual_expression": sstr(self.residual_expression),
            "observable_status": self.observable_status,
            "absorption_target": self.absorption_target,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class ExactContribution:
    """One exact-in-e contribution to the binary timing residual."""

    source_term: str
    source_label: str
    rate_dt_expression: sp.Expr
    rate_dE_expression: sp.Expr
    integrated_residual: sp.Expr
    secular_slope: sp.Expr
    absorbed_coordinate_time_residual: sp.Expr
    periodic_residual: sp.Expr
    periodic_basis: sp.Expr
    periodic_basis_coefficient: sp.Expr
    notes: str

    def to_dictionary_terms(self, expansion_order: str) -> tuple[DictionaryTerm, ...]:
        secular = DictionaryTerm(
            term=f"{self.source_term}_secular",
            expansion_order=expansion_order,
            harmonic="integrated residual: E secular; absorbed coordinate-time trend",
            observable_status="absorbed",
            absorption_target="spin phase / spin frequency",
            notes=(
                f"{self.notes} d_tau_dE={sstr(self.rate_dE_expression)}; "
                f"absorbed_slope={sstr(self.secular_slope)}"
            ),
        )

        if is_zero(self.periodic_residual):
            periodic_status = "null"
            periodic_target = "none"
            periodic_note = "periodic residual cancels at exact-in-e leading order"
        else:
            periodic_status = "degenerate"
            periodic_target = "Einstein delay gamma"
            periodic_note = (
                "periodic residual is proportional to the Einstein-delay sin(E) basis"
            )

        periodic = DictionaryTerm(
            term=f"{self.source_term}_einstein_sinE",
            expansion_order=expansion_order,
            harmonic="periodic residual: sin(E) Einstein-delay basis",
            observable_status=periodic_status,
            absorption_target=periodic_target,
            notes=(
                f"{periodic_note}; residual={sstr(self.periodic_residual)}; "
                f"basis_coefficient={sstr(self.periodic_basis_coefficient)}"
            ),
        )
        return (secular, periodic)

    def to_jsonable(self) -> dict[str, str]:
        return {
            "source_term": self.source_term,
            "source_label": self.source_label,
            "d_tau_dt": sstr(self.rate_dt_expression),
            "d_tau_dE": sstr(self.rate_dE_expression),
            "integrated_residual": sstr(self.integrated_residual),
            "secular_slope_coordinate_time": sstr(self.secular_slope),
            "absorbed_coordinate_time_residual": sstr(
                self.absorbed_coordinate_time_residual
            ),
            "periodic_residual_after_absorption": sstr(self.periodic_residual),
            "periodic_basis": sstr(self.periodic_basis),
            "periodic_basis_coefficient": sstr(self.periodic_basis_coefficient),
            "notes": self.notes,
        }


@dataclass(frozen=True)
class ExactBinaryClockProjection:
    """Exact-in-e leading-order binary clock projection."""

    radius: sp.Expr
    inverse_radius: sp.Expr
    dt_dE: sp.Expr
    coordinate_time: sp.Expr
    alpha_a: sp.Expr
    delta_alpha_a: sp.Expr
    potential_u_ext: sp.Expr
    velocity_squared_a: sp.Expr
    einstein_delay_basis: sp.Expr
    contributions: tuple[ExactContribution, ...]

    def to_jsonable(self) -> dict[str, object]:
        return {
            "radius": sstr(self.radius),
            "inverse_radius": sstr(self.inverse_radius),
            "dt_dE": sstr(self.dt_dE),
            "coordinate_time": sstr(self.coordinate_time),
            "alpha_a": sstr(self.alpha_a),
            "delta_alpha_a": sstr(self.delta_alpha_a),
            "potential_u_ext": sstr(self.potential_u_ext),
            "velocity_squared_a": sstr(self.velocity_squared_a),
            "einstein_delay_basis": sstr(self.einstein_delay_basis),
            "contributions": [
                contribution.to_jsonable()
                for contribution in self.contributions
            ],
        }


@dataclass(frozen=True)
class BinaryClockExpansion:
    """Symbolic binary clock-rate pieces at the declared truncation."""

    eccentricity_order: int
    radius: sp.Expr
    inverse_radius: sp.Expr
    alpha_a: sp.Expr
    delta_alpha_a: sp.Expr
    potential_u_ext: sp.Expr
    velocity_squared_a: sp.Expr
    rate_terms: tuple[RateTerm, ...]

    def to_jsonable(self) -> dict[str, object]:
        return {
            "eccentricity_order": self.eccentricity_order,
            "radius": sstr(self.radius),
            "inverse_radius": sstr(self.inverse_radius),
            "alpha_a": sstr(self.alpha_a),
            "delta_alpha_a": sstr(self.delta_alpha_a),
            "potential_u_ext": sstr(self.potential_u_ext),
            "velocity_squared_a": sstr(self.velocity_squared_a),
            "rate_terms": [term.to_jsonable() for term in self.rate_terms],
        }


def sstr(expr: sp.Expr) -> str:
    """Return a stable one-line SymPy string."""

    return sp.sstr(sp.factor(expr))


def is_zero(expr: sp.Expr) -> bool:
    """Return True when SymPy can reduce an expression to zero."""

    return sp.simplify(expr) == 0


def symbols() -> dict[str, sp.Symbol]:
    """Return the symbols used by the binary clock engine."""

    names = "G c m_A m_B M a e E n t s_A zeta_1 zeta_2"
    return {name: sp.Symbol(name) for name in names.split()}


def truncate_eccentricity(expr: sp.Expr, e: sp.Symbol, order: int) -> sp.Expr:
    """Expand and truncate an expression through O(e**order)."""

    if order < 0:
        raise ValueError("eccentricity order must be non-negative")
    return sp.series(expr, e, 0, order + 1).removeO().expand()


def split_rate_harmonics(expr: sp.Expr, eccentric_anomaly: sp.Symbol) -> HarmonicSplit:
    """Split a rate expression into constant, cos(E), and residual pieces."""

    cos_e = sp.cos(eccentric_anomaly)
    expanded = sp.expand(expr)
    constant = sp.expand(expanded.subs(cos_e, 0))
    cos_e_coefficient = sp.expand(expanded.coeff(cos_e, 1))
    remainder = sp.simplify(expanded - constant - cos_e_coefficient * cos_e)
    return HarmonicSplit(
        constant=constant,
        cos_e=cos_e_coefficient,
        remainder=remainder,
    )


def make_rate_term(
    name: str,
    label: str,
    expression: sp.Expr,
    eccentric_anomaly: sp.Symbol,
    notes: str,
) -> RateTerm:
    """Create a rate term with its harmonic split."""

    return RateTerm(
        name=name,
        label=label,
        expression=sp.expand(expression),
        split=split_rate_harmonics(expression, eccentric_anomaly),
        notes=notes,
    )


def make_exact_contribution(
    source_term: str,
    source_label: str,
    rate_dt_expression: sp.Expr,
    dt_dE: sp.Expr,
    coordinate_time: sp.Expr,
    eccentric_anomaly: sp.Symbol,
    einstein_delay_basis: sp.Expr,
    notes: str,
) -> ExactContribution:
    """Project one exact d tau/dt contribution through eccentric anomaly."""

    rate_dE_expression = sp.cancel(rate_dt_expression * dt_dE)
    rate_dE_expression = sp.expand(rate_dE_expression)
    integrated_residual = sp.expand(sp.integrate(rate_dE_expression, eccentric_anomaly))
    secular_e_coefficient = sp.expand(integrated_residual).coeff(eccentric_anomaly, 1)
    n = symbols()["n"]
    secular_slope = sp.simplify(secular_e_coefficient * n)
    absorbed_coordinate_time_residual = sp.expand(secular_slope * coordinate_time)
    periodic_residual = sp.simplify(
        integrated_residual - absorbed_coordinate_time_residual
    )
    periodic_basis_coefficient = (
        sp.simplify(periodic_residual / einstein_delay_basis)
        if not is_zero(periodic_residual)
        else sp.Integer(0)
    )

    return ExactContribution(
        source_term=source_term,
        source_label=source_label,
        rate_dt_expression=rate_dt_expression,
        rate_dE_expression=rate_dE_expression,
        integrated_residual=integrated_residual,
        secular_slope=secular_slope,
        absorbed_coordinate_time_residual=absorbed_coordinate_time_residual,
        periodic_residual=periodic_residual,
        periodic_basis=einstein_delay_basis,
        periodic_basis_coefficient=periodic_basis_coefficient,
        notes=notes,
    )


def derive_binary_clock_expansion(eccentricity_order: int = 1) -> BinaryClockExpansion:
    """Derive symbolic binary clock-rate pieces for body A at O(e^1).

    The carrier orbit is Newtonian. The clock expansion is kept at leading
    post-Keplerian clock order O(c^-2). This first real engine intentionally
    supports only the declared O(e^1) binary MVP.
    """

    if eccentricity_order != 1:
        raise ValueError("the current binary dictionary engine is fixed to O(e^1)")

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
    radius = a * (1 - e * cos_e)
    inverse_radius = truncate_eccentricity(1 / radius, e, eccentricity_order)
    potential_u_ext = sp.expand(G * m_b * inverse_radius)

    barycentric_factor_a = (m_b / total_mass) ** 2
    relative_velocity_squared = G * total_mass * (2 * inverse_radius - 1 / a)
    velocity_squared_a = truncate_eccentricity(
        barycentric_factor_a * relative_velocity_squared,
        e,
        eccentricity_order,
    )

    alpha_a = 1 + zeta_1 * s_a + zeta_2 * s_a**2
    delta_alpha_a = alpha_a - 1

    kinetic_rate_term = truncate_eccentricity(
        -velocity_squared_a / (2 * c**2),
        e,
        eccentricity_order,
    )
    baseline_potential_rate_term = truncate_eccentricity(
        -potential_u_ext / c**2,
        e,
        eccentricity_order,
    )
    zeta1_rate_term = truncate_eccentricity(
        -(zeta_1 * s_a) * potential_u_ext / c**2,
        e,
        eccentricity_order,
    )
    zeta2_rate_term = truncate_eccentricity(
        -(zeta_2 * s_a**2) * potential_u_ext / c**2,
        e,
        eccentricity_order,
    )

    rate_terms = (
        make_rate_term(
            name="baseline_kinetic",
            label="-v_A^2/(2 c^2)",
            expression=kinetic_rate_term,
            eccentric_anomaly=eccentric_anomaly,
            notes="Baseline special-relativistic clock term on the GR carrier.",
        ),
        make_rate_term(
            name="baseline_potential",
            label="-U_ext/c^2",
            expression=baseline_potential_rate_term,
            eccentric_anomaly=eccentric_anomaly,
            notes="Baseline GR external-potential clock term.",
        ),
        make_rate_term(
            name="zeta1_sA_clock",
            label="-zeta_1 s_A U_ext/c^2",
            expression=zeta1_rate_term,
            eccentric_anomaly=eccentric_anomaly,
            notes="Linear sensitivity correction in the decoupled clock sector.",
        ),
        make_rate_term(
            name="zeta2_sA2_clock",
            label="-zeta_2 s_A^2 U_ext/c^2",
            expression=zeta2_rate_term,
            eccentric_anomaly=eccentric_anomaly,
            notes="Quadratic sensitivity correction in the decoupled clock sector.",
        ),
    )

    return BinaryClockExpansion(
        eccentricity_order=eccentricity_order,
        radius=radius,
        inverse_radius=inverse_radius,
        alpha_a=alpha_a,
        delta_alpha_a=delta_alpha_a,
        potential_u_ext=potential_u_ext,
        velocity_squared_a=velocity_squared_a,
        rate_terms=rate_terms,
    )


def derive_binary_clock_exact_e() -> ExactBinaryClockProjection:
    """Derive the exact-in-e leading-order binary timing projection.

    The carrier orbit remains Newtonian and the only non-GR ingredient remains
    the decoupled clock EFT. The exact step is the anomaly-time Jacobian
    ``dt/dE = (1 - e cos(E)) / n``.
    """

    sym = symbols()
    G = sym["G"]
    c = sym["c"]
    m_b = sym["m_B"]
    total_mass = sym["M"]
    a = sym["a"]
    e = sym["e"]
    eccentric_anomaly = sym["E"]
    n = sym["n"]
    s_a = sym["s_A"]
    zeta_1 = sym["zeta_1"]
    zeta_2 = sym["zeta_2"]

    cos_e = sp.cos(eccentric_anomaly)
    sin_e = sp.sin(eccentric_anomaly)
    radius = a * (1 - e * cos_e)
    inverse_radius = 1 / radius
    dt_dE = (1 - e * cos_e) / n
    coordinate_time = (eccentric_anomaly - e * sin_e) / n
    einstein_delay_basis = sin_e

    potential_u_ext = G * m_b * inverse_radius
    barycentric_factor_a = (m_b / total_mass) ** 2
    relative_velocity_squared = G * total_mass * (2 * inverse_radius - 1 / a)
    velocity_squared_a = sp.cancel(barycentric_factor_a * relative_velocity_squared)

    alpha_a = 1 + zeta_1 * s_a + zeta_2 * s_a**2
    delta_alpha_a = alpha_a - 1

    rate_specs = (
        (
            "baseline_kinetic",
            "-v_A^2/(2 c^2)",
            sp.cancel(-velocity_squared_a / (2 * c**2)),
            "Baseline special-relativistic clock term on the GR carrier.",
        ),
        (
            "baseline_potential",
            "-U_ext/c^2",
            sp.cancel(-potential_u_ext / c**2),
            "Baseline GR external-potential clock term.",
        ),
        (
            "zeta1_sA_clock",
            "-zeta_1 s_A U_ext/c^2",
            sp.cancel(-(zeta_1 * s_a) * potential_u_ext / c**2),
            "Linear sensitivity correction in the decoupled clock sector.",
        ),
        (
            "zeta2_sA2_clock",
            "-zeta_2 s_A^2 U_ext/c^2",
            sp.cancel(-(zeta_2 * s_a**2) * potential_u_ext / c**2),
            "Quadratic sensitivity correction in the decoupled clock sector.",
        ),
    )

    contributions = tuple(
        make_exact_contribution(
            source_term=name,
            source_label=label,
            rate_dt_expression=rate_dt_expression,
            dt_dE=dt_dE,
            coordinate_time=coordinate_time,
            eccentric_anomaly=eccentric_anomaly,
            einstein_delay_basis=einstein_delay_basis,
            notes=notes,
        )
        for name, label, rate_dt_expression, notes in rate_specs
    )

    return ExactBinaryClockProjection(
        radius=radius,
        inverse_radius=inverse_radius,
        dt_dE=dt_dE,
        coordinate_time=coordinate_time,
        alpha_a=alpha_a,
        delta_alpha_a=delta_alpha_a,
        potential_u_ext=potential_u_ext,
        velocity_squared_a=velocity_squared_a,
        einstein_delay_basis=einstein_delay_basis,
        contributions=contributions,
    )


def project_rate_term(rate_term: RateTerm) -> tuple[RateProjection, ...]:
    """Map one rate term into timing-residual structures.

    At this truncation, a constant rate coefficient maps to secular pulse-phase
    structure, while a cos(E) rate coefficient maps to the standard sin(E)
    Einstein-delay residual harmonic.
    """

    sym = symbols()
    eccentric_anomaly = sym["E"]
    n = sym["n"]
    t = sym["t"]

    projections: list[RateProjection] = []
    constant_coefficient = rate_term.split.coefficient("constant")
    if not is_zero(constant_coefficient):
        projections.append(
            RateProjection(
                source_term=rate_term.name,
                source_label=rate_term.label,
                rate_harmonic="constant",
                rate_coefficient=constant_coefficient,
                residual_structure="secular spin phase/frequency",
                residual_expression=sp.expand(constant_coefficient * t),
                observable_status="absorbed",
                absorption_target="spin phase / spin frequency",
                notes=(
                    f"{rate_term.notes} The constant rate piece is a secular "
                    "clock rescaling at this order."
                ),
            )
        )

    cos_e_coefficient = rate_term.split.coefficient("cos(E)")
    if not is_zero(cos_e_coefficient):
        projections.append(
            RateProjection(
                source_term=rate_term.name,
                source_label=rate_term.label,
                rate_harmonic="cos(E)",
                rate_coefficient=cos_e_coefficient,
                residual_structure="sin(E) Einstein-delay harmonic",
                residual_expression=sp.expand(
                    cos_e_coefficient * sp.sin(eccentric_anomaly) / n
                ),
                observable_status="degenerate",
                absorption_target="Einstein delay gamma",
                notes=(
                    f"{rate_term.notes} The periodic piece has the standard "
                    "Einstein-delay timing shape at O(e^1)."
                ),
            )
        )

    return tuple(projections)


def build_rate_projections(
    expansion: BinaryClockExpansion | None = None,
    eccentricity_order: int = 1,
) -> tuple[RateProjection, ...]:
    """Build all timing projections for the binary dictionary."""

    if expansion is None:
        expansion = derive_binary_clock_expansion(eccentricity_order)

    projections: list[RateProjection] = []
    for rate_term in expansion.rate_terms:
        projections.extend(project_rate_term(rate_term))
    return tuple(projections)


def unexpected_rate_harmonic_remainders(
    expansion: BinaryClockExpansion | None = None,
    eccentricity_order: int = 1,
) -> dict[str, str]:
    """Return nonzero harmonic remainders beyond constant and cos(E)."""

    if expansion is None:
        expansion = derive_binary_clock_expansion(eccentricity_order)

    return {
        term.name: sstr(term.split.remainder)
        for term in expansion.rate_terms
        if not is_zero(term.split.remainder)
    }


def expansion_order_label(eccentricity_order: int = 1) -> str:
    return (
        "O(c^-2), Newtonian carrier orbit, "
        f"eccentricity through O(e^{eccentricity_order})"
    )


def exact_e_expansion_order_label() -> str:
    return "O(c^-2), Newtonian carrier orbit, exact in eccentricity e"


def build_dictionary_terms(eccentricity_order: int = 1) -> list[DictionaryTerm]:
    """Build the populated binary dictionary table from symbolic projections."""

    expansion = derive_binary_clock_expansion(eccentricity_order)
    projections = build_rate_projections(expansion, eccentricity_order)
    order = expansion_order_label(eccentricity_order)
    return [projection.to_dictionary_term(order) for projection in projections]


def build_exact_e_dictionary_terms(
    projection: ExactBinaryClockProjection | None = None,
) -> list[DictionaryTerm]:
    """Build the exact-in-e dictionary table."""

    if projection is None:
        projection = derive_binary_clock_exact_e()

    order = exact_e_expansion_order_label()
    rows: list[DictionaryTerm] = []
    for contribution in projection.contributions:
        rows.extend(contribution.to_dictionary_terms(order))
    return rows


def dictionary_rows(eccentricity_order: int = 1) -> list[dict[str, str]]:
    """Return the table rows as dictionaries in stable column order."""

    return [term.to_row() for term in build_dictionary_terms(eccentricity_order)]


def exact_e_dictionary_rows(
    projection: ExactBinaryClockProjection | None = None,
) -> list[dict[str, str]]:
    """Return exact-in-e table rows as dictionaries in stable column order."""

    return [
        term.to_row()
        for term in build_exact_e_dictionary_terms(projection)
    ]


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
