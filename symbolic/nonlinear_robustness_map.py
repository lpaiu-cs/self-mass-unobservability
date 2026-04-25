from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import product

from component_separability_audit import SeparabilityDesign, classify_separability_design


@dataclass(frozen=True)
class RobustnessRow:
    p: float
    eccentricity: float
    rho: Fraction
    damping: float
    channel: str
    projection_mode: str
    harmonic_cutoff: int
    usable_count: int
    usable_harmonics: tuple[int, ...]
    generated_required: int
    component_verdict: str
    bracketed: bool
    verdict: str


@dataclass(frozen=True)
class RobustnessSummary:
    total: int
    verdict_counts: dict[str, int]
    channel_counts: dict[str, dict[str, int]]
    positive_fraction: float


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _projection_configs() -> tuple[tuple[str, str, int, bool], ...]:
    return (
        ("acceleration", "acceleration", 1, False),
        ("range", "range-fixed-kappa", 2, False),
        ("range", "range-free-kappa", 2, True),
    )


def _is_bracketed(design: SeparabilityDesign) -> bool:
    lower = any(harmonic < design.rho for harmonic in design.usable_harmonics)
    upper = any(harmonic > design.rho for harmonic in design.usable_harmonics)
    return lower and upper


def classify_robustness_row(
    p: float,
    eccentricity: float,
    rho: Fraction,
    damping: float,
    channel: str,
    projection_mode: str,
    projection_nuisance: int,
    include_projection_nuisance: bool,
    harmonic_cutoff: int = 6,
    relative_cutoff: float = 1.0e-3,
    sideband_order: int = 1,
    samples: int = 1024,
) -> RobustnessRow:
    design = classify_separability_design(
        f"p={p} e={eccentricity} rho={_format_fraction(rho)} {projection_mode}",
        p=p,
        eccentricity=eccentricity,
        rho=rho,
        damping=damping,
        channel=channel,
        harmonic_cutoff=harmonic_cutoff,
        criterion="relative",
        relative_cutoff=relative_cutoff,
        sideband_order=sideband_order,
        projection_nuisance=projection_nuisance,
        include_range_projection_nuisance=include_projection_nuisance,
        samples=samples,
    )
    bracketed = _is_bracketed(design)
    if design.verdict != "component-separable-and-budget-ready":
        verdict = design.verdict
    elif not bracketed:
        verdict = "budget-ready-without-resonance-bracket"
    else:
        verdict = "robust-positive"
    return RobustnessRow(
        p=p,
        eccentricity=eccentricity,
        rho=rho,
        damping=damping,
        channel=channel,
        projection_mode=projection_mode,
        harmonic_cutoff=harmonic_cutoff,
        usable_count=len(design.usable_harmonics),
        usable_harmonics=design.usable_harmonics,
        generated_required=design.generated_budget_required,
        component_verdict=design.verdict,
        bracketed=bracketed,
        verdict=verdict,
    )


def robustness_grid(
    p_values: tuple[float, ...] = (1.5, 2.0, 3.0),
    eccentricities: tuple[float, ...] = (0.1, 0.3, 0.6),
    rho_values: tuple[Fraction, ...] = (Fraction(4, 3), Fraction(3, 2), Fraction(2, 1)),
    damping_values: tuple[float, ...] = (0.1, 0.2),
    harmonic_cutoff: int = 6,
    relative_cutoff: float = 1.0e-3,
    samples: int = 1024,
) -> tuple[RobustnessRow, ...]:
    rows: list[RobustnessRow] = []
    for p, eccentricity, rho, damping, config in product(
        p_values,
        eccentricities,
        rho_values,
        damping_values,
        _projection_configs(),
    ):
        channel, projection_mode, projection_nuisance, include_projection_nuisance = config
        rows.append(
            classify_robustness_row(
                p,
                eccentricity,
                rho,
                damping,
                channel,
                projection_mode,
                projection_nuisance,
                include_projection_nuisance,
                harmonic_cutoff=harmonic_cutoff,
                relative_cutoff=relative_cutoff,
                samples=samples,
            )
        )
    return tuple(rows)


def summarize_rows(rows: tuple[RobustnessRow, ...]) -> RobustnessSummary:
    verdict_counts = Counter(row.verdict for row in rows)
    channel_counter: dict[str, Counter[str]] = {}
    for row in rows:
        channel_counter.setdefault(row.projection_mode, Counter())[row.verdict] += 1
    positive = verdict_counts.get("robust-positive", 0)
    return RobustnessSummary(
        total=len(rows),
        verdict_counts=dict(verdict_counts),
        channel_counts={
            key: dict(value)
            for key, value in sorted(channel_counter.items())
        },
        positive_fraction=0.0 if not rows else positive / len(rows),
    )


def representative_rows(rows: tuple[RobustnessRow, ...]) -> tuple[RobustnessRow, ...]:
    representatives: list[RobustnessRow] = []
    seen: set[str] = set()
    for row in rows:
        if row.verdict not in seen:
            representatives.append(row)
            seen.add(row.verdict)
    return tuple(representatives)


def _format_row(row: RobustnessRow) -> str:
    harmonics = ",".join(str(item) for item in row.usable_harmonics) or "none"
    return (
        f"p={row.p:g} e={row.eccentricity:g} rho={_format_fraction(row.rho)} "
        f"delta={row.damping:g} {row.projection_mode} | "
        f"usable={row.usable_count}/{row.generated_required} [{harmonics}] | "
        f"bracket={row.bracketed} | component={row.component_verdict} | "
        f"verdict={row.verdict}"
    )


def audit_report() -> str:
    rows = robustness_grid()
    summary = summarize_rows(rows)
    lines = [
        "Nonlinear second-order robustness map",
        "",
        "Grid:",
        "- p in {1.5, 2, 3}",
        "- e in {0.1, 0.3, 0.6}",
        "- rho in {4/3, 3/2, 2}",
        "- delta in {0.1, 0.2}",
        "- projection in {acceleration, range-fixed-kappa, range-free-kappa}",
        "- H=6, eta=1e-3, M=1",
        "",
        "Verdict counts:",
    ]
    for verdict, count in sorted(summary.verdict_counts.items()):
        lines.append(f"- {verdict}: {count}/{summary.total}")
    lines.append(f"- robust-positive fraction: {summary.positive_fraction:.3f}")
    lines.append("")
    lines.append("Channel counts:")
    for channel, counts in summary.channel_counts.items():
        detail = ", ".join(f"{verdict}={count}" for verdict, count in sorted(counts.items()))
        lines.append(f"- {channel}: {detail}")
    lines.append("")
    lines.append("Representative rows:")
    for row in representative_rows(rows):
        lines.append(f"- {_format_row(row)}")
    return "\n".join(lines)


def main() -> None:
    print(audit_report())


if __name__ == "__main__":
    main()
