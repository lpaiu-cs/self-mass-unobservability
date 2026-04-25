from __future__ import annotations

from dataclasses import dataclass


DELTA_MAX = 4
LIVE_STATEMENT = (
    'the former live bottleneck "prove or refute irreducible family-envelope closure beyond the audited scalar/vector/STF classes" is now resolved positively'
)


@dataclass(frozen=True)
class IrrepFamilyEntry:
    class_id: str
    group: str
    resolution_state: str
    parity: str
    irreducible_character: str
    resolution_mechanism: str
    current_theorem_role: str


@dataclass(frozen=True)
class IrrepFamilySummary:
    delta_max: int
    live_statement: str
    envelope_closes_on_audited_classes: bool
    first_open_nonstf_family: str | None
    entries: tuple[IrrepFamilyEntry, ...]


def irrep_family_entries() -> tuple[IrrepFamilyEntry, ...]:
    return (
        IrrepFamilyEntry(
            class_id="Scalar",
            group="scalar",
            resolution_state="already audited",
            parity="even",
            irreducible_character="rank-0 ordinary tensor",
            resolution_mechanism="audited directly as R0a and R0b",
            current_theorem_role="audited irreducible family class",
        ),
        IrrepFamilyEntry(
            class_id="Vector",
            group="vector",
            resolution_state="already audited",
            parity="even",
            irreducible_character="rank-1 ordinary tensor",
            resolution_mechanism="audited directly as R1",
            current_theorem_role="audited irreducible family class",
        ),
        IrrepFamilyEntry(
            class_id="STF2",
            group="stf_rank2",
            resolution_state="already audited",
            parity="even",
            irreducible_character="rank-2 STF",
            resolution_mechanism="audited directly as the rank-2 STF special case",
            current_theorem_role="audited irreducible family class",
        ),
        IrrepFamilyEntry(
            class_id="STFge3",
            group="stf_rankL_ge3",
            resolution_state="already audited",
            parity="even",
            irreducible_character="rank-L STF with L >= 3",
            resolution_mechanism="covered by the audited STF self-witness threshold theorem",
            current_theorem_role="audited irreducible family theorem class",
        ),
        IrrepFamilyEntry(
            class_id="TraceDesc",
            group="mixed_symmetry_or_nonstf",
            resolution_state="absorbed by trace reduction",
            parity="even",
            irreducible_character="ordinary tensor with explicit delta traces",
            resolution_mechanism="reduced to lower-rank scalar/vector/STF classes",
            current_theorem_role="not a genuinely new primitive family",
        ),
        IrrepFamilyEntry(
            class_id="PseudoOdd",
            group="mixed_symmetry_or_nonstf",
            resolution_state="excluded by parity/nonspin assumptions",
            parity="odd or pseudo",
            irreducible_character="odd epsilon-dual sector",
            resolution_mechanism="excluded by A2 after irrep decomposition",
            current_theorem_role="outside the theorem domain",
        ),
        IrrepFamilyEntry(
            class_id="MixedEvenDual",
            group="mixed_symmetry_or_nonstf",
            resolution_state="absorbed by trace reduction",
            parity="even",
            irreducible_character="even-dual mixed-symmetry sector",
            resolution_mechanism="dualized to an ordinary lower-rank tensor and then reduced to scalar/vector/STF classes",
            current_theorem_role="not a genuinely new primitive family",
        ),
    )


def irrep_family_summary(delta_max: int = DELTA_MAX) -> IrrepFamilySummary:
    if delta_max != 4:
        raise ValueError("The irreducible family census is only fixed at Delta <= 4.")
    entries = irrep_family_entries()
    open_entries = tuple(
        entry for entry in entries if entry.resolution_state == "still genuinely open"
    )
    return IrrepFamilySummary(
        delta_max=delta_max,
        live_statement=LIVE_STATEMENT,
        envelope_closes_on_audited_classes=not open_entries,
        first_open_nonstf_family=None if not open_entries else open_entries[0].class_id,
        entries=entries,
    )


def irrep_family_report(delta_max: int = DELTA_MAX) -> str:
    summary = irrep_family_summary(delta_max=delta_max)
    lines = [
        "key\tvalue",
        f"delta_max\t{summary.delta_max}",
        f"live_statement\t{summary.live_statement}",
        f"envelope_closes_on_audited_classes\t{str(summary.envelope_closes_on_audited_classes).lower()}",
        f"first_open_nonstf_family\t{summary.first_open_nonstf_family or 'none'}",
        "",
        "class_id\tgroup\tresolution_state\tparity\tirreducible_character\tresolution_mechanism\tcurrent_theorem_role",
    ]
    for entry in summary.entries:
        lines.append(
            "\t".join(
                (
                    entry.class_id,
                    entry.group,
                    entry.resolution_state,
                    entry.parity,
                    entry.irreducible_character,
                    entry.resolution_mechanism,
                    entry.current_theorem_role,
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(irrep_family_report())


if __name__ == "__main__":
    main()
