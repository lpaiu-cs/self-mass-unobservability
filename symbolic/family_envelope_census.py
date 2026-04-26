from __future__ import annotations

from dataclasses import dataclass

from irrep_family_census import PACKAGE_STATUS, irrep_family_summary


DELTA_MAX = 4


@dataclass(frozen=True)
class FamilyEnvelopeEntry:
    class_id: str
    envelope_state: str
    tensor_rank: str
    parity: str
    derivative_character: str
    smallest_expected_witness_type: str
    current_theorem_role: str
    rank_order: int


@dataclass(frozen=True)
class FamilyEnvelopeSummary:
    delta_max: int
    package_status: str
    envelope_closed: bool
    first_open_class: str | None
    entries: tuple[FamilyEnvelopeEntry, ...]


def family_envelope_entries() -> tuple[FamilyEnvelopeEntry, ...]:
    return (
        FamilyEnvelopeEntry(
            class_id="Scalar",
            envelope_state="audited",
            tensor_rank="0",
            parity="even",
            derivative_character="bare or derivative-only local family",
            smallest_expected_witness_type="S / dotS2",
            current_theorem_role="audited irreducible family class",
            rank_order=0,
        ),
        FamilyEnvelopeEntry(
            class_id="Vector",
            envelope_state="audited",
            tensor_rank="1",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="V2 / EVV",
            current_theorem_role="audited irreducible family class",
            rank_order=1,
        ),
        FamilyEnvelopeEntry(
            class_id="STF2",
            envelope_state="audited",
            tensor_rank="2 STF",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="X2 / EX",
            current_theorem_role="audited irreducible family class",
            rank_order=2,
        ),
        FamilyEnvelopeEntry(
            class_id="STFge3",
            envelope_state="audited",
            tensor_rank="L >= 3 STF",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="Y2 / class-limited mixed pattern",
            current_theorem_role="audited irreducible family theorem class",
            rank_order=3,
        ),
        FamilyEnvelopeEntry(
            class_id="TraceDesc",
            envelope_state="absorbed by trace reduction",
            tensor_rank="ordinary tensor with traces",
            parity="even",
            derivative_character="local descendant sector",
            smallest_expected_witness_type="reduces to lower-rank audited class",
            current_theorem_role="not a genuinely new primitive family",
            rank_order=90,
        ),
        FamilyEnvelopeEntry(
            class_id="PseudoOdd",
            envelope_state="excluded by explicit assumption",
            tensor_rank="epsilon-dual pseudo sector",
            parity="odd or pseudo",
            derivative_character="local irrep sector",
            smallest_expected_witness_type="outside current domain",
            current_theorem_role="excluded by A2",
            rank_order=91,
        ),
        FamilyEnvelopeEntry(
            class_id="MixedEvenDual",
            envelope_state="absorbed by trace reduction",
            tensor_rank="even-dual mixed-symmetry sector",
            parity="even",
            derivative_character="local irrep sector",
            smallest_expected_witness_type="reduces to STF + traces",
            current_theorem_role="not a genuinely new primitive family",
            rank_order=92,
        ),
        FamilyEnvelopeEntry(
            class_id="State",
            envelope_state="excluded by explicit assumption",
            tensor_rank="not a local tensor family",
            parity="any",
            derivative_character="orbital-timescale state variable",
            smallest_expected_witness_type="chi_A-type loophole",
            current_theorem_role="excluded by A4; tracked as loophole class",
            rank_order=99,
        ),
        FamilyEnvelopeEntry(
            class_id="Nonlocal",
            envelope_state="excluded by explicit assumption",
            tensor_rank="not a local primitive family",
            parity="any",
            derivative_character="hereditary or nonlocal kernel",
            smallest_expected_witness_type="retarded-kernel loophole",
            current_theorem_role="excluded by A3; tracked as loophole class",
            rank_order=100,
        ),
    )


def family_envelope_summary(delta_max: int = DELTA_MAX) -> FamilyEnvelopeSummary:
    if delta_max != 4:
        raise ValueError("The current family-envelope census is only fixed at Delta <= 4.")
    irrep_summary = irrep_family_summary(delta_max=delta_max)
    return FamilyEnvelopeSummary(
        delta_max=delta_max,
        package_status=PACKAGE_STATUS,
        envelope_closed=irrep_summary.envelope_closes_on_audited_classes,
        first_open_class=irrep_summary.first_open_nonstf_family,
        entries=family_envelope_entries(),
    )


def family_envelope_report(delta_max: int = DELTA_MAX) -> str:
    summary = family_envelope_summary(delta_max=delta_max)
    lines = [
        "key\tvalue",
        f"delta_max\t{summary.delta_max}",
        f"package_status\t{summary.package_status}",
        f"envelope_closed\t{str(summary.envelope_closed).lower()}",
        f"first_open_class\t{summary.first_open_class or 'none'}",
        "",
        (
            "class_id\tenvelope_state\ttensor_rank\tparity\tderivative_character\t"
            "smallest_expected_witness_type\tcurrent_theorem_role"
        ),
    ]
    for entry in summary.entries:
        lines.append(
            "\t".join(
                (
                    entry.class_id,
                    entry.envelope_state,
                    entry.tensor_rank,
                    entry.parity,
                    entry.derivative_character,
                    entry.smallest_expected_witness_type,
                    entry.current_theorem_role,
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(family_envelope_report())


if __name__ == "__main__":
    main()
