from __future__ import annotations

from dataclasses import dataclass


DELTA_MAX = 4
LIVE_BOTTLENECK = (
    "the omitted rank-4 contraction EEQ and the resulting high-rank exhaustiveness patch before any move to Reven6+"
)


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
    live_bottleneck: str
    envelope_closed: bool
    smallest_unaudited_class: str | None
    entries: tuple[FamilyEnvelopeEntry, ...]


def family_envelope_entries() -> tuple[FamilyEnvelopeEntry, ...]:
    return (
        FamilyEnvelopeEntry(
            class_id="R2",
            envelope_state="audited",
            tensor_rank="2 STF",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="X2 / EX",
            current_theorem_role="audited uniqueness obstruction and threshold class",
            rank_order=2,
        ),
        FamilyEnvelopeEntry(
            class_id="R0a",
            envelope_state="audited",
            tensor_rank="0",
            parity="even",
            derivative_character="bare source allowed",
            smallest_expected_witness_type="S",
            current_theorem_role="audited uniqueness obstruction and threshold class",
            rank_order=0,
        ),
        FamilyEnvelopeEntry(
            class_id="R0b",
            envelope_state="audited",
            tensor_rank="0",
            parity="even",
            derivative_character="derivative-only or shift-symmetric",
            smallest_expected_witness_type="dotS2 / DtS_E2",
            current_theorem_role="audited uniqueness obstruction and threshold class",
            rank_order=0,
        ),
        FamilyEnvelopeEntry(
            class_id="R1",
            envelope_state="audited",
            tensor_rank="1",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="V2 / EVV",
            current_theorem_role="audited uniqueness obstruction and threshold class",
            rank_order=1,
        ),
        FamilyEnvelopeEntry(
            class_id="Rodd+",
            envelope_state="audited",
            tensor_rank="odd >= 3",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="T2 / ETT",
            current_theorem_role="audited uniqueness obstruction and threshold class; enlarged audited-set composition now re-closed",
            rank_order=3,
        ),
        FamilyEnvelopeEntry(
            class_id="Reven4+",
            envelope_state="audited",
            tensor_rank="even >= 4",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="Q2 / EEQ and EQQ",
            current_theorem_role="audited uniqueness obstruction and threshold class; enlarged audited-set composition now re-closed, but exhaustive bookkeeping patch still active",
            rank_order=4,
        ),
        FamilyEnvelopeEntry(
            class_id="Rodd5+",
            envelope_state="audited",
            tensor_rank="odd >= 5",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="U2 / EUU",
            current_theorem_role="audited uniqueness obstruction and threshold class; enlarged audited-set composition now re-closed",
            rank_order=5,
        ),
        FamilyEnvelopeEntry(
            class_id="Reven6+",
            envelope_state="still unaudited",
            tensor_rank="even >= 6",
            parity="even",
            derivative_character="unsuppressed local family",
            smallest_expected_witness_type="rank-6 self-contraction type",
            current_theorem_role="next smallest unaudited family gate after the closed Rodd5+ composition layer",
            rank_order=6,
        ),
        FamilyEnvelopeEntry(
            class_id="Podd",
            envelope_state="excluded by explicit assumption",
            tensor_rank="any",
            parity="odd",
            derivative_character="any",
            smallest_expected_witness_type="not in MVP domain",
            current_theorem_role="excluded by A2",
            rank_order=99,
        ),
        FamilyEnvelopeEntry(
            class_id="Spin",
            envelope_state="excluded by explicit assumption",
            tensor_rank="any",
            parity="any",
            derivative_character="spin- or orientation-carrying",
            smallest_expected_witness_type="not in MVP domain",
            current_theorem_role="excluded by A2",
            rank_order=99,
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
            rank_order=99,
        ),
    )


def family_envelope_summary(delta_max: int = DELTA_MAX) -> FamilyEnvelopeSummary:
    if delta_max != 4:
        raise ValueError("The current family-envelope census is only fixed at Delta <= 4.")

    entries = family_envelope_entries()
    unaudited = tuple(
        entry for entry in entries if entry.envelope_state == "still unaudited"
    )
    smallest = min(unaudited, key=lambda item: (item.rank_order, item.class_id), default=None)
    return FamilyEnvelopeSummary(
        delta_max=delta_max,
        live_bottleneck=LIVE_BOTTLENECK,
        envelope_closed=not unaudited,
        smallest_unaudited_class=None if smallest is None else smallest.class_id,
        entries=entries,
    )


def family_envelope_report(delta_max: int = DELTA_MAX) -> str:
    summary = family_envelope_summary(delta_max=delta_max)
    lines = [
        "key\tvalue",
        f"delta_max\t{summary.delta_max}",
        f"live_bottleneck\t{summary.live_bottleneck}",
        f"envelope_closed\t{str(summary.envelope_closed).lower()}",
        f"smallest_unaudited_class\t{summary.smallest_unaudited_class or 'none'}",
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
