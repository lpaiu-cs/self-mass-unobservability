from __future__ import annotations

from dataclasses import dataclass

from high_rank_family_enumerator import high_rank_audit_summary
from r3_sector_delta4 import r3_summary
from r4_sector_delta4 import r4_summary
from r5_sector_delta4 import r5_summary
from r6_sector_delta4 import r6_summary


@dataclass(frozen=True)
class StfRankPatternEntry:
    rank: int
    first_self_witness: str
    first_mixed_witness_layer: str
    threshold_type: str
    current_audited_status: str


@dataclass(frozen=True)
class StfRankPatternSummary:
    attempted_tower_theorem_holds: bool
    failure_rank: int | None
    failure_reason: str | None
    entries: tuple[StfRankPatternEntry, ...]


def stf_rank_pattern_entries() -> tuple[StfRankPatternEntry, ...]:
    r3 = r3_summary()
    r4 = r4_summary()
    r5 = r5_summary()
    r6 = r6_summary()
    rank4_exhaustive = high_rank_audit_summary(4)

    # The rank-4 mixed layer opens with BOTH cubics. Before the sector
    # completion EEQ surfaced only through the exhaustive audit's omission
    # set; after the completion it is the sector's own first mixed witness.
    rank4_layer = "EQQ"
    if r4.first_mixed_witness == "EEQ" or "EEQ" in rank4_exhaustive.omitted_from_manual:
        rank4_layer = "EEQ and EQQ"

    return (
        StfRankPatternEntry(
            rank=3,
            first_self_witness=r3.first_self_witness or "T2",
            first_mixed_witness_layer=r3.first_mixed_witness or "ETT",
            threshold_type="self-only sharp",
            current_audited_status="supports attempted odd-rank STF pattern",
        ),
        StfRankPatternEntry(
            rank=4,
            first_self_witness=r4.first_self_witness or "Q2",
            first_mixed_witness_layer=rank4_layer,
            threshold_type="self-only sharp threshold, but mixed formula exception",
            current_audited_status="explicit exception to the attempted L>=3 STF tower theorem",
        ),
        StfRankPatternEntry(
            rank=5,
            first_self_witness=r5.first_self_witness or "U2",
            first_mixed_witness_layer=r5.first_mixed_witness or "EUU",
            threshold_type="self-only sharp",
            current_audited_status="supports attempted odd-rank STF pattern",
        ),
        StfRankPatternEntry(
            rank=6,
            first_self_witness=r6.first_self_witness or "Z2",
            first_mixed_witness_layer=(
                ", ".join(r6.first_mixed_layer_labels)
                if r6.first_mixed_layer_labels
                else (r6.first_mixed_witness or "EZZ")
            ),
            threshold_type="self-only sharp",
            current_audited_status="supports the isolated-rank4 even-rank exception picture within the audited STF sector",
        ),
    )


def stf_rank_pattern_summary() -> StfRankPatternSummary:
    entries = stf_rank_pattern_entries()
    failure_entry = next(
        (entry for entry in entries if "exception" in entry.current_audited_status),
        None,
    )
    return StfRankPatternSummary(
        attempted_tower_theorem_holds=failure_entry is None,
        failure_rank=None if failure_entry is None else failure_entry.rank,
        failure_reason=(
            None
            if failure_entry is None
            else "rank-4 STF admits the additional mixed cubic EEQ, so the first mixed layer is not universally EYY"
        ),
        entries=entries,
    )


def stf_rank_pattern_report() -> str:
    summary = stf_rank_pattern_summary()
    lines = [
        "key\tvalue",
        f"attempted_tower_theorem_holds\t{str(summary.attempted_tower_theorem_holds).lower()}",
        f"failure_rank\t{summary.failure_rank if summary.failure_rank is not None else 'none'}",
        f"failure_reason\t{summary.failure_reason or 'none'}",
        "",
        "rank\tfirst_self_witness\tfirst_mixed_witness_layer\tthreshold_type\tcurrent_audited_status",
    ]
    for entry in summary.entries:
        lines.append(
            "\t".join(
                (
                    str(entry.rank),
                    entry.first_self_witness,
                    entry.first_mixed_witness_layer,
                    entry.threshold_type,
                    entry.current_audited_status,
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(stf_rank_pattern_report())


if __name__ == "__main__":
    main()
