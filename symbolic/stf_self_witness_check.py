from __future__ import annotations

from dataclasses import dataclass

from high_rank_family_enumerator import family_symbol, generated_family_classes
from r3_sector_delta4 import r3_summary
from r4_sector_delta4 import r4_summary
from r5_sector_delta4 import r5_summary
from r6_sector_delta4 import r6_summary


@dataclass(frozen=True)
class StfSelfWitnessEntry:
    rank: int
    family_symbol: str
    first_self_witness: str
    linear_scalar_exists: bool
    mixed_quadratic_witness_exists: bool
    first_mixed_witness_layer: str
    threshold_at_delta4: str
    current_audited_status: str


@dataclass(frozen=True)
class StfSelfWitnessSummary:
    theorem_holds: bool
    failure_rank: int | None
    failure_reason: str | None
    entries: tuple[StfSelfWitnessEntry, ...]


def _summary_for_rank(rank: int):
    return {
        3: r3_summary,
        4: r4_summary,
        5: r5_summary,
        6: r6_summary,
    }[rank]()


def _mixed_layer_label(rank: int) -> str:
    if rank == 4:
        return "EEQ and EQQ"
    summary = _summary_for_rank(rank)
    if rank == 6 and getattr(summary, "first_mixed_layer_labels", ()):
        return ", ".join(summary.first_mixed_layer_labels)
    return summary.first_mixed_witness


def _has_linear_scalar(rank: int) -> bool:
    symbol = family_symbol(rank)
    return any(item.signature == (symbol,) for item in generated_family_classes(rank))


def _has_mixed_quadratic(rank: int) -> bool:
    symbol = family_symbol(rank)
    return any(item.signature == ("E", symbol) for item in generated_family_classes(rank))


def _audited_status(rank: int) -> str:
    if rank == 4:
        return "audited threshold instance with separate mixed-pattern exception"
    if rank == 6:
        return "audited threshold instance; mixed layer reverts to the EYY-first audited pattern"
    return "audited threshold instance"


def stf_self_witness_entries() -> tuple[StfSelfWitnessEntry, ...]:
    entries: list[StfSelfWitnessEntry] = []
    for rank in (3, 4, 5, 6):
        symbol = family_symbol(rank)
        summary = _summary_for_rank(rank)
        entries.append(
            StfSelfWitnessEntry(
                rank=rank,
                family_symbol=symbol,
                first_self_witness=summary.first_self_witness,
                linear_scalar_exists=_has_linear_scalar(rank),
                mixed_quadratic_witness_exists=_has_mixed_quadratic(rank),
                first_mixed_witness_layer=_mixed_layer_label(rank),
                threshold_at_delta4="w_Y >= 3",
                current_audited_status=_audited_status(rank),
            )
        )
    return tuple(entries)


def stf_self_witness_summary() -> StfSelfWitnessSummary:
    entries = stf_self_witness_entries()
    failure_entry = next(
        (
            entry
            for entry in entries
            if entry.linear_scalar_exists or entry.mixed_quadratic_witness_exists or entry.first_self_witness == ""
        ),
        None,
    )
    return StfSelfWitnessSummary(
        theorem_holds=failure_entry is None,
        failure_rank=None if failure_entry is None else failure_entry.rank,
        failure_reason=(
            None
            if failure_entry is None
            else "a linear or mixed-quadratic STF witness appeared below the quadratic norm layer"
        ),
        entries=entries,
    )


def stf_self_witness_report() -> str:
    summary = stf_self_witness_summary()
    lines = [
        "key\tvalue",
        f"theorem_holds\t{str(summary.theorem_holds).lower()}",
        f"failure_rank\t{summary.failure_rank if summary.failure_rank is not None else 'none'}",
        f"failure_reason\t{summary.failure_reason or 'none'}",
        "",
        "rank\tfamily_symbol\tfirst_self_witness\tlinear_scalar_exists\tmixed_quadratic_witness_exists\tfirst_mixed_witness_layer\tthreshold_at_delta4\tcurrent_audited_status",
    ]
    for entry in summary.entries:
        lines.append(
            "\t".join(
                (
                    str(entry.rank),
                    entry.family_symbol,
                    entry.first_self_witness,
                    str(entry.linear_scalar_exists).lower(),
                    str(entry.mixed_quadratic_witness_exists).lower(),
                    entry.first_mixed_witness_layer,
                    entry.threshold_at_delta4,
                    entry.current_audited_status,
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(stf_self_witness_report())


if __name__ == "__main__":
    main()
