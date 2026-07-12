from __future__ import annotations

from high_rank_family_enumerator import (
    available_block_names,
    generated_family_classes,
    high_rank_audit_summary,
    manual_class_key,
    manual_survivor_classes,
)


def high_rank_diff_report() -> str:
    lines = [
        "rank\tfamily\tstatus\tlabel\tweight\tsignature\tkey",
    ]
    for rank in (3, 4, 5):
        family = {3: "T", 4: "Q", 5: "U"}[rank]
        generated = {item.key: item for item in generated_family_classes(rank)}
        # manual classes outside the enumerator's block universe (DtE after
        # the rank-4 completion) cannot be keyed; they are reported below as
        # manual_only with an out-of-universe marker.
        universe = set(available_block_names(rank))
        manual_items = manual_survivor_classes(rank)
        out_of_universe = [item for item in manual_items if not set(item.signature) <= universe]
        manual = {
            manual_class_key(rank, item): item
            for item in manual_items
            if set(item.signature) <= universe
        }

        for key in sorted(generated):
            item = generated[key]
            status = "matched" if key in manual else "omitted_from_manual"
            lines.append(
                "\t".join(
                    (
                        str(rank),
                        family,
                        status,
                        item.label,
                        str(item.weight),
                        ",".join(item.signature),
                        repr(item.key),
                    )
                )
            )

        for key in sorted(manual):
            if key in generated:
                continue
            item = manual[key]
            lines.append(
                "\t".join(
                    (
                        str(rank),
                        family,
                        "manual_only",
                        item.label,
                        str(item.weight),
                        ",".join(item.signature),
                        repr(key),
                    )
                )
            )

        for item in sorted(out_of_universe, key=lambda entry: entry.label):
            lines.append(
                "\t".join(
                    (
                        str(rank),
                        family,
                        "manual_only",
                        item.label,
                        str(item.weight),
                        ",".join(item.signature),
                        "out-of-universe",
                    )
                )
            )

    lines.append("")
    lines.append("rank\tfamily\tgenerated_count\tmanual_count\tmatched_count\tomitted_from_manual\tmanual_only")
    for rank in (3, 4, 5):
        summary = high_rank_audit_summary(rank)
        lines.append(
            "\t".join(
                (
                    str(summary.rank),
                    summary.family_name,
                    str(summary.generated_count),
                    str(summary.manual_count),
                    str(summary.matched_count),
                    ",".join(summary.omitted_from_manual) if summary.omitted_from_manual else "-",
                    ",".join(summary.manual_only) if summary.manual_only else "-",
                )
            )
        )
    return "\n".join(lines)


def main() -> None:
    print(high_rank_diff_report())


if __name__ == "__main__":
    main()
