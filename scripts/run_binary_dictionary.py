"""Emit the first binary clock dictionary table stub."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from smu_clock.binary_pk_dictionary import (  # noqa: E402
    TABLE_COLUMNS,
    derive_binary_clock_expansion,
    dictionary_rows,
    validate_rows,
)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TABLE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate the Request 8 binary clock dictionary table stub."
    )
    parser.add_argument(
        "--eccentricity-order",
        type=int,
        default=1,
        help="Truncate the symbolic eccentricity expansion through O(e^N).",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=REPO_ROOT / "outputs" / "tsv" / "binary_pk_dictionary_stub.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=REPO_ROOT / "outputs" / "json" / "binary_pk_dictionary_stub.json",
        help="Output JSON path.",
    )
    args = parser.parse_args()

    expansion = derive_binary_clock_expansion(args.eccentricity_order)
    rows = dictionary_rows(args.eccentricity_order)
    validate_rows(rows)

    payload = {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "free_fall_sector": "GR",
            "propagation_sector": "GR",
            "clock_sector": "decoupled EFT only",
            "post_keplerian_order": "leading clock O(c^-2)",
            "carrier_orbit": "Newtonian Kepler binary",
            "eccentricity_order": args.eccentricity_order,
            "status": "skeleton; not a final theorem or tied-vs-decoupled verdict",
        },
        "symbolic_expansion": expansion.to_jsonable(),
        "rows": rows,
    }

    write_tsv(args.tsv, rows)
    write_json(args.json, payload)

    print(f"wrote {args.tsv}")
    print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
