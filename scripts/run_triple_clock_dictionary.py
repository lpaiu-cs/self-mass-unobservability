"""Emit the leading hierarchical-triple clock dictionary."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from smu_clock.binary_pk_dictionary import TABLE_COLUMNS  # noqa: E402
from smu_clock.triple_clock_dictionary import payload  # noqa: E402


TSV_PATH = REPO_ROOT / "outputs" / "tsv" / "triple_clock_dictionary_leading.tsv"
JSON_PATH = REPO_ROOT / "outputs" / "json" / "triple_clock_dictionary_leading.json"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TABLE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    data = payload()
    rows = data["rows"]
    if not isinstance(rows, list):
        raise TypeError("payload rows must be a list")

    write_tsv(TSV_PATH, rows)
    write_json(JSON_PATH, data)

    print(f"wrote {TSV_PATH}")
    print(f"wrote {JSON_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
