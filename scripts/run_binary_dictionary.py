"""Emit the first populated binary clock dictionary table."""

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
    build_rate_projections,
    derive_binary_clock_expansion,
    derive_binary_clock_exact_e,
    dictionary_rows,
    exact_e_dictionary_rows,
    unexpected_rate_harmonic_remainders,
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
        description="Generate the Request 8 binary clock dictionary table."
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
        default=REPO_ROOT / "outputs" / "tsv" / "binary_pk_dictionary.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=REPO_ROOT / "outputs" / "json" / "binary_pk_dictionary.json",
        help="Output JSON path.",
    )
    parser.add_argument(
        "--exact-tsv",
        type=Path,
        default=REPO_ROOT / "outputs" / "tsv" / "binary_pk_dictionary_exact_e.tsv",
        help="Exact-in-e output TSV path.",
    )
    parser.add_argument(
        "--exact-json",
        type=Path,
        default=REPO_ROOT / "outputs" / "json" / "binary_pk_dictionary_exact_e.json",
        help="Exact-in-e output JSON path.",
    )
    args = parser.parse_args()

    expansion = derive_binary_clock_expansion(args.eccentricity_order)
    projections = build_rate_projections(expansion, args.eccentricity_order)
    rows = dictionary_rows(args.eccentricity_order)
    validate_rows(rows)

    exact_projection = derive_binary_clock_exact_e()
    exact_rows = exact_e_dictionary_rows(exact_projection)
    validate_rows(exact_rows)

    payload = {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "free_fall_sector": "GR",
            "propagation_sector": "GR",
            "clock_sector": "decoupled EFT only",
            "post_keplerian_order": "leading clock O(c^-2)",
            "carrier_orbit": "Newtonian Kepler binary",
            "eccentricity_order": args.eccentricity_order,
            "status": "first populated O(e^1) binary table; not a final theorem",
        },
        "symbolic_expansion": expansion.to_jsonable(),
        "unexpected_rate_harmonic_remainders": unexpected_rate_harmonic_remainders(
            expansion,
            args.eccentricity_order,
        ),
        "projections": [projection.to_jsonable() for projection in projections],
        "rows": rows,
    }
    exact_payload = {
        "schema": list(TABLE_COLUMNS),
        "assumptions": {
            "free_fall_sector": "GR",
            "propagation_sector": "GR",
            "clock_sector": "decoupled EFT only",
            "post_keplerian_order": "leading clock O(c^-2)",
            "carrier_orbit": "Newtonian Kepler binary",
            "eccentricity_order": "exact in e",
            "time_anomaly_jacobian": "dt/dE = (1 - e cos(E)) / n",
            "einstein_delay_basis": "sin(E)",
            "status": "exact-in-e leading-order binary table; not a final theorem",
        },
        "exact_projection": exact_projection.to_jsonable(),
        "rows": exact_rows,
    }

    write_tsv(args.tsv, rows)
    write_json(args.json, payload)
    write_tsv(args.exact_tsv, exact_rows)
    write_json(args.exact_json, exact_payload)

    print(f"wrote {args.tsv}")
    print(f"wrote {args.json}")
    print(f"wrote {args.exact_tsv}")
    print(f"wrote {args.exact_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
