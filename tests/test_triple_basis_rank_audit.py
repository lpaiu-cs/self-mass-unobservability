from __future__ import annotations

import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from smu_clock.triple_basis_rank_audit import (  # noqa: E402
    leading_outer_dipole_harmonic_vector,
    rows_as_dicts,
    validate_rows,
)


class TripleBasisRankAuditTests(unittest.TestCase):
    def test_outer_dipole_vector_contains_combination_frequencies(self) -> None:
        vector = leading_outer_dipole_harmonic_vector()

        self.assertIn("cos(lambda_in - lambda_out)", vector)
        self.assertIn("cos(2 lambda_in - lambda_out)", vector)
        self.assertIn("cos(lambda_in - 2 lambda_out)", vector)

    def test_audit_rows_are_classified_and_valid(self) -> None:
        rows = rows_as_dicts()

        validate_rows(rows)
        self.assertTrue(rows)
        self.assertTrue(
            all(row["observable_status"] == "degenerate" for row in rows)
        )


if __name__ == "__main__":
    unittest.main()
