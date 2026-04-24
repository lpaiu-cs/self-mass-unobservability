from __future__ import annotations

import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from smu_clock.triple_prior_identifiability import (  # noqa: E402
    fisher_rank,
    marginal_clock_information,
    rows_as_dicts,
    validate_rows,
)


class TriplePriorIdentifiabilityTests(unittest.TestCase):
    def test_no_prior_is_rank_deficient(self) -> None:
        self.assertEqual(fisher_rank(0.0), 1)
        self.assertEqual(marginal_clock_information(0.0), 0.0)

    def test_carrier_prior_breaks_rank(self) -> None:
        self.assertEqual(fisher_rank(0.1), 2)
        self.assertGreater(marginal_clock_information(0.1), 0.0)
        self.assertGreater(marginal_clock_information(100.0), 0.9)

    def test_rows_are_valid(self) -> None:
        rows = rows_as_dicts()
        validate_rows(rows)
        self.assertEqual(rows[0]["verdict"], "no-go")
        self.assertEqual(rows[-1]["verdict"], "motivates-runtime")


if __name__ == "__main__":
    unittest.main()
