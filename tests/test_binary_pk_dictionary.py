from __future__ import annotations

import sys
import unittest
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from smu_clock.binary_pk_dictionary import (  # noqa: E402
    ALLOWED_STATUSES,
    TABLE_COLUMNS,
    build_rate_projections,
    derive_binary_clock_expansion,
    derive_binary_clock_exact_e,
    dictionary_rows,
    exact_e_dictionary_rows,
    truncate_eccentricity,
    unexpected_rate_harmonic_remainders,
)


class BinaryPkDictionaryTests(unittest.TestCase):
    def test_no_independent_rate_harmonic_beyond_constant_and_cos_e(self) -> None:
        expansion = derive_binary_clock_expansion(eccentricity_order=1)

        self.assertEqual(unexpected_rate_harmonic_remainders(expansion), {})
        for term in expansion.rate_terms:
            self.assertLessEqual(
                set(term.split.nonzero_harmonics()),
                {"constant", "cos(E)"},
            )

    def test_every_dictionary_row_is_classified(self) -> None:
        rows = dictionary_rows(eccentricity_order=1)

        self.assertTrue(rows)
        for row in rows:
            self.assertEqual(set(row), set(TABLE_COLUMNS))
            self.assertIn(row["observable_status"], ALLOWED_STATUSES)
            for column in TABLE_COLUMNS:
                self.assertTrue(str(row[column]).strip())

    def test_zeta_harmonics_match_potential_carrier_structure(self) -> None:
        expansion = derive_binary_clock_expansion(eccentricity_order=1)
        terms = {term.name: term for term in expansion.rate_terms}
        s_a = sp.Symbol("s_A")
        zeta_1 = sp.Symbol("zeta_1")
        zeta_2 = sp.Symbol("zeta_2")

        potential = terms["baseline_potential"].split
        zeta1 = terms["zeta1_sA_clock"].split
        zeta2 = terms["zeta2_sA2_clock"].split

        self.assertEqual(
            sp.simplify(zeta1.constant - zeta_1 * s_a * potential.constant),
            0,
        )
        self.assertEqual(
            sp.simplify(zeta1.cos_e - zeta_1 * s_a * potential.cos_e),
            0,
        )
        self.assertEqual(
            sp.simplify(zeta2.constant - zeta_2 * s_a**2 * potential.constant),
            0,
        )
        self.assertEqual(
            sp.simplify(zeta2.cos_e - zeta_2 * s_a**2 * potential.cos_e),
            0,
        )

    def test_exact_e_rows_are_classified(self) -> None:
        rows = exact_e_dictionary_rows()

        self.assertTrue(rows)
        for row in rows:
            self.assertEqual(set(row), set(TABLE_COLUMNS))
            self.assertIn(row["observable_status"], ALLOWED_STATUSES)
            for column in TABLE_COLUMNS:
                self.assertTrue(str(row[column]).strip())

    def test_exact_e_periodic_terms_are_einstein_delay_basis_only(self) -> None:
        exact = derive_binary_clock_exact_e()
        E = sp.Symbol("E")

        for contribution in exact.contributions:
            residual = contribution.periodic_residual
            coefficient = contribution.periodic_basis_coefficient
            self.assertEqual(
                sp.simplify(residual - coefficient * sp.sin(E)),
                0,
            )
            self.assertFalse(sp.simplify(coefficient).has(sp.sin(E), sp.cos(E)))

    def test_exact_e_truncates_to_existing_oe1_projection(self) -> None:
        oe1 = derive_binary_clock_expansion(eccentricity_order=1)
        oe1_projections = {
            (projection.source_term, projection.rate_harmonic): projection
            for projection in build_rate_projections(oe1)
        }
        exact = derive_binary_clock_exact_e()
        e = sp.Symbol("e")
        t = sp.Symbol("t")

        for contribution in exact.contributions:
            exact_secular = truncate_eccentricity(
                contribution.secular_slope * t,
                e,
                1,
            )
            oe1_secular = oe1_projections[
                (contribution.source_term, "constant")
            ].residual_expression
            self.assertEqual(sp.simplify(exact_secular - oe1_secular), 0)

            exact_periodic = truncate_eccentricity(
                contribution.periodic_residual,
                e,
                1,
            )
            oe1_periodic = oe1_projections[
                (contribution.source_term, "cos(E)")
            ].residual_expression
            self.assertEqual(sp.simplify(exact_periodic - oe1_periodic), 0)


if __name__ == "__main__":
    unittest.main()
