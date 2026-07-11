from __future__ import annotations

import sys
from pathlib import Path

SYMBOLIC_ROOT = Path(__file__).resolve().parents[1]
if str(SYMBOLIC_ROOT) not in sys.path:
    sys.path.insert(0, str(SYMBOLIC_ROOT))

from family_envelope_census import family_envelope_summary
from hereditary_kernel_demo import hereditary_kernel_summary
from irrep_family_census import irrep_family_summary
from nonanalytic_jet_demo import nonanalytic_jet_summary
from stateful_counterexample_demo import stateful_counterexample_summary
from weight_spectrum_demo import weight_spectrum_summary


def main() -> None:
    family = family_envelope_summary()
    assert family.envelope_closed
    assert family.first_open_class is None

    irrep = irrep_family_summary()
    assert irrep.envelope_closes_on_audited_classes
    assert irrep.first_open_nonstf_family is None

    weight = weight_spectrum_summary()
    assert weight.local_weight_spectrum_finiteness_suffices
    assert weight.a8_stronger_than_necessary
    assert weight.sharp_a8_failure_mode

    nonanalytic = nonanalytic_jet_summary()
    assert nonanalytic.locality_kept_in_all_cases
    assert nonanalytic.finite_family_operator_closure_kept_in_all_cases
    assert nonanalytic.smallest_local_nonanalytic_counterexample
    assert nonanalytic.broken_layer

    hereditary = hereditary_kernel_summary()
    assert hereditary.analyticity_kept_in_all_cases
    assert hereditary.finite_primitive_family_envelope_kept_in_all_cases
    assert hereditary.smallest_genuinely_hereditary_counterexample
    assert hereditary.broken_layer

    stateful = stateful_counterexample_summary()
    assert stateful.locality_kept_in_all_cases
    assert stateful.analyticity_kept_in_all_cases
    assert stateful.finite_primitive_family_envelope_kept_in_all_cases
    assert stateful.smallest_local_finite_state_counterexample
    assert stateful.finite_state_augmented_collapse_survives
    assert stateful.broken_layer

    print("symbolic smoke checks passed")


if __name__ == "__main__":
    main()
