PYTHON ?= python

.PHONY: chi-response chi-two-frequency chi-basis-audit frequency-sweep forcing-dictionary projection-audit triple-shared-tau-bridge triple-gr-carrier-inventory triple-projection-nuisance-gate triple-projection-manifold-gate named-timing-model-projection-audit nutimo-runtime-worthiness-pilot nutimo-external-handoff-packet nonlinear-sideband symbolic-check

chi-response:
	$(PYTHON) symbolic/chi_relaxation_response.py

chi-two-frequency:
	$(PYTHON) symbolic/chi_two_frequency_response.py

chi-basis-audit:
	$(PYTHON) symbolic/chi_basis_audit.py

frequency-sweep:
	$(PYTHON) symbolic/frequency_sweep_distinguishability.py

forcing-dictionary:
	$(PYTHON) symbolic/forcing_observable_dictionary.py

projection-audit:
	$(PYTHON) symbolic/projection_channel_audit.py

triple-shared-tau-bridge:
	$(PYTHON) symbolic/triple_shared_tau_bridge.py

triple-gr-carrier-inventory:
	$(PYTHON) symbolic/triple_gr_carrier_inventory.py

triple-projection-nuisance-gate:
	$(PYTHON) symbolic/triple_projection_nuisance_gate.py

triple-projection-manifold-gate:
	$(PYTHON) symbolic/triple_projection_manifold_gate.py

named-timing-model-projection-audit:
	$(PYTHON) symbolic/named_timing_model_projection_audit.py

nutimo-runtime-worthiness-pilot:
	$(PYTHON) symbolic/nutimo_runtime_worthiness_pilot.py

nutimo-external-handoff-packet:
	$(PYTHON) symbolic/nutimo_external_handoff_packet.py

nonlinear-sideband:
	$(PYTHON) symbolic/nonlinear_sideband_test.py

symbolic-check:
	$(PYTHON) symbolic/checks/test_symbolic.py
