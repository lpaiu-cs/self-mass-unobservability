PYTHON ?= python

.PHONY: chi-response chi-two-frequency chi-basis-audit frequency-sweep forcing-dictionary projection-audit triple-shared-tau-bridge triple-gr-carrier-inventory symbolic-check

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

symbolic-check:
	$(PYTHON) symbolic/checks/test_symbolic.py
