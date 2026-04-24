PYTHON ?= python

.PHONY: chi-response chi-two-frequency chi-basis-audit frequency-sweep forcing-dictionary symbolic-check

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

symbolic-check:
	$(PYTHON) symbolic/checks/test_symbolic.py
