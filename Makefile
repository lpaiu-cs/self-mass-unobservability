PYTHON ?= python

.PHONY: chi-response chi-two-frequency symbolic-check

chi-response:
	$(PYTHON) symbolic/chi_relaxation_response.py

chi-two-frequency:
	$(PYTHON) symbolic/chi_two_frequency_response.py

symbolic-check:
	$(PYTHON) symbolic/checks/test_symbolic.py
