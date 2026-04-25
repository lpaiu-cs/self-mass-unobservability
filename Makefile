PYTHON ?= python

.PHONY: worldline-expand sensitivity-expand enumerate-basis enumerate-contractions survivor-rank primitive-attack eb-sector eb-rank normal-form-reduce nonlinear-comparator shared-tau-ratio sample-budget orbital-harmonic-budget symbolic-check legacy-request1 legacy-request2 legacy-request7

worldline-expand:
	$(PYTHON) symbolic/worldline_expand.py

sensitivity-expand:
	$(PYTHON) symbolic/sensitivity_expand.py

enumerate-basis:
	$(PYTHON) symbolic/enumerate_basis.py

enumerate-contractions:
	$(PYTHON) symbolic/enumerate_contractions_delta4.py

survivor-rank:
	$(PYTHON) symbolic/survivor_rank_check.py

primitive-attack:
	$(PYTHON) symbolic/primitive_family_attack.py

eb-sector:
	$(PYTHON) symbolic/eb_sector_delta4.py

eb-rank:
	$(PYTHON) symbolic/eb_survivor_rank_check.py

normal-form-reduce:
	$(PYTHON) symbolic/normal_form_reduce.py

nonlinear-comparator:
	$(PYTHON) symbolic/nonlinear_comparator_audit.py

shared-tau-ratio:
	$(PYTHON) symbolic/shared_tau_ratio_audit.py

sample-budget:
	$(PYTHON) symbolic/sample_budget_audit.py

orbital-harmonic-budget:
	$(PYTHON) symbolic/orbital_harmonic_budget_audit.py

symbolic-check:
	$(PYTHON) symbolic/checks/test_symbolic.py

legacy-request1:
	cd request1 && $(PYTHON) request1_com_decoupling.py

legacy-request2:
	cd request2 && $(PYTHON) request2_internal_structure.py

legacy-request7:
	cd request7 && $(PYTHON) request7_joint_consistency_scaffold.py
