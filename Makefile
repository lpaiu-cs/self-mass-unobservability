PYTHON ?= python

.PHONY: worldline-expand sensitivity-expand enumerate-basis normal-form-reduce symbolic-check legacy-request1 legacy-request2 legacy-request7

worldline-expand:
	$(PYTHON) symbolic/worldline_expand.py

sensitivity-expand:
	$(PYTHON) symbolic/sensitivity_expand.py

enumerate-basis:
	$(PYTHON) symbolic/enumerate_basis.py

normal-form-reduce:
	$(PYTHON) symbolic/normal_form_reduce.py

symbolic-check:
	$(PYTHON) symbolic/checks/test_symbolic.py

legacy-request1:
	cd request1 && $(PYTHON) request1_com_decoupling.py

legacy-request2:
	cd request2 && $(PYTHON) request2_internal_structure.py

legacy-request7:
	cd request7 && $(PYTHON) request7_joint_consistency_scaffold.py
