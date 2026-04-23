PYTHON ?= $(shell command -v python3.12 2>/dev/null || command -v python3)

.PHONY: request1 request2 request3 request5-phaseA request6 request7 core

request1:
	$(PYTHON) request1_com_decoupling.py

request2:
	$(PYTHON) request2_internal_structure.py

request3:
	$(PYTHON) request3_llr_mock.py

request5-phaseA:
	$(PYTHON) request5_j0337_phaseA.py

request6:
	$(PYTHON) request6_clock_sector.py

request7:
	$(PYTHON) request7_joint_consistency_scaffold.py

core: request1 request2 request3 request5-phaseA request6 request7
