PYTHON ?= $(shell command -v python3.12 2>/dev/null || command -v python3)

.PHONY: request1 request2 request3 request5-phaseA request6 request7 paper-tex paper-pdf core

request1:
	cd request1 && $(PYTHON) request1_com_decoupling.py

request2:
	cd request2 && $(PYTHON) request2_internal_structure.py

request3:
	cd request3 && $(PYTHON) request3_llr_mock.py

request5-phaseA:
	cd request5 && $(PYTHON) request5_j0337_phaseA.py

request6:
	cd request6 && $(PYTHON) request6_clock_sector.py

request7:
	cd request7 && $(PYTHON) request7_joint_consistency_scaffold.py

paper-tex:
	$(MAKE) -C paper tex

paper-pdf:
	$(MAKE) -C paper pdf

core: request1 request2 request3 request5-phaseA request6 request7
