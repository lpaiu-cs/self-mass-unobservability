PYTHON ?= python

.PHONY: tex pdf check clean

tex:
	$(PYTHON) build_figures.py
	$(PYTHON) build_manuscript.py

pdf: tex
	$(PYTHON) -c "from pathlib import Path; Path('build').mkdir(exist_ok=True)"
	latexmk -pdf -interaction=nonstopmode -halt-on-error -output-directory=build main.tex

check:
	$(PYTHON) ../verification/verify_unified_paper.py

clean:
	latexmk -C -output-directory=build main.tex
