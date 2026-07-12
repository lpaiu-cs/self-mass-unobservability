PYTHON ?= $(shell command -v python3.12 2>/dev/null || command -v python3 || command -v python)
PAPER_A_TEX := paper-A-collapse-theorem.tex
PAPER_B_TEX := paper-B-dynamic-sep-limit.tex

.PHONY: tex pdf paper-a-tex paper-a-pdf paper-b-tex paper-b-pdf clean

tex:
	$(PYTHON) build_manuscript.py

pdf: tex
	mkdir -p build
	latexmk -pdf -interaction=nonstopmode -halt-on-error -output-directory=build main.tex

paper-a-tex:
	$(PYTHON) build_paper_a.py

paper-a-pdf: paper-a-tex
	mkdir -p build
	latexmk -pdf -interaction=nonstopmode -halt-on-error -output-directory=build $(PAPER_A_TEX)

paper-b-tex:
	$(PYTHON) build_paper_b.py

paper-b-pdf: paper-b-tex
	mkdir -p build
	latexmk -pdf -interaction=nonstopmode -halt-on-error -output-directory=build $(PAPER_B_TEX)

clean:
	latexmk -C -output-directory=build main.tex
	latexmk -C -output-directory=build $(PAPER_A_TEX)
	latexmk -C -output-directory=build $(PAPER_B_TEX)
