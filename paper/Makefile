PYTHON ?= $(shell command -v python3.12 2>/dev/null || command -v python3)

.PHONY: tex pdf clean

tex:
	$(PYTHON) build_manuscript.py

pdf: tex
	mkdir -p build
	latexmk -pdf -interaction=nonstopmode -halt-on-error -output-directory=build main.tex

clean:
	latexmk -C -output-directory=build main.tex
