# Paper Draft

This directory keeps the manuscript source and the current LaTeX draft.

- `manuscript.md`: prose-first draft source copied into the repository.
- `build_manuscript.py`: converts the markdown draft into `main.tex` (requires `pandoc`; `pip install pypandoc-binary` provides a bundled binary).
- `main.tex`: generated LaTeX draft committed to the repository for review.
- `paper-A-collapse-theorem.md`: Paper A theorem-track manuscript source.
- `build_paper_a.py`: converts Paper A into `paper-A-collapse-theorem.tex` without requiring `pandoc`.
- `references.bib`: citation metadata for the selected references.
- `Makefile`: local build helpers.

Typical workflow:

```bash
make tex
make pdf
make paper-a-tex
make paper-a-pdf
```

or from the repository root:

```bash
make paper-tex
make paper-pdf
```

The LaTeX draft is intentionally conservative. It keeps the current prose,
equations, and section structure. Citation metadata for the selected
references is provided in `references.bib`; wiring the prose into a full
BibTeX-driven `\cite` pass is still left for the journal-style revision.
