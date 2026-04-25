# Paper Draft

This directory keeps the manuscript source and the current LaTeX draft.

- `manuscript.md`: prose-first draft source copied into the repository.
- `dynamic-loophole-outline.md`: separate paper-facing outline for the A4
  dynamic-loophole branch after the M16 analytic phase-boundary result.
- `dynamic-loophole-manuscript.md`: prose-first dynamic-loophole manuscript
  package built from the outline without overwriting the older draft.
- `build_manuscript.py`: converts the markdown draft into `main.tex`.
- `main.tex`: generated LaTeX draft committed to the repository for review.
- `Makefile`: local build helpers.

Typical workflow:

```bash
make tex
make pdf
```

or from the repository root:

```bash
make paper-tex
make paper-pdf
```

The LaTeX draft is intentionally conservative. It keeps the current prose,
equations, and section structure, but it does not yet convert the selected
references into a full BibTeX-driven citation pass.
