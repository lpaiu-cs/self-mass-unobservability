# GRG submission edition - 28 September 2026

Target: General Relativity and Gravitation, Research article.

The PRD-labelled package was a preparation snapshot, not a journal submission. This GRG edition is the submitted journal format when a separate receipt confirms submission. Preparation alone does not establish receipt.

- `manuscript-source.zip`: editable main manuscript, bibliography and embedded figure.
- `manuscript.pdf`: locally compiled 14-page main manuscript.
- `ESM_1.pdf`: 31-page Online Resource 1 for publication with the paper.
- `cover-letter.pdf`: 1-page editorial letter.
- `manifest.json`: SHA-256 bindings for this edition.

The abstract is 230 whitespace-separated words. The edition adds five keywords, the author-confirmed Seoul/Republic of Korea address, no-external-funding declaration and author contribution statement. It uses Springer supplementary-material terminology and a GRG cover letter. Main Sections 2-6 and the complete scientific body of the supplement are byte-identical to the corrected preparation sources. The original package and its tag remain preserved.

Validation: main, supplement and letter compiled with installed Tectonic; no undefined citations/references or overfull boxes; all 46 pages rendered and visually inspected. The four-file source ZIP recompiles independently to the same page text and rendered pixels on every page (PDF byte hashes differ due to build metadata). The built-in standalone cover compiler failed because its platform directories were unavailable; installed Tectonic succeeded. No numerical experiments were rerun for this formatting step.

**Imported from prior work.** Scientific scope is unchanged: finite-frequency identifiability results and a conditional J0337 application, with no empirical SEP bound or detected relaxation. Previous derivative uncertainty and the uncomputed thermal tail remain explicitly stated.

Rebuild the TeX from the markdown with `paper/build_manuscript.py` using absolute source/output arguments and `PANDOC` pointing to the installed Pandoc. Compile from this directory with Tectonic; the ZIP contains only the main manuscript dependencies. The supplementary source, bibliography and four figures are also retained here for editorial revisions.
