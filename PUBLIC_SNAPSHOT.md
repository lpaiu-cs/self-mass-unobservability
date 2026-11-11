# Public snapshot of the self-mass-unobservability repository

This public repository receives snapshots of the author's working tree. The working history could not be hosted here because it carried about 50 GB of binary numerical arrays (GitHub rejects files larger than 100 MB).

- Each snapshot commit names, in its message, the commit of the working history it was taken from.
- Included: every tracked file except binary arrays (`*.npz`, `*.npy`, `*.gz`) and files larger than 50 MB. `PUBLIC_SNAPSHOT_OMITTED.tsv` in each snapshot lists the omitted paths with their sizes and git blob ids.
- On 11 November 2026, after the research line was closed, the author discarded the arrays. They are no longer available. `docs/discarded-arrays-2026-11-11.tsv` lists the 28,290 discarded files (50.17 GB) by path, size and git blob id, and the manifests (`paper/revision-manifest.json` and `outputs/**/*-manifest.json`) keep their SHA-256 digests. Thirty-two small arrays that were already part of the public history remain.
- The working history was then rewritten without the arrays. Commits that contained arrays received new identifiers, and `docs/commit-id-map-2026-11-11.tsv` maps each old identifier to its replacement. Commit identifiers cited in the paper, the notes and snapshot messages before 11 November 2026 are the old ones. The rewritten working history is not public; it is available from the author on request.
- Most notes under `notes/` are written in Korean.

The corrected submission snapshot is identified by the tag `prd-submission-2026-09-28` (a preparation version, not a claim that a journal submission occurred). Its main paper and Supplemental Material are in `paper/`; their PDFs are in `output/pdf/`. `paper/submission-manifest.json` binds the submission files, and `paper/revision-manifest.json` retains the research-input hashes.

With NumPy installed, `python verification/replay_public_inference.py` replays the six-coefficient omnibus and physical lag sections using the public `outputs/research-completion/public-inference.json`. This does not require the arrays. The broader `verify_unified_paper.py` and the run and audit scripts that read arrays can no longer be run; the compact replay does not certify timing derivatives or replace the full stellar/timing calculations.

The obsolete `clock-timing-dictionary` and `collapse-theorem-dynamic-visibility` branch histories are retained under `archive/2026-09-28/` tags. The retired `lpaiu/dynamic-chi-observable` tip is already an ancestor of public main. Local working checkouts are unaffected by public branch cleanup.
