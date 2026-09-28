# Public snapshot of the self-mass-unobservability repository

The author's working repository holds about 51.6 GB in about 71,600 files, most of it binary numerical arrays. GitHub cannot host that history (it rejects files larger than 100 MB), so this public repository receives snapshots of the working tree.

- Each snapshot commit names, in its message, the commit of the full history it was taken from.
- Included: every tracked file except binary arrays (`*.npz`, `*.npy`, `*.gz`) and files larger than 50 MB.
- Omitted: those arrays, about 50 GB. The manifests that bind them (`paper/revision-manifest.json` and `outputs/**/*-manifest.json`) list their SHA-256 digests. `PUBLIC_SNAPSHOT_OMITTED.tsv` in each snapshot lists every omitted path with its size and git blob id.
- Commit identifiers cited in the paper and in the notes refer to the full history, which is not public because of its size. The full history and the omitted arrays are available from the author on request.
- Most notes under `notes/` are written in Korean.

The corrected submission snapshot is identified by the tag `prd-submission-2026-09-28` (a preparation version, not a claim that a journal submission occurred). Its main paper and Supplemental Material are in `paper/`; their PDFs are in `output/pdf/`. `paper/submission-manifest.json` binds the submission files, and `paper/revision-manifest.json` retains the research-input hashes.

With NumPy installed, `python verification/replay_public_inference.py` replays the six-coefficient omnibus and physical lag sections using the public `outputs/research-completion/public-inference.json`. This does not require omitted arrays. The broader `verify_unified_paper.py` does require them; the compact replay does not certify timing derivatives or replace the full stellar/timing calculations.

The obsolete `clock-timing-dictionary` and `collapse-theorem-dynamic-visibility` branch histories are retained under `archive/2026-09-28/` tags. The retired `lpaiu/dynamic-chi-observable` tip is already an ancestor of public main. Local working checkouts are unaffected by public branch cleanup.
