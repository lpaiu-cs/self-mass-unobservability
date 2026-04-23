# Data Directory Guide

The `data/` tree contains mirrored public inputs, vendored third-party code, and bounded hand-off workspaces used by the Request 4 and Request 5 scouts.

## `data/request4_llr`

Assets for the real-data LLR branch:

- `apollo_legacy_text_release_2026-04-22/`: mirrored APOLLO legacy normal-point release used by the ingest scaffold.
- `crd_monthly_ensemble_2026-04-23/`: mirrored public monthly CRD payloads and rejected-payload bookkeeping.
- `ilrs_pivot_scout/`: ILRS support code and pivot-scout artifacts.
- `mlrs_handshake_lab/`: bounded MLRS replay / interface-probe workspace.
- `mlrs_promotion_audit_2026-04-23/`: CRD-to-MLRS promotion audit working tree.
- `pep_hand_off_2026-04-23/`: bounded PEP acquisition/build hand-off snapshot.
- `skyfield_cache/`: local ephemeris cache, ignored by Git.

## `data/request5_j0337`

Assets for the strong-field `J0337` Phase B scout:

- `nutimo_public_release_2026-04-23/2020/`: mirrored 2020 public Nutimo release.
- `nutimo_public_release_2026-04-23/2025/`: mirrored 2025 public release, including unpacked and build-probe trees.

## Working Assumption

These directories are kept in Git because the memo trail refers to specific mirrored public releases and bounded reproduction probes. If you later want a slimmer public repository, the first candidates to externalize are the mirrored tarballs, unpacked third-party trees, and large derived workspaces; keep manifests and checksums in Git and move the bulky payloads to Zenodo, GitHub Releases, or Git LFS.
