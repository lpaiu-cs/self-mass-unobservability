# Data Directory Guide

The `data/` tree contains mirrored public inputs, vendored third-party code, and bounded hand-off workspaces used by the Request 4 and Request 5 scouts.

## `data/request4_llr`

Assets for the real-data LLR branch:

- `apollo_legacy_text_release_2026-04-22/`: mirrored APOLLO legacy normal-point release used by the ingest scaffold.
- `crd_monthly_ensemble_2026-04-23/`: mirrored public monthly CRD payloads and rejected-payload bookkeeping.
- `ilrs_pivot_scout/`: ILRS support code and pivot-scout artifacts.
- `mlrs_handshake_lab/`: bounded MLRS replay / interface-probe workspace.
- `mlrs_promotion_audit_2026-04-23/`: CRD-to-MLRS promotion audit snapshot. The checked-in summary artifacts stay in Git, while the split/work subtree is now treated as a local regenerated workspace.
- `pep_hand_off_2026-04-23/`: bounded PEP acquisition/build hand-off snapshot.
- `skyfield_cache/`: local ephemeris cache, ignored by Git.

## `data/request5_j0337`

Assets for the strong-field `J0337` Phase B scout:

- `nutimo_public_release_2026-04-23/2020/`: mirrored 2020 public Nutimo release.
- `nutimo_public_release_2026-04-23/2025/`: mirrored 2025 public release. The pinned archive and curated top-level mirror stay in Git; transient `unpacked/` and `build_probe/` trees are now ignored.

## Working Assumption

The public repository keeps pinned releases, manifests, memos, scripts, and the source-level patches needed to understand the bounded reproduction story. Re-generated workspaces and extracted third-party duplicates are ignored once they are no longer needed as first-class repo artifacts. If you later want an even slimmer public repository, the next candidates to externalize are the mirrored tarballs and large vendored hand-off labs; keep manifests and checksums in Git and move the bulky payloads to Zenodo, GitHub Releases, or Git LFS.
