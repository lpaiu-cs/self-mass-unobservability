# Data Directory Guide

The `data/` tree contains mirrored public inputs, vendored third-party code, and bounded hand-off workspaces used by the Request 4 and Request 5 scouts.

## `data/request4_llr`

Assets for the real-data LLR branch:

- `apollo_legacy_text_release_2026-04-22/`: APOLLO ingest manifest. Public raw payloads are re-downloadable and ignored by Git.
- `crd_monthly_ensemble_2026-04-23/`: monthly CRD fetch manifest. Public raw and rejected payloads are re-downloadable and ignored by Git.
- `ilrs_pivot_scout/`: pivot-scout manifests and notes. Vendored ILRS downloads are public inputs and are ignored by Git.
- `mlrs_handshake_lab/`: bounded MLRS replay / interface-probe workspace. Public JPL ephemeris mirrors used by the lab stay local and are ignored by Git.
- `mlrs_promotion_audit_2026-04-23/`: CRD-to-MLRS promotion audit snapshot. The checked-in summary artifacts stay in Git, while the split/work subtree is now treated as a local regenerated workspace.
- `pep_hand_off_2026-04-23/`: bounded PEP acquisition/build hand-off snapshot.
- `skyfield_cache/`: local ephemeris cache, ignored by Git.

## `data/request5_j0337`

Assets for the strong-field `J0337` Phase B scout:

- `nutimo_public_release_2026-04-23/`: local mirror root for public Nutimo releases. The payloads are publicly downloadable and ignored by Git; the tracked registry for this scout is [`request5_j0337_phaseB_public_inputs_manifest.tsv`](../request5_j0337_phaseB_public_inputs_manifest.tsv).

## Working Assumption

The public repository keeps manifests, memos, scripts, summaries, and the source-level patches needed to understand the bounded reproduction story. Publicly downloadable payloads are intentionally excluded from Git once their fetch manifests or registries are recorded. Re-generated workspaces and extracted third-party duplicates are also ignored once they are no longer needed as first-class repo artifacts.
