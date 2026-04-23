# Request 4: Lunar CRD Monthly Ensemble

## Scope

This stage widens the **public real-data CRD ensemble** for Request 4.

It does **not** create 67 new MLRS replay cases. The MLRS bundled end-to-end
replay ensemble is still the same two public sample cases. What changes here is
the size of the observation-side normal-point archive that is already pinned in
the workspace.

## Why This Stage Exists

After the APOLLO-only bespoke branch hit its stop rule, the healthy Request 4
path narrowed to:

1. canonical CRD data ingest, and
2. bounded MLRS hand-off / observation-operator work.

The bundled MLRS release only gives two runnable public sample cases. That was
enough to establish software-interface viability, but not enough for any real
ensemble-level statement about coverage, stations, targets, or sample diversity.

So this stage acquires a larger public lunar CRD normal-point ensemble from the
official EDC mirror.

## Implementation

The acquisition and parse logic lives in
[request4_llr_crd_monthly_ensemble.py](request4_llr_crd_monthly_ensemble.py#L1).

The scan covers:

- six lunar targets: `apollo11`, `apollo14`, `apollo15`, `luna17`, `luna21`,
  `nglr1`
- fifteen target-month candidates per target:
  `2025-01` through `2026-03`

The script now uses the EDC HTTPS mirror directly:

- `https://edc.dgfi.tum.de/pub/slr/data/npt_crd_v2/...`

This choice matters. Direct CDDIS `GET` requests in this workspace can still
fall through to Earthdata login HTML even when `HEAD` returns `200`, whereas
the EDC mirror returns the actual CRD payloads for the public files used here.

The script also validates payload truth rather than trusting `HEAD` alone:

- files whose first bytes do not look like `h1 crd` are treated as rejected
  payloads,
- the tracked fetch registry is
  [data/request4_llr/crd_monthly_ensemble_2026-04-23/fetch_manifest.json](data/request4_llr/crd_monthly_ensemble_2026-04-23/fetch_manifest.json#L1),
- raw and rejected payload directories are kept only as local cache because
  they are publicly downloadable from the upstream mirror.

## Outputs

- candidate/file manifest:
  [request4_llr_crd_monthly_ensemble_files.tsv](request4_llr_crd_monthly_ensemble_files.tsv#L1)
- parsed normal-point table:
  [request4_llr_crd_monthly_ensemble_normal_points.tsv](request4_llr_crd_monthly_ensemble_normal_points.tsv#L1)
- summary:
  [request4_llr_crd_monthly_ensemble_summary.json](request4_llr_crd_monthly_ensemble_summary.json#L1)
- figure:
  [request4_llr_crd_monthly_ensemble_summary.svg](request4_llr_crd_monthly_ensemble_summary.svg)

## Result

The widened public lunar CRD ensemble now contains:

- `67` valid monthly CRD payload files out of `90` target-month candidates
- `961` lunar normal points
- coverage from `2025-01-03T00:37:50Z` to `2026-02-27T22:24:50.366086Z`
- six targets:
  `apollo11`, `apollo14`, `apollo15`, `luna17`, `luna21`, `nglr1`
- four stations:
  `GRSM`, `APOL`, `WETL`, `MATM`

Target-level counts:

- `apollo15`: `14` files, `441` normal points
- `apollo11`: `12` files, `118` normal points
- `apollo14`: `11` files, `104` normal points
- `luna17`: `11` files, `111` normal points
- `luna21`: `11` files, `147` normal points
- `nglr1`: `8` files, `40` normal points

The rejected target-months are informative too. Most of the misses cluster in:

- `2026-03` across every target,
- `2026-01` and `2026-02` for all targets except `apollo15`,
- a small set of isolated 2025 misses such as
  `apollo14_202509`, `luna17_202509`, `luna21_202509`,
  and `nglr1_202501`, `nglr1_202502`, `nglr1_202506`, `nglr1_202507`.

## Interpretation

This materially improves Request 4 on the **data-side ensemble**:

- the workspace no longer has only two lunar public samples total,
- it now has a real monthly CRD archive slice spanning six targets and four
  stations,
- and that is enough to support ensemble-level observation-side calibration and
  sample-selection work.

But the limit is just as important:

- this does **not** mean there are now 67 public MLRS replay cases,
- the bounded MLRS end-to-end replay ensemble is still the same two bundled
  sample cases,
- so the new asset is a widened CRD observation ensemble, not a widened MLRS
  runnable ensemble.

That distinction should remain explicit in all later Request 4 claims.

## Project Meaning

The practical effect is that Request 4 is now stronger in a specific way:

- the weak-field branch now has both a bounded MLRS observation-operator path
  and a broader real lunar CRD normal-point ensemble,
- so the next step no longer needs to argue from only two public sample cases
  on the data side,
- but it still needs to respect that the actual MLRS replay and insertion
  tests remain two-case bounded experiments.

The next healthy use of this ensemble is:

1. operator-level response calibration across more real normal-point months,
2. station/target/time diversity audit for weak-field sensitivity,
3. and only then any attempt at a narrow weak-field likelihood layer.

What it is **not** yet is a final LLR estimator or a published-Nordtvedt-grade
weak-field posterior.
