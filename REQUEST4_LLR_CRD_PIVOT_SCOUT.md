# Request 4: CRD Pivot Scout

## Scope

This stage is the first concrete move away from the APOLLO-only bespoke branch
and toward the canonical Request 4 path.

It does **not** claim a weak-field fit. It answers a narrower question:

- what parts of the `CRD/ILRS canonical path` are publicly reachable right now,
- what public software can actually be pinned in the workspace,
- and what can already be parsed without Earthdata-authenticated archive access.

## Why This Stage Exists

The APOLLO-only surrogate branch has now hit its stop rule. That means Request 4
should no longer grow by adding more surrogate terms to the bespoke APOLLO path.

The next healthy move is to pivot toward:

1. the canonical CRD data path, and/or
2. an existing LLR-oriented estimation stack.

This scout is the bridge between those statements and actual reproducible files
in the workspace.

## Implementation

The scout lives in
[request4_llr_crd_pivot_scout.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_crd_pivot_scout.py:1).

It does four things:

1. pins and probes four public ILRS tarballs:
   - CRD v2 sample code
   - CRD v2 test data
   - `orbitNP`
   - `MLRS Lunar Code`
2. probes the canonical CDDIS CRD root and records that it currently requires
   Earthdata-authenticated access,
3. parses the public lunar CRD sample normal-point files bundled with the MLRS
   release,
4. writes a machine-readable scout summary and a TSV of parsed lunar sample
   normal points.

## Outputs

- parsed lunar CRD sample table:
  [request4_llr_crd_lunar_sample_normal_points.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_crd_lunar_sample_normal_points.tsv:1)
- scout summary:
  [request4_llr_crd_pivot_summary.json](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_crd_pivot_summary.json:1)
- scout figure:
  [request4_llr_crd_pivot_summary.svg](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_crd_pivot_summary.svg)
- pinned asset manifest:
  [data/request4_llr/ilrs_pivot_scout/manifest.json](/Users/lpaiu/vs/lab/self-mass-unobservability/data/request4_llr/ilrs_pivot_scout/manifest.json:1)

## Current Result

The scout establishes three practical facts.

First, the canonical CRD archive is real but not anonymously reachable in this
workspace:

- `https://cddis.nasa.gov/archive/slr/data/npt_crd_v2/` currently returns
  `HTTP 401`, consistent with Earthdata Login gating.

Second, the public ILRS support material is accessible and enough to build a
canonical bridge now:

- CRD sample code tarball: public
- CRD test-data tarball: public
- `orbitNP` tarball: public
- `MLRS Lunar Code` tarball: public

Third, the MLRS bundle already contains public lunar CRD examples, and those
can be parsed now without any authenticated archive access.

## Interpretation

The current scout sharpens the Request 4 pivot into a concrete split:

- `CRD canonical ingest` is feasible now and should remain the data-format
  backbone,
- but `CDDIS live data access` needs Earthdata-authenticated download,
- and the nearest public LLR software hand-off point is the `MLRS Lunar Code`
  bundle, which is promising but still requires local build modernization.

That means the next credible Request 4 actions are now narrower and clearer:

1. keep the CRD parser/ingest path alive,
2. decide whether to modernize the public MLRS lunar stack locally,
3. or else connect this CRD ingest layer to a more mature external LLR
   estimator/codebase.

Either way, Request 4 is no longer blocked on "where is the data format?".
It is now blocked on authenticated archive access and/or estimator integration.
