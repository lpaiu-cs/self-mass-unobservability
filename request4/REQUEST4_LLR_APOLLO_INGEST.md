# Request 4: APOLLO Real-Data Ingest Scaffold

## Scope

This stage does **not** claim a real LLR parameter fit yet.

It implements the first real-data step that Request 4 was missing:

- pin one concrete public normal-point release,
- download it into the workspace,
- parse it into a machine-readable table,
- and generate a minimal QC/summary product.

That moves Request 4 from "specification only" to an actual reproducible
real-data ingest scaffold.

## Chosen Public Release

The release pinned here is the public APOLLO legacy-text normal-point archive
served from the APOLLO normal-point page:

- source page: [APOLLO Normal Points](https://tmurphy.physics.ucsd.edu/apollo/norm_pts.html)
- access date pinned in the workspace manifest: `2026-04-22`

The page states that the newer CRD format is mirrored through the NASA CDDIS
archive, but the legacy text batches remain directly downloadable and are much
easier to pin reproducibly in this workspace. The scaffold therefore targets
the legacy APOLLO text batches `a` through `f`, which together cover
`2006-04-07` through `2020-12-27`.

The official APOLLO page also documents the field layout of the legacy normal
point record and provides the reflector/station coding used by these files.

For format context, the relevant ILRS pages are:

- [ILRS formats index](https://ilrs.gsfc.nasa.gov/data_and_products/formats/index.html)
- [Historic ILRS normal-point format intro](https://ilrs.gsfc.nasa.gov/data_and_products/data/npt/npt_formatintro.html)

## Implementation

The downloader/parser lives in
[request4_llr_apollo_ingest.py](request4_llr_apollo_ingest.py#L1).

What it does:

1. downloads the APOLLO source page snapshot,
2. downloads the six public batch files:
   - `group_a`
   - `group_b`
   - `group_c_replacement`
   - `group_d_complete`
   - `group_e_complete`
   - `group_f`
3. stores them under the pinned release directory
   `data/request4_llr/apollo_legacy_text_release_2026-04-22/`,
4. computes `sha256` hashes for reproducibility,
5. parses each fixed-width normal-point line into a structured table,
6. writes a TSV and JSON/SVG QC summary.

The parser converts the two-way light time to a one-way geometric range using

```math
R_{\rm one-way} = \frac{c\,t_{\rm round-trip}}{2}.
```

It also maps APOLLO reflector IDs to:

- `Apollo11`
- `Luna17`
- `Apollo14`
- `Apollo15`
- `Luna21`

## Outputs

- parsed table:
  [request4_llr_apollo_normal_points.tsv](request4_llr_apollo_normal_points.tsv#L1)
- release summary:
  [request4_llr_apollo_release_summary.json](request4_llr_apollo_release_summary.json#L1)
- release summary figure:
  [request4_llr_apollo_release_summary.svg](request4_llr_apollo_release_summary.svg)
- pinned raw release:
  [data/request4_llr/apollo_legacy_text_release_2026-04-22/manifest.json](../data/request4_llr/apollo_legacy_text_release_2026-04-22/manifest.json#L1)

## Current Ingest Result

The executed ingest currently yields:

- total normal points: `3846`
- coverage: `2006-04-07T06:26:30Z` to `2020-12-27T06:58:30Z`
- single station ID throughout the release: `70610`
- reflector counts:
  - `Apollo15`: `1792`
  - `Apollo11`: `752`
  - `Apollo14`: `744`
  - `Luna17`: `300`
  - `Luna21`: `258`
- all current records carry quality grade `B`

## One Archive Inconsistency Already Exposed

The scaffold also surfaced one archive-level mismatch that should be carried
forward explicitly rather than hidden:

- the APOLLO page describes `group_b` as `506` normal points spanning
  `2010-12-01` through `2012-04-06`,
- the file currently downloadable from that page parses to `514` records and
  extends to `2012-04-07T07:44:05Z`.

This does **not** look like a parser failure: the raw file itself contains 514
legacy-format records. So at the current stage the correct interpretation is
that the ingest scaffold has uncovered a public-archive inconsistency between
the page description and the downloadable file.

That is exactly the sort of issue the data-ingest stage is supposed to expose
before any science fit is trusted.

## What This Resolves

This stage resolves one of the explicit blockers listed in
[REQUEST4_LLR_REALDATA.md](REQUEST4_LLR_REALDATA.md#L1):

- the workspace now **does** contain a pinned public LLR normal-point release,
- the workspace now **does** contain a validated ingest/parser scaffold,
- and the release manifest makes the data pinning reproducible.

## What It Still Does Not Resolve

This is still not a real Request 4 fit.

The remaining missing pieces are still the actual estimation stack:

- a GR baseline force model and residual model,
- partial derivatives or a Bayesian estimator,
- station/reflector bias estimation,
- troposphere / Earth-orientation surrogates,
- lunar tide and SRP/thermal nuisance estimation.

So this stage should be read as:

- **real-data ingest achieved**,
- **real-data inference not yet achieved**.

## Practical Consequence

After freezing Request 6, this is the correct next move for the tied-vs-
decoupled program. The next concrete escalation for Request 4 is no longer
"find public data" but

1. attach this parsed release to a baseline GR residual model, and
2. then add the EFT remapping for `delta_SEP`.

That next baseline step is now scaffolded explicitly in
[REQUEST4_LLR_BASELINE_SCAFFOLD.md](REQUEST4_LLR_BASELINE_SCAFFOLD.md#L1).

That scaffold has now been pushed one step further into an actual surrogate
baseline probe; see
[REQUEST4_LLR_BASELINE_FIT.md](REQUEST4_LLR_BASELINE_FIT.md#L1).
