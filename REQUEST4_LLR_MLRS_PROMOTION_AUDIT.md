# Request 4 MLRS Promotion Audit

This audit tests whether representative public monthly lunar `CRD v2` cases can
be promoted into the existing `MLRS` replay path without broad surgery.

Pipeline:

1. Fetch monthly `fr2` from the public `EDC` mirror.
2. Split the monthly aggregate with public `crd_split` into single-pass files.
3. Convert a chosen pass from `CRD v2` to `CRD v1` with public `crd_conv`.
4. Apply a bounded `MLRS-compatible` normalization:
   - keep only the first five non-`na` config IDs in `c0`
   - drop unsupported `c6/c7/41/42` lines
   - lowercase the `h1-h4` tags to match the bundled style
5. Stage the normalized pass as `.frcal` in the `MLRS` lab copy and run `bin/ldb_crd`.

The current question is no longer whether the parser can survive public CRD.
The question is whether representative public passes can produce a usable
`93` residual stream downstream of `recalc`.

Success criterion:

- `recalc` passes, and
- the resulting `.frfil` contains nonzero `93` records so that `Poisson` is operating on a real residual stream.

Failure criterion:

- parser/recalc failure, or
- `recalc` passes but produces `93 = 0`, which means the case reaches the replay path only as a structural input and not yet as a usable observation-operator case.

## Current Result

Representative shortlist promotion currently shows two bounded failure signatures.

1. Some monthly `np2` rows do not have a usable public `fr2` payload at the
   matching public full-rate URL.
2. For rows that do have usable public `fr2`, the public monthly file can be
   split and normalized enough to get through parser-level failure and through
   `recalc`, but the resulting promoted case still reaches `Poisson` with
   `93 = 0`, so no usable residual stream is formed.

So the present verdict is:

- `MLRS` public-monthly transferability is real up to the parser/recalc layer,
- but the weak-field posterior path is still not closed for diverse public
  passes,
- because a public or semi-public `raw frd -> usable frcal / 93 stream` layer
  is still missing.
