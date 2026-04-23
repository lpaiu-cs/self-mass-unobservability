# Request 4: MLRS `frcal/93` Layer Scout

This note closes the last bounded search in the `Request 4` MLRS branch:

- bundled two-case `MLRS` replay already works,
- bounded recalc/state/weak-field seams already work,
- public monthly `CRD` promotion now reaches parser/recalc,
- but promoted public passes still do not produce a usable residual-bearing
  `93` stream.

So the remaining question was narrow:

> Is there a public or semi-public `raw frd -> usable frcal / 93 stream` layer
> that closes the MLRS path into a real weak-field posterior branch?

## Inputs

- Local reproducible scout:
  [request4_llr_mlrs_frcal_layer_scout.py](request4_llr_mlrs_frcal_layer_scout.py)
- Machine summary:
  [request4_llr_mlrs_frcal_layer_summary.json](request4_llr_mlrs_frcal_layer_summary.json)
- Visual summary:
  [request4_llr_mlrs_frcal_layer_summary.svg](request4_llr_mlrs_frcal_layer_summary.svg)
- Upstream promotion evidence:
  [REQUEST4_LLR_MLRS_PROMOTION_AUDIT.md](REQUEST4_LLR_MLRS_PROMOTION_AUDIT.md)

## Evidence Chain

### 1. Official CRD documentation puts calibration/filtering before distribution

The official CRD implementation notes say that one valid station path is:

- convert acquisition data to `CRD`,
- then proceed with calibration, filtering, and normal pointing using `CRD` as
  the station-native format,
- and that this was the path used by `MLRS`.

They also state that some stations use intermediate processing files or
databases that are richer than the distributed format, and only later write the
distributed `CRD` output.

Source:

- [ILRS CRD v1.01 manual](https://ilrs.gsfc.nasa.gov/docs/crd_v1.01.pdf)
- [ILRS CRD overview](https://ilrs.gsfc.nasa.gov/data_and_products/formats/crd.html)

### 2. Official CRD documentation says 9x station records are normally stripped before transmittal

The CRD format examples explicitly note that `91`, `92`, and `93` are
user-defined / station-defined records and that they are normally stripped off
by the station before transmission.

That matters because the current public-monthly promotion failure is exactly:

- parser/recalc survives,
- but the promoted public cases still carry `93 = 0`,
- so `Poisson` never gets the residual-bearing stream seen in the bundled
  reference cases.

Source:

- [ILRS CRD v1.01 manual](https://ilrs.gsfc.nasa.gov/docs/crd_v1.01.pdf)

### 3. The public MLRS bundle expects already-calibrated `.frcal` input

The bundled `MLRS` manual says:

- `Poisson` input is a standard, range-calibrated `CRD` full-rate file,
- `npt_test.sh` uses `.frcal` inputs,
- and those `.frcal` inputs are fully calibrated full-rate files.

So the public `MLRS` test path does not start from generic raw public monthly
`fr2`; it starts from a richer, already-calibrated pass-level input.

Local source:

- [MLRS Lunar Prediction and Normal Point Manual-v1.0.doc](data/request4_llr/mlrs_handshake_lab/MLRS%20Lunar%20Prediction%20and%20Normal%20Point%20Manual-v1.0.doc)

### 4. The public `ldb_crd` script says the calibration code is removed

The shipped `MLRS` batch script still consumes `.frcal`, but the script itself
contains the note:

- `MLRS-specific calibration code removed`

That is exactly the missing layer suggested by the promotion audit.

Local source:

- [ldb_crd](data/request4_llr/mlrs_handshake_lab/bin/ldb_crd)

### 5. Archive access is real, but archive access is not the blocker anymore

The official `CDDIS` guidance says archive download over `https` requires
`Earthdata Login`. That explains the earlier `401` and still matters
operationally.

But it is no longer the core blocker for `Request 4`, because:

- public `ILRS/EDC` monthly `CRD` observation files are already in the
  workspace,
- representative public monthly cases already reach parser/recalc,
- and the failure now occurs later, at the missing residual-bearing
  promotion layer.

Source:

- [CDDIS Archive Access](https://www.earthdata.nasa.gov/centers/cddis-daac/archive-access)

## Bounded Verdict

The bounded search result is now:

- no public or semi-public `raw frd -> usable frcal / 93 stream` layer was
  found,
- the public `MLRS` path remains valid as a bundled two-case observation-operator
  lab,
- but the public-monthly-to-posterior transfer path is still missing the
  residual-bearing promotion layer,
- and the best explanation is that this layer is station-internal, removed from
  the public bundle, or both.

## Project Consequence

`Request 4` stop rule is now triggered.

This means:

- do not keep growing `MLRS` as an open-ended bespoke estimator,
- keep the existing `MLRS` artifacts as evidence that a bounded local
  observation operator exists,
- and hand the final weak-field posterior branch to a more mature external LLR
  estimator / codebase unless a new semi-public calibrated-pass layer appears.

## Current Status Line

The most accurate one-line status is now:

> `MLRS` survives as a realistic bounded weak-field observation-operator
> laboratory, but the public `raw frd -> usable frcal / 93 stream` promotion
> layer is not exposed, so `Request 4` cannot close to a public-MLRS weak-field
> posterior path in the current workspace.
