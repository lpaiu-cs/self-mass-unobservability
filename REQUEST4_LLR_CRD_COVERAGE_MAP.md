# Request 4: CRD Coverage Map And Representative Cases

## Scope

This stage does not build a weak-field posterior.

It does two narrower things on top of the new monthly lunar CRD ensemble:

1. converts the `67` valid monthly CRD files into a file-level coverage map,
2. produces a small representative-case shortlist for the next bounded
   `MLRS` runnable-promotion attempts.

That is the right next layer after the ensemble acquisition, because the
current bottleneck is now **transferability**, not raw sample count.

## Inputs

- monthly file manifest:
  [request4_llr_crd_monthly_ensemble_files.tsv](request4_llr_crd_monthly_ensemble_files.tsv#L1)
- parsed normal-point rows:
  [request4_llr_crd_monthly_ensemble_normal_points.tsv](request4_llr_crd_monthly_ensemble_normal_points.tsv#L1)

## Implementation

The analysis lives in
[request4_llr_crd_coverage_map.py](request4_llr_crd_coverage_map.py#L1).

For each valid monthly file it computes:

- dominant station and station mix,
- row count,
- coverage span,
- median `num_ranges`,
- median `bin_rms_ps`,
- median `bin_pmm_ps`,
- median `return_rate_hz`,
- median one-way range.

It then forms a representative shortlist in two steps:

1. keep the strongest file for each `(dominant_station, target)` combo,
2. from those combo winners, build a greedy diversity shortlist with a hard
   `row_count >= 5` floor.

Very thin combo winners are not promoted into the shortlist; they are emitted
separately as **coverage gaps**.

## Outputs

- file-level coverage table:
  [request4_llr_crd_monthly_coverage.tsv](request4_llr_crd_monthly_coverage.tsv#L1)
- representative-case shortlist:
  [request4_llr_crd_representative_cases.tsv](request4_llr_crd_representative_cases.tsv#L1)
- summary:
  [request4_llr_crd_coverage_summary.json](request4_llr_crd_coverage_summary.json#L1)
- figure:
  [request4_llr_crd_coverage_summary.svg](request4_llr_crd_coverage_summary.svg)

## Result

The file-level map confirms that the widened CRD ensemble is not homogeneous.

The dominant-station structure is sharply uneven:

- `GRSM` dominates `447` normal points,
- `APOL` contributes `324`,
- `WETL` contributes `171`,
- `MATM` contributes only `19`.

The target structure is also uneven:

- `apollo15` is the clear anchor with `441` normal points across `14` files,
- `luna21` has `147` points across `11` files,
- `apollo11`, `apollo14`, `luna17` sit in the `104–118` range,
- `nglr1` is sparse with only `40` points across `8` files.

That unevenness matters because it means a naive weak-field likelihood would
mostly inherit the geometry and station mix of a small number of monthly files.

## Representative Shortlist

The current runnable-promotion shortlist is:

1. `apollo15_202503.np2` `[WETL/apollo15]` rows=`123`
2. `luna21_202506.np2` `[GRSM/luna21]` rows=`9`
3. `apollo15_202511.np2` `[APOL/apollo15]` rows=`32`
4. `luna17_202502.np2` `[GRSM/luna17]` rows=`51`
5. `apollo15_202502.np2` `[GRSM/apollo15]` rows=`99`
6. `luna21_202502.np2` `[WETL/luna21]` rows=`71`
7. `nglr1_202503.np2` `[APOL/nglr1]` rows=`11`
8. `apollo11_202502.np2` `[GRSM/apollo11]` rows=`45`

This shortlist should be read narrowly:

- it is **not** an insertion-result table,
- it is a bounded list of next candidate months for attempted runnable
  promotion and operator-gain measurement.

## Coverage Gaps

Two dominant-station/target combos are present only in too-thin form and are
therefore flagged as gaps rather than shortlist members:

- `MATM / apollo15`: `apollo15_202601.np2`, rows=`2`
- `APOL / luna17`: `luna17_202511.np2`, rows=`3`

These are useful because they identify where the current monthly archive slice
is still too thin to support a stable runnable-promotion experiment.

## Interpretation

This is the first stage that makes the transferability problem concrete.

The project now has:

- a bounded `MLRS` observation-operator path that works on two bundled sample
  cases,
- and a much wider public CRD observation ensemble that shows how uneven the
  real station/target/month landscape is.

So the next Request 4 move should no longer be "try a posterior anyway".
It should be:

1. take a few cases from this shortlist,
2. see whether they can be promoted to bounded `MLRS` runnable cases,
3. measure `operator gain` there,
4. and only then decide whether a narrow weak-field likelihood is honest.

If the promoted cases stay explainable by a low-dimensional
station/target/time/geometry response layer, Request 4 remains a live main
weak-field branch. If not, the right outcome is to stop at the bounded
observation-operator methodology and hand off the final estimator to a more
mature external stack.
