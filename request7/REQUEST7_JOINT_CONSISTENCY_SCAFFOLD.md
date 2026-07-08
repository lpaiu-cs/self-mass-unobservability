# Request 7: Provisional Joint Consistency Scaffold

This is the next meaningful local step after the current freezes:

- `Request 4` local `MLRS/PEP` branch stops at external environment hand-off,
- `Request 5 Phase B` has public inputs but not a closed local runtime stack,
- `Request 6` is frozen as a local clock-sector audit.

So the highest-value in-workspace result is not another hand-off memo.
It is a **provisional joint scaffold** that combines the living local branches:

- `Request 5 Phase A` as the strong-field free-fall anchor,
- `Request 6` as the local clock-sector audit,
- and `Request 3` only as a weak-field scale reference, not as evidence.

## Artifacts

- scaffold script:
  [request7_joint_consistency_scaffold.py](request7_joint_consistency_scaffold.py)
- machine summary:
  [request7_joint_consistency_summary.json](request7_joint_consistency_summary.json)
- visual summary:
  [request7_joint_consistency_summary.svg](request7_joint_consistency_summary.svg)
- scenario table:
  [request7_joint_consistency_scenarios.tsv](request7_joint_consistency_scenarios.tsv)
- reference tied posterior:
  [request7_joint_consistency_reference_tied_posterior.tsv](request7_joint_consistency_reference_tied_posterior.tsv)

## Model Structure

### Decoupled model

The decoupled scaffold uses:

- free-fall parameters: `sigma_1`, `sigma_2`
- clock parameters: local basis `eta_*`, `kappa_*`

with likelihood

```math
\mathcal L_{\rm dec}
=
\mathcal L_{\rm J0337\,Phase\,A}(\sigma_1,\sigma_2)
\times
\mathcal L_{\rm clock}(\eta_*,\kappa_*).
```

### Tied model

The tied scaffold enforces

```math
\zeta_i = \sigma_i
```

and then maps into the local clock basis at the `Request 6` expansion point
`s_*`:

```math
\eta_* = \sigma_1 s_* + \sigma_2 s_*^2,
\qquad
\kappa_* = \sigma_1 + 2 s_* \sigma_2 .
```

The tied likelihood is then

```math
\mathcal L_{\rm tied}
=
\mathcal L_{\rm J0337\,Phase\,A}(\sigma_1,\sigma_2)
\times
\mathcal L_{\rm clock}\!\left(\eta_*(\sigma),\kappa_*(\sigma)\right).
```

## Important Scope Constraint

`Request 3` is **not** used as a weak-field data likelihood here.

It is retained only as a scale sanity check, through its mock coefficient

```math
A_{\cos D} \approx 13.1\,\sigma_1\ {\rm m},
```

which is stored as reference metadata in the summary.

That is the correct boundary because `Request 3` is a toy injection-recovery
study, not a real-data weak-field evidence branch.

## Clock Likelihood Surrogate

The clock branch is represented by a local Gaussian surrogate in
`(\eta_*, \kappa_*)` derived from the `Request 6` two-system clock-only
posterior at the current expansion point `s_*`. Two variants are swept:

- `grid_matched`: moments matched directly to the `Request 6` grid posterior.
  Its `kappa_*` width (`~2.0e-2`) inherits the flat `zeta` prior box of that
  grid and therefore **overstates** the clock slope information (the Fisher
  lever-arm audit of the same two systems gives `|kappa_*|_95 ~= 8.42`).
- `eta_only_kappa_free` (**reference variant**): `eta_*` moments from the same
  grid posterior (the data-driven band constraint), `kappa_*` width taken from
  the Fisher lever-arm audit (`sigma = 8.42/1.96 ~= 4.3`), correlation
  dropped. This removes the `zeta`-prior-box artifact from the surrogate.

This is a deliberate approximation.

It is justified here because:

- `Request 6` was already shown to be local-basis dominated,
- the current purpose is a provisional consistency scaffold,
- not a final clock-sector TOA analysis.

A natural future upgrade is to moment-match the final
`B1913 + low-side` combined posterior (whose slope information is
approximately data-driven, `|kappa_*|_95 ~= 4.6e-2`) instead of the two-system
stage; that is left to the next revision because it also requires re-deriving
the expansion point consistently.

## Prior Sensitivity Sweep

The scaffold explicitly sweeps:

- `J0337` bound choice:
  `optimistic` / `conservative`
- EOS surrogate prior:
  `low`, `wide`, `high`
- clock prior width in local basis:
  `tight`, `medium`, `wide`

That sweep is not optional decoration.
It is the correct way to show which directions are data-supported and which
remain prior-dominated.

## Current Numerical Reading

The executed scaffold gives a very consistent picture.

For the reference scenario
`optimistic / wide / medium / eta_only_kappa_free`:

- the pre-clock projected tension is only
  `0.591` Mahalanobis
  (`0.592` under the `grid_matched` surrogate — the tension is carried by the
  `eta_*` direction, so the surrogate choice barely moves it),
- the tied posterior projects to
  `|eta_*|_95 ≈ 1.99e-6` (data-limited direction),
  `|kappa_*|_95 ≈ 5.62e-5` (box-conditional: it inherits the open `Phase A`
  degeneracy ridge and its `sigma` prior box, see the `request5`
  `box_scale_check`),
- while the standalone `Request 6` clock audit is much broader:
  `|eta_*|_95 ≈ 1.48e-3` (data-driven),
  with the clock slope direction essentially unconstrained at the two-system
  stage (Fisher `|kappa_*|_95 ~= 8.42`; the earlier quoted `3.59e-2` was set
  by the `zeta` grid box, not by the data).

So, at the current local-data level, the tied image of the `J0337 Phase A`
posterior sits deep inside the allowed clock local-basis region.
In the data-limited `eta_*` direction the tied projection is narrower than
the standalone clock audit by a factor of about `7.5e2`.

That means the present clock branch is **not** creating a strong inconsistency
with the strong-field free-fall branch.
The current joint scaffold therefore reads primarily as

- `J0337 Phase A` setting the scale,
- `Request 6` remaining compatible but weak,
- and the tied-vs-decoupled Bayes factor being strongly affected by prior
  volume rather than by a sharp data clash.

Across the full sweep (2 surrogates x 2 bounds x 3 EOS priors x 3 clock
priors = 36 scenarios):

- `log10 BF(decoupled/tied)` stays between about `-2.41` and `-0.49`
  (`[-2.41, -0.82]` under `grid_matched`, `[-1.49, -0.49]` under
  `eta_only_kappa_free`),
- the span is driven almost entirely by the chosen clock prior box: the clock
  surrogate likelihood is nearly constant over the `J0337` `sigma` grid, so
  the `J0337` factor cancels in the ratio and the BF reduces to an
  Occam-volume ratio (only a handful of distinct BF values occur across all
  scenarios),
- while the projected pre-clock tension and the `J0337`-averaged clock
  likelihood stay nearly unchanged (the tension is structurally stable
  because the symmetric `J0337` posterior always projects to a near-zero
  mean).

That is exactly the signature of a provisional scaffold in which

- the identifiable tension direction is weak,
- the clock branch is local and broad,
- and model comparison is informative mainly as a prior-sensitivity audit,
  **not** as a data-driven discriminant.

## Current Reading

The scaffold should be read as an **interim scientific result**, not as the
final tied-vs-decoupled verdict.

What it answers:

- given the currently available local branches,
- how compatible is the `Request 6` local clock audit with the strong-field
  `J0337 Phase A` posterior under a tied relation?

What it does **not** answer:

- the final `LLR + J0337 + clock` verdict,
- a real weak-field posterior,
- or a full `J0337` TOA refit.

## Status Line

The most accurate one-line status is now:

> `Request 7` is a provisional joint consistency scaffold that combines the
> living local branches (`J0337 Phase A` + clock local audit) under tied and
> decoupled models, while keeping `Request 3` as scale reference only and
> leaving the final verdict to future external weak-field / full-TOA branches.
