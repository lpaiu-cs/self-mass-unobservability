# Request 5 Phase A: PSR J0337+1715 Published-Bound Translation

## Scope

This is the **Phase A** implementation from the repaired plan:

- no raw TOA fit,
- no three-body timing code,
- only translation of published `Delta` bounds into the EFT parameters

```math
\Delta_{\rm SMU}
=
\sigma_1 (s_{\rm NS} - s_{\rm WD})
+ \sigma_2 (s_{\rm NS}^2 - s_{\rm WD}^2).
```

The code uses the two most recent published `J0337+1715` SEP bounds as input hypotheses, from Voisin et al. (2025), A&A 693, A143 ([arXiv:2411.10066](https://arxiv.org/abs/2411.10066)):

- optimistic: `|Delta| < 1.5e-6` at `95%` (circum-ternary planet timing-noise hypothesis)
- conservative: `|Delta| < 2.3e-6` at `95%` (achromatic red-noise timing-noise hypothesis)

and interprets each as a Gaussian upper limit with width `sigma_Delta = Delta_95 / 1.96`. Earlier published limits (Archibald et al. 2018, `|Delta| < 2.6e-6`; Voisin et al. 2020) are compatible with and superseded by these inputs.

Implementation: [request5_j0337_phaseA.py](request5_j0337_phaseA.py)

## Model

The posterior is built on a grid in `(sigma_1, sigma_2)` and marginalized over a neutron-star binding-energy prior:

```math
s_{\rm NS} \sim {\rm Uniform}(s_{\rm min}, s_{\rm max}),
\qquad
s_{\rm WD}=10^{-4}\ \text{fixed}.
```

Three EOS-surrogate priors are compared:

- `low`: `s_NS in [0.10, 0.15]`
- `wide`: `s_NS in [0.10, 0.20]`
- `high`: `s_NS in [0.15, 0.20]`

For a chosen bound model,

```math
{\cal L}(\sigma_1,\sigma_2,s_{\rm NS})
\propto
\exp\left[
-\frac{\Delta_{\rm SMU}(\sigma_1,\sigma_2,s_{\rm NS})^2}{2\sigma_\Delta^2}
\right],
```

and the published-bound posterior is obtained by averaging over the EOS prior.

## Outputs

- summary JSON: [request5_j0337_phaseA_summary.json](request5_j0337_phaseA_summary.json)
- optimistic posterior table: [request5_j0337_phaseA_posterior_optimistic.tsv](request5_j0337_phaseA_posterior_optimistic.tsv)
- conservative posterior table: [request5_j0337_phaseA_posterior_conservative.tsv](request5_j0337_phaseA_posterior_conservative.tsv)
- summary figure: [request5_j0337_phaseA_summary.svg](request5_j0337_phaseA_summary.svg)

## Results

The **data-driven headline numbers are the conditional slices** (wide prior `s_NS in [0.10, 0.20]`):

- optimistic, fixing `sigma_2 = 0`:
  `|sigma_1|_95 ~ 1.15e-5`
- conservative, fixing `sigma_2 = 0`:
  `|sigma_1|_95 ~ 1.75e-5`
- optimistic, fixing `sigma_1 = 0`:
  `|sigma_2|_95 ~ 1.00e-4`
- conservative, fixing `sigma_1 = 0`:
  `|sigma_2|_95 ~ 1.53e-4`

Those conditional values are the right comparison to the back-of-the-envelope statements in the research note. In particular, the optimistic `sigma_1` slice lands at the expected `~1e-5` level.

The fully marginalized `2D` posterior also produces joint quantiles
(`|sigma_1|_95 ~ 4.7e-5`, `|sigma_2|_95 ~ 3.57e-4` for both bounds), but these
are **box-conditional, not standalone constraints**: the published `|Delta|`
bound constrains one EOS-smeared linear combination of `(sigma_1, sigma_2)`,
leaving a slowly decaying degeneracy ridge, and under the flat prior box the
joint marginal quantiles track the prior box. Two facts in the summary JSON
make this explicit:

- the joint quantiles are identical for the optimistic and conservative
  inputs, i.e. they do not respond to the data;
- the `box_scale_check` diagnostic reruns the same posterior with the prior
  box enlarged `x2` and `x4`, and the joint marginals grow by `x3.96/x4.04`
  (essentially linearly with the box), while the conditional slices stay flat
  (`x1.04/x0.93`).

Earlier revisions of this memo quoted the joint marginal quantiles as results;
they should not be used that way.

## EOS Sensitivity

The EOS prior matters, but only at the factor-of-order-unity level in Phase A:

- low prior `[0.10, 0.15]` gives the tightest `sigma_1` slice,
- high prior `[0.15, 0.20]` gives the tightest `sigma_2` slice,
- the wide prior sits between them.

This is exactly what the scaling suggests:

```math
\Delta_{\rm SMU} \sim \sigma_1 s_{\rm NS} + \sigma_2 s_{\rm NS}^2,
```

so larger `s_NS` helps the quadratic term more strongly than the linear one.

## Interpretation

Phase A confirms two points from the repaired plan:

- A published `Delta` bound already gives a usable strong-field posterior without raw TOAs.
- Once both `sigma_1` and `sigma_2` float simultaneously, the result is ridge-like; the clean one-parameter limits require conditional slices or a full Phase B TOA refit.

So this stage is a legitimate **strong-field posterior translator**, but not yet a competitive timing analysis.

## Role In The Project After Freezing Request 6

With Request 6 now effectively frozen as a support/local-audit branch, the
J0337 line becomes one of the main routes to the tied-vs-decoupled verdict.

That means this document should be read in two layers:

- as a Phase A translator from published strong-field bounds into
  `(sigma_1, sigma_2)`,
- and as the current strong-field anchor for the free-fall side of the final
  `LLR + J0337 + clock` consistency program.

In other words, Request 5 is no longer just a side calculation that supports
Request 6. It is now part of the main decision branch.

That main-branch interpretation is now reinforced by a new bounded scout:
`REQUEST5_J0337_PHASEB_PUBLIC_INPUTS.md` shows that public `TOA`, `par`, and
`Nutimo` code releases for `J0337` are in fact available and locally mirrored.
So the remaining gate for Phase B is no longer "do public inputs exist?" but
"can the local workspace close the Nutimo dependency/runtime stack?".

That runtime question has also been checked in bounded form:
`REQUEST5_J0337_PHASEB_BUILD_FEASIBILITY.md` shows that the current host stops
on Apple clang `-fopenmp` and then missing `Boost/Tempo2` stack, so Phase B is
currently public-input-ready but not locally runnable here.

## Practical Consequence

If the project is prioritizing the final tied-vs-decoupled judgement, the next
high-value escalation is not more clock-only source chasing. It is either:

1. a stronger J0337 implementation,
2. a real-data LLR implementation,
3. or the eventual joint combination of those free-fall branches with the now
   fixed Request 6 clock-sector audit.
