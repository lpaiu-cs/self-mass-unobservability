# Request 6: Corrected Clock-Sector Leave-One-Out Fit

This stage replaces the earlier toy ansatz with a version that is actually consistent with the clock EFT written in the project note.

## What Was Corrected

The earlier implementation used

```math
\gamma_{\rm obs} = (1+\zeta_1 s+\zeta_2 s^2)\gamma_{\rm GR},
```

which rescales the full Einstein delay. That is **not** the direct consequence of the project EFT

```math
\frac{d\tau_p}{dt}
=
1-\frac{v_p^2}{2c^2}
-\left(1+\eta_p\right)\frac{Gm_c}{rc^2},
\qquad
\eta_p = \zeta_1 s_p + \zeta_2 s_p^2.
```

If only the gravitational redshift part is modified while the special-relativistic Doppler part remains GR-like, then the observable correction is

```math
\gamma_{\rm SMU}
=
\left(\frac{P_b}{2\pi}\right)^{1/3}(T_\odot M)^{2/3}
X_c\left(1+X_c+\eta_p\right)e
=
\gamma_{\rm GR}\left[1+\frac{\eta_p}{1+X_c}\right].
```

That is the model now implemented in code.

## Leave-One-Out Structure

The second correction is methodological.

The old toy fit inserted masses that had themselves been quoted from GR timing solutions. That made the clock-sector test partially circular. The revised code infers masses from **non-gamma PK information only**:

- `PSR B1534+12`: from `dotomega`, `s`, and `r`
- `PSR J0737-3039A/B`: from theory-independent mass ratio `R`, plus `s` and `r`

Only after that leave-one-out mass inference does the code compare the measured `gamma` to the EFT-modified prediction.

This still is not a full timing-package analysis, but it is no longer reusing `gamma` to define the masses that `gamma` is supposed to test.

## Data Sources

- `PSR B1534+12`: Fonseca et al. (2014), [arXiv:1402.4836](https://arxiv.org/abs/1402.4836)
- `PSR J0737-3039A/B`: Kramer et al. (2021), [arXiv:2112.06795](https://arxiv.org/abs/2112.06795)
- The Double Pulsar mass ratio `R = 1.0714(11)` is taken from Kramer et al. (2006) as quoted in the 2021 paper.

## Numerical Result

Running `python3 request6_clock_sector.py` now gives:

```text
observable model: potential_only

clock-only:
|zeta1|_95 = 3.40e-2
|zeta2|_95 = 2.375e-1
eta(sbar=0.134)|_95 = 1.482e-3

tied model:
|zeta1|_95 = 4.70e-5
|zeta2|_95 = 3.567e-4

Bayes factor:
B(clock/tied) ~= 1.60e-1
```

The important point is that the present sample does **not** separately identify `\zeta_1` and `\zeta_2`. The exact per-system compressed observable is

```math
\delta_{\gamma,i} = \frac{\eta(s_i)}{1+X_{c,i}}.
```

Because the current systems have very similar `s_i` and only modestly different `X_{c,i}`, the present likelihood is **well approximated** by a single effective combination

```math
\eta(\bar s) = \zeta_1 \bar s + \zeta_2 \bar s^2
```

near the self-gravity value actually covered by the current systems. So `\eta(\bar s)` is a representative compression of the true identifiable direction, not the exact fundamental observable.

In the present run:

- global `\bar s ~= 0.134`
- spread across systems `~= 7.87e-4`
- relative spread `~= 0.586%`
- slope-information proxy is correspondingly tiny
- clock-only posterior for `eta(\bar s)`:
  `eta(0.134) = (3.83 +- 6.48)e-4`, with `|eta(0.134)|_95 ~= 1.48e-3`

It is cleaner to reparameterize around `s_* = \bar s`:

```math
\eta_* := \zeta_1 s_* + \zeta_2 s_*^2,
\qquad
\kappa_* := \left.\frac{d\eta}{ds}\right|_{s_*} = \zeta_1 + 2 s_* \zeta_2.
```

The present data mostly constrain `\eta_*`, weakly constrain `\kappa_*`, and are nearly blind to curvature beyond that.

So the current sample barely spans the `s_p` direction at all. That is why `(\zeta_1,\zeta_2)` separately remain broad even though the local combination near `\bar s` is much better constrained.

System-by-system GR leave-one-out diagnostics are:

- `PSR B1534+12`: `gamma_obs/gamma_GR - 1 ~= 3.39e-3`
- `PSR J0737-3039A/B`: `gamma_obs/gamma_GR - 1 ~= 1.31e-4`

At the clock-only posterior mean, the pulls are moderate once predictive uncertainty from the non-gamma mass inference is included:

- `PSR B1534+12`: pull `~= +1.42`
- `PSR J0737-3039A/B`: pull `~= -0.31`

## Interpretation

This version succeeds at the two conceptual repairs that mattered most:

- the `gamma` observable now matches the stated clock EFT,
- the masses are inferred in a leave-one-out way from non-gamma observables.

So this is now a valid **예비검사 / preliminary check** of the corrected clock-sector EFT.

It is still not publication-grade for three reasons:

- the current systems probe almost the same `s_p`, so they mainly constrain `eta(\bar s)` rather than the full `(\zeta_1,\zeta_2)` surface,
- source-specific EOS posteriors for the self-gravity fraction are replaced by a surrogate mass-conditioned prior,
- full timing covariances and TOA-level nuisance models are not included.

Therefore the current `clock-only` limits should be read as a corrected phenomenological **preliminary check**, not as a final validation or falsification of the clock-sector principle.

There is one more caveat for the tied model. In the clock-only run, keeping the free-fall sector at GR makes the leave-one-out PK reconstruction internally consistent. In the tied case, importing a strong-field prior from Request 5 while still using GR-form non-gamma PK reconstruction makes the comparison partly hybrid. The current `B(clock/tied)` is therefore not a final self-consistent physical verdict; it is only a provisional comparison between a broad clock-only likelihood and a tied model dominated by the imported strong-field prior.

This point matters for project strategy. Request 6 should now be read as a
**local clock-sector audit** around the currently sampled `s_*`, not as the
place where the tied model lives or dies. The tied hypothesis is a joint
free-fall-plus-clock statement (`\zeta_i = \sigma_i`), so its actual verdict
has to come from joint consistency with the free-fall sector constraints
already being built in the LLR and J0337 branches. Request 6 can sharpen or
weaken the clock-only branch, but it should not by itself be sold as the final
discriminator between tied and decoupled physics.

## Status After `B1913` Plus Both Low-Side Follow-Ups

That next-stage source program has now been carried out in the current staging
framework.

- `PSR B1913+16` was added through the nuisance-aware extension documented in
  `REQUEST6_B1913_COVARIANCE.md`, which moved the original audit from
  `|kappa_*|_95 ~= 8.42` down to `4.88e-2`.
- `PSR J1141-6545` was then added as the physically cleaner low-side surrogate,
  using the independent scintillation inclination for non-`gamma` mass
  inference and keeping the WD-driven `xdot` / spin-orbit sector explicit as a
  caveat.
- `PSR J1906+0746` was added as the more precise but more nuisance-sensitive
  low-side surrogate, with an external-inclination branch and an explicitly
  correlated `xdot`-conditioned branch.

The net result is that the low-side bottleneck survives.

Using the same original local basis point `s* = 0.134196`, the updated
Fisher-style audit now gives

- baseline after `B1913`: `|kappa_*|_95 = 4.877e-2`
- `+ J1141`: `4.874e-2`
- `+ J1906`: `4.864e-2`
- `+ both`: `4.861e-2`

So adding both low-side systems only improves the current staged audit by about
`0.3%` relative to the `B1913` baseline.

This is consistent with the source-specific effective rows found in the current
surrogate implementation:

- `J1141`: genuine low-side source, but `sigma_delta = 1.787e-2` is far too
  broad to drive the slope direction
- `J1906`: better nominal precision with `sigma_delta = 6.163e-3`, but the
  `xdot`-conditioned branch is effectively suppressed by the `gamma`
  likelihood, so the source does not become the hoped-for low-side rescue

At the combined `B1913 + J1141 + J1906` stage, the front-page local-amplitude
summary is still modestly better constrained than the slope direction:

- local basis point: `s_ref = 0.155832`
- `|eta_*|_95 = 2.106e-3`
- `|kappa_*|_95 = 4.628e-2`

So the present Request 6 branch remains what it was already becoming after the
earlier conceptual fixes: a legitimate local clock-sector audit, not the place
where the tied-vs-decoupled question is finally settled.

## What Has To Come Next

The next upgrade path is now narrower and clearer:

1. do not expect more of the same PK-summary additions to rescue the slope
   direction; `B1913 + J1141 + J1906` already show the present ceiling,
2. if Request 6 is to be pushed further, the next meaningful upgrade is not
   “one more source” but a source-specific covariance-aware or TOA-level
   treatment of the low-side nuisance sectors,
3. if that stronger treatment still fails to drive `|\kappa_*|_{95}` toward
   `1e-2`, then Request 6 should be written up as a support section / local
   clock-sector audit rather than as the project's main novelty engine,
4. the actual tied-vs-decoupled verdict should remain in the joint
   free-fall-plus-clock consistency analysis with the LLR and J0337 branches.

## Files

- `request6_clock_sector.py`
- `request6_clock_sector.ipynb`
- `request6_clock_sector.executed.ipynb`
- `request6_clock_sector_summary.json`
- `request6_clock_sector_summary.svg`
- `request6_clock_sector_posterior_clock_only.tsv`
- `request6_clock_sector_posterior_tied_optimistic.tsv`
- `request6_clock_sector_posterior_tied_conservative.tsv`
- `REQUEST6_B1913_COVARIANCE.md`
- `REQUEST6_LOW_SIDE_EXTENSIONS.md`
