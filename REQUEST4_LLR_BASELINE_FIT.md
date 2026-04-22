# Request 4: APOLLO Baseline Surrogate Fit

## Scope

This stage is the first actual residual fit on the pinned APOLLO release, but
it is still **not** a production LLR estimator.

The goal is narrower:

- put a concrete nominal Earth-Moon model under the APOLLO-only branch,
- add the nuisance scaffold from `REQUEST4_LLR_BASELINE_SCAFFOLD.md`,
- and test whether a small bespoke surrogate can plausibly close into a stable
  weak-field branch.

That is exactly the point where the project either justifies a deeper Request 4
build-out or triggers the pivot/stop rule.

## Implementation

The fitter is in
[request4_llr_apollo_baseline_fit.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_apollo_baseline_fit.py:1).

Because the default workstation `astropy/scipy` stack is ABI-broken against the
current NumPy, this stage uses a local runtime:

- local venv: `.venv-request4`
- ephemeris library: `skyfield`
- nominal ephemeris: `DE421`

The nominal model is deliberately limited:

1. compute geocentric Earth-Moon center distance from `DE421`,
2. subtract it from the observed one-way APOLLO range,
3. fit the remaining residuals with a weighted linear surrogate containing:
   - batch indicators,
   - reflector indicators,
   - meteo / observational surrogates,
   - annual / monthly harmonic blocks,
   - solar-day and sidereal-day harmonic blocks,
   - annual sidebands on the main lunar harmonics.

This is still far short of a real LLR stack. It does **not** include full
station geometry, topocentric light time, reflector orientation, full
troposphere, Earth orientation, or relativistic light-time modeling.

## Outputs

- residual series:
  [request4_llr_apollo_baseline_residuals.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_apollo_baseline_residuals.tsv:1)
- fitted coefficients:
  [request4_llr_apollo_baseline_coefficients.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_apollo_baseline_coefficients.tsv:1)
- summary:
  [request4_llr_apollo_baseline_fit_summary.json](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_apollo_baseline_fit_summary.json:1)
- summary figure:
  [request4_llr_apollo_baseline_fit_summary.svg](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_apollo_baseline_fit_summary.svg)

## Current Result

On the pinned APOLLO release, the current surrogate produces:

- nominal `DE421` geocentric Earth-Moon-center residual:
  - `rms = 6.0039e6 m`
  - `wrms = 6.1755e6 m`
- after the APOLLO-only surrogate layer:
  - `rms = 5.0476e5 m`
  - `wrms = 4.0617e5 m`
  - `median |residual| = 2.5733e5 m`
  - `p95 |residual| = 9.7472e5 m`

So the surrogate improves the nominal residual scale by about:

- `11.9x` in RMS
- `15.2x` in WRMS

but it still leaves residuals at `O(10^5-10^6 m)`.

That is the key outcome. The APOLLO-only bespoke branch does **not** close into
a genuinely narrow baseline estimator on this surrogate architecture.

## Interpretation

The fit should be read as a **probe of estimator closability**, not as a
physics verdict.

There are only two outcomes that matter here:

- if the residual layer collapses to something genuinely narrow, Request 4 can
  be grown further inside this bespoke branch;
- if the residual layer stays at large scales despite the nuisance scaffold,
  the correct move is to pivot to an existing LLR estimator or to the
  CRD/ILRS canonical path.

That is the actual decision this stage is meant to inform.

The current result lands in the second branch. It does not say the weak-field
question is impossible. It says this self-built APOLLO-only surrogate is no
longer the healthy place to keep spending effort if the goal is a defensible
final weak-field verdict.
