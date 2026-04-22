# Request 3: LLR Mock Injection-Recovery

## Scope

This stage implements the **mock-data** part of Request 3, not the production LLR refit from Request 4.

The code uses:

- a planar Sun-Earth-Moon barycentric integrator,
- Newtonian gravity plus a pairwise EIH-inspired toy `1PN` correction,
- passive factors

```math
f_i = 1 + \sigma_1 s_i + \sigma_2 s_i^2,
```

- daily synthetic normal points over 720 days,
- nuisance surrogates for station/reflector bias, solar-radiation pressure, thermal expansion, lunar Love response, and initial-state-vector error.

The implementation lives in [request3_llr_mock.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request3_llr_mock.py).

## 1. Analytic Benchmark

For the weak-field Earth-Moon sector, define

```math
\Delta_{\rm SEP}
=
\sigma_1 (s_E-s_M) + \sigma_2 (s_E^2-s_M^2).
```

Using the standard Nordtvedt weak-field coefficient, the leading synodic range signal is

```math
\delta r_{\cos D} \approx A_N \Delta_{\rm SEP},
```

with

```math
A_N = \frac{13.1\ {\rm m}}{s_E-s_M}.
```

For the self-gravity fractions used in the script,

```math
s_E = 4.64\times 10^{-10},
\qquad
s_M = 1.88\times 10^{-11},
```

this becomes

```math
\delta r_{\cos D}
\approx
(13.1\ {\rm m})\,\sigma_1
+ (6.32468\times 10^{-9}\ {\rm m})\,\sigma_2.
```

The tiny coefficient of `sigma_2` is the main reason weak-field LLR is almost blind to the quadratic term.

## 2. Mock Model

The dynamics are integrated directly for the `sigma` sector. The nuisance sector is then injected and recovered with a linear template model:

```text
range_mock
=
range_dyn(sigma_1, sigma_2)
+ beta_bias * 1
+ beta_srp * SRP(t)
+ beta_thermal * TH(t)
+ beta_love * K2(t)
+ beta_dx * T_dx(t)
+ beta_dy * T_dy(t)
+ beta_dvx * T_dvx(t)
+ beta_dvy * T_dvy(t)
+ noise.
```

This separation is deliberate:

- `sigma` is tested through the actual toy integrator,
- nuisance parameters are treated as local surrogates that are linearized about the GR baseline,
- the fit then profiles over nuisance coefficients on a `(sigma_1, sigma_2)` grid.

So this is a clean **injection-recovery sensitivity study**, not an attempt to reproduce modern LLR analysis pipelines.

## 3. Numerical Outputs

The script writes:

- summary JSON: [request3_llr_mock_summary.json](/Users/lpaiu/vs/lab/self-mass-unobservability/request3_llr_mock_summary.json)
- residual tables:
  [request3_llr_mock_residuals_null.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request3_llr_mock_residuals_null.tsv),
  [request3_llr_mock_residuals_injected.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request3_llr_mock_residuals_injected.tsv)
- posterior tables:
  [request3_llr_mock_posterior_null.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request3_llr_mock_posterior_null.tsv),
  [request3_llr_mock_posterior_injected.tsv](/Users/lpaiu/vs/lab/self-mass-unobservability/request3_llr_mock_posterior_injected.tsv)
- summary figure: [request3_llr_mock_summary.svg](/Users/lpaiu/vs/lab/self-mass-unobservability/request3_llr_mock_summary.svg)

The two scenarios currently implemented are:

1. `null`: `sigma_1 = sigma_2 = 0`
2. `injected`: `sigma_1 = 5\times 10^{-4}`, `sigma_2 = 0`

with `5 mm` mock normal-point noise in both cases.

## 4. Results

From the executed run:

- Null case:
  `sigma_1 = (1.82 \pm 6.24)\times 10^{-5}`
- Injected case:
  `sigma_1 = (4.42 \pm 0.64)\times 10^{-4}`

So the mock pipeline recovers the injected linear signal at the right order and keeps the null case statistically consistent with zero.

For the synodic `cos D` amplitude:

- analytic benchmark for injected case:
  `6.55 mm`
- toy integrator injected amplitude:
  `5.81 mm`
- recovered fit amplitude:
  `5.45 mm`

This is the right consistency target for Request 3: the toy integrator lands in the same mm-scale regime as the weak-field analytic estimate, while remaining visibly less precise than real LLR.

## 5. What the Posterior Means

The `sigma_2` posterior is broad and partly edge-hugging. That is expected here, not a bug:

- weak-field LLR only sees `sigma_2` through `s_E^2-s_M^2`,
- the corresponding synodic coefficient is only `6.3e-9 m` per unit `sigma_2`,
- so a mock study with cm-to-mm nuisance structure has almost no real leverage on `sigma_2`.

Therefore this Request 3 implementation supports the research note's main point:

- LLR mock data can see `sigma_1`,
- LLR mock data are nearly blind to `sigma_2`,
- strong-field systems must carry the quadratic term.

## 6. Limitations

This code is intentionally not a real LLR estimator.

- The `1PN` force law is a compact toy, not a certified production ephemeris.
- The nuisance sector is surrogate-based.
- The geometry is planar and station/reflector modeling is compressed into basis functions.
- The posterior is for sensitivity exploration, not publication-grade inference.

That is acceptable for Request 3 because Request 4 is exactly where the full LLR pipeline belongs.
