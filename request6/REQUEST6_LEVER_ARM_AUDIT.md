# Request 6 Follow-up: Lever-Arm Audit

This follow-up does not add a new gravity constraint by itself. It answers the next practical question:

What kind of additional source is needed before the current clock-sector fit can separate the local amplitude `eta_*` from the local slope `kappa_*`?

## Setup

The corrected Request 6 fit showed that the present sample is concentrated near

```math
s_* \simeq 0.1342
```

with only a tiny weighted spread in `s_p`. That means the current data mainly constrain the local combination

```math
\eta_* = \eta(s_*),
\qquad
\kappa_* = \left.\frac{d\eta}{ds}\right|_{s_*},
```

but with strong anisotropy: `eta_*` is seen, `kappa_*` is barely seen.

This audit linearizes the observable

```math
\delta_{\gamma,i} \equiv \frac{\eta(s_i)}{1+X_{c,i}}
\approx
\frac{\eta_* + \kappa_*(s_i-s_*)}{1+X_{c,i}}
```

and uses a Gaussian Fisher estimate to see how the `kappa_*` direction improves when either

- one hypothetical new system is added, or
- a symmetric pair of new systems is added at `s_* \pm |\Delta s|`.

## Current Baseline

Using the present two-system sample:

- `s_* = 0.1341955`
- current linearized Gaussian `|eta_*|_95 ~= 3.31e-3`
- current linearized Gaussian `|kappa_*|_95 ~= 8.42`
- slope-information proxy `~= 4.10e-1`

So the main bottleneck is not sample size by itself. It is the missing `s`-lever arm.

## What The Audit Scans

For a single hypothetical added system, the code scans:

- candidate self-gravity fraction `s in [0.08, 0.20]`
- dimensionless compressed precision `sigma_delta in [5e-5, 2e-3]`
- fixed illustrative `X_c = 0.5`

and computes how much the local slope bound `|kappa_*|_95` would improve.

It also scans a symmetric-pair design, where two new systems are placed at `s_* \pm |\Delta s|` with the same `sigma_delta`.

## Interpretation

The useful output is not a list of pulsars. It is a design rule:

- if a new source lands too close to `s_*`, it hardly helps no matter how many near-duplicate systems are added,
- if a new source is well separated in `s` and still has competitive `sigma_delta`, `kappa_*` improves quickly.
- with the current baseline, **one** additional source is not enough to reach `|kappa_*|_95 < 1e-2` anywhere in the scanned range,
- a **symmetric pair** can reach that regime once both `|\Delta s|` and `sigma_delta` are sufficiently good.

That is exactly the logic behind the weighted spread criterion

```math
\mathcal I_{\rm slope}
\propto
\sum_i
\frac{(s_i-\bar s_w)^2}{(1+X_{c,i})^2 \sigma_{\delta,i}^2}.
```

The audit therefore turns the qualitative claim
"we need genuinely different `s_p`"
into a quantitative target for source selection.

Representative results from the current scan are:

- single candidate, `sigma_delta = 1e-3`, best at `s = 0.200`:
  `|kappa_*|_95 ~= 4.86e-2`
- single candidate, `sigma_delta = 1e-4`, best at `s = 0.200`:
  `|kappa_*|_95 ~= 1.89e-2`

So a single extra source, even at the edge of the scanned `s` range, still does not push the local slope below `1e-2`.

For a symmetric pair:

- pair with `sigma_delta = 2e-4`, best at `|\Delta s| ~= 0.054`:
  `|kappa_*|_95 ~= 7.65e-3`
- pair with `sigma_delta = 1e-4`, best at `|\Delta s| ~= 0.054`:
  `|kappa_*|_95 ~= 3.83e-3`

So the audit says the next useful upgrade is not "one more similar binary." It is either:

1. one very high-leverage outlying system with better-than-current precision, or
2. more realistically, a pair of systems that bracket the current `s_*` with `|\Delta s|` of a few `10^{-2}` and competitive `delta_gamma` precision.

## Files

- `request6_lever_arm_audit.py`
- `request6_lever_arm_audit.ipynb`
- `request6_lever_arm_audit_summary.json`
- `request6_lever_arm_audit_summary.svg`
