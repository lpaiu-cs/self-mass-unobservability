# Request 10.7b: External Joint-Fit Experiment (Shared-tau_chi Real-Data Upper Limit)

Status: Counterexample candidate. This note pre-registers the real external
runtime experiment granted **conditionally** by the Request 10.7a verdict
(`runtime-motivated AND conditional`, no stop rule triggered). It implements
promotion conditions C1-C2 of `../request10_external/decision_summary.md`.

Status: Note. Pre-registration discipline: this file was written on 2026-07-10
**before** the first dynamic-chi fit against the real residual vector. Only
baseline summary statistics (wrms, EFAC, drift structure) and the 10.7a
synthetic-injection artifacts were inspected beforehand. No real-data
`beta`-hat value was computed before this note was committed to disk.

## Target And Data

Status: Imported from prior work. Transfer under test (shared across the three
carriers):

```text
G(z) = c_Y + beta / (1 + tau_chi z),   z_k = i Omega_k,
carriers: Omega_in, Omega_out, |Omega_in - Omega_out| of PSR J0337+1715
          (P_i = 1.6293990 d, P_o = 327.25513 d; Nancay 2013-2021).
```

Status: Imported from prior work. Data vector: the real baseline residuals
`res0` (12474 TOAs, us), TOA errors `errs` (us, EFAC 1.1225 fixed as in the
release), internal epochs `toas` (days) of the 10.7a baseline
`parfile-planetGR-max-bestfit` run. Unit-drive convention `Lambda_k F_k = 1`
(pilot caveat (c)): `beta` and `c_Y` are quoted in us of common pre-transfer
drive amplitude; physical drive amplitudes rescale by known factors.

## Pre-Registered Method

Status: Counterexample candidate. Joint linearized model (weighted metric
`W = diag(1/errs^2)`):

```text
r(t) = J dtheta + a0*1 + sum_k Re[(c_Y + beta/(1 + i Omega_k tau_chi)) e^{i Omega_k t}] + n
```

- `J` = the 10.7a central-difference finite Jacobian (12474 x 28, us per
  parameter unit), column-normalized; `a0*1` = explicit constant column
  (offset guard). Nuisance block = 29 columns.
- **Joint fitting (C2):** at each `tau_chi` on the grid, all 31 amplitudes
  (29 nuisance + `c_Y` + `beta`) are fit simultaneously by weighted least
  squares. Implementation uses Frisch-Waugh-Lovell block elimination (project
  BOTH the data and the signal columns off the nuisance span, then solve the
  2x2 problem); this is algebraically identical to the joint solve and is NOT
  the sequential project-then-fit that biased tau-hat toward ~112 d in 10.7a
  Stage 3.
- `tau_chi` grid: quoted window `[10, 327] d` (C1), 49 log-spaced points plus
  anchors `{26, 52, 104} d`. Diagnostic-only extension `[1, 1000] d` is
  computed but not quoted as a limit.
- Noise scale: quoted `sigma_beta` is inflated by `s = sqrt(chi2_red)` of the
  29-parameter nuisance-only GLS fit. **Stop rule:** if `s > 1.75`
  (`chi2_red > 3`), stop and require explicit red-noise modeling before any
  quoted limit.

Status: Counterexample candidate. Detection criterion (must pass BOTH, else
upper limit):

- Global statistic `Z = max |beta-hat/sigma_beta|` over the quoted window;
  null distribution calibrated by >= 500 white-noise simulations
  (`N(0, (s*errs)^2)`) through the identical pipeline; require global
  `p < 0.003`.
- Off-carrier controls: repeat the full fit with all three carriers scaled by
  `f in {0.87, 0.93, 1.07, 1.13}` (no signal expected); require control `|z|`
  consistent with the simulated null.

Status: Counterexample candidate. Upper limit definition (flat-prior Gaussian
posterior on `beta`, marginal over all other amplitudes):

```text
u95(tau): Phi((u - beta_hat)/sigma_beta) - Phi((-u - beta_hat)/sigma_beta) = 0.95
```

(equals 1.96 sigma_beta when beta_hat = 0). Deliverable: the curve
`u95(tau_chi)` over the quoted window, its minimum, and `beta_hat +/- sigma`
at the anchors.

## Pre-Registered Comparator (C2)

Status: Counterexample candidate. The complex degree-1 per-carrier comparator
`A_k = c0 + c1 (i Omega_k)` (4 real parameters) is fit jointly with the same
29-column nuisance block. Per the Request 10.3 counting boundary
(`K_required = N + 2`, three carriers = exactly the `N = 1` boundary), this is
the alternative that three samples cannot decisively separate. Decision rule:
report `chi2(shared-tau, 3 par)` vs `chi2(complex-deg1, 4 par)`; the shared-tau
branch is called *favored* only if it fits at least as well with fewer
parameters. A comparator win does NOT void the upper limit: `u95(tau_chi)`
stands as a bound on the shared-tau branch regardless.

## Pre-Registered Validation

Status: Counterexample candidate.

- Injection-recovery through the FULL joint pipeline: inject
  `beta = u95(tau)` columns at `tau in {26, 52, 104} d` into the real
  residuals; require `|beta-hat - beta_inj| < 2 sigma_beta` at the injected
  `tau` (coverage calibration).
- Sequential-refit bias reproduction (diagnostic): rerun 10.7a-style
  project-then-fit on one injection and confirm the documented tau-hat
  saturation, demonstrating the joint fit removes it.
- Linearization check (WSL runtime, ~3 min/eval): apply the fitted nuisance
  displacement `dtheta-hat` of the `tau = 52 d` configuration to the live
  Nutimo interface; require
  `|| res(p0 + dtheta-hat) - res0 - J dtheta-hat ||_W / || J dtheta-hat ||_W < 0.05`.

## Outputs

Status: Note. Return artifacts, committed next to the 10.7a set under
`../request10_external/`:

| artifact | content |
| --- | --- |
| `joint_fit_upper_limit.json` | full pre-registered pipeline output (limits, detection stats, comparator, validations) |
| `beta_limit_curve.tsv` | `tau_chi, beta_hat, sigma_beta, u95` over quoted + diagnostic grids |
| `joint_fit_summary.md` | one-page verdict using Request 10.7 labels |

## Stop Rules

Status: Proven. Stop (no quoted limit) if any of: `s > 1.75`; a required 10.7a
artifact is missing or inconsistent (carrier frequencies must match
`carrier_projection_rank.json` to 1e-9 rad/day); the injection-recovery
coverage check fails; the linearization check fails.
