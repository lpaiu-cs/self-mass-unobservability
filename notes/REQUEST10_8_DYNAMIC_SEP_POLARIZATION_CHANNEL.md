# Request 10.8: Dynamic Free-Fall SEP-Polarization Channel (Design + Feasibility Gates)

Status: Counterexample candidate. The 10.7b-e chain closed the CLOCK-sector
dynamic-chi channel with a first quoted limit whose physical depth saturates
at order unity (`beta_phys < 0.36-0.59`). The structural diagnosis is a double
curse: the clock drive is weak (`e U/c^2 ~ 1e-10`, us-scale Einstein delays),
and the one strong-drive carrier (outer, `D_out = 545 us`) is owned by the
timing model's outer-orbit parameters (absorption cost ~1e3). This note
designs the channel that evades both: the same one-state dynamic-chi transfer
coupled to the FREE-FALL sector, where the static analogue (`Delta`) is
measured at the 1e-6 level on the same data and the polarization signature is
non-degenerate with orbital-element fitting.

Status: Note. Paper-B restructure decision recorded here: the 10.7 clock
chain becomes the methods + first-limit section; the 10.8b real-data outcome
is the headline candidate.

## Theorem-track connection

Status: Imported from prior work. This channel is the observable content of
the `F-A4+` state-augmented salvage: inside the static free-fall theorem
domain every parity-even signal collapses into finitely many static
sensitivity coordinates (already repackaged by the timing fit); at fixed
order the ONLY qualitatively new free-fall observable the theorem permits is
the finite local state data `(beta_ff, tau_chi)` of an orbital-timescale
internal state. 10.8 is the first direct instrumentation of exactly that
data.

## Semantics anchor (verified in the release source)

Status: Proven. Nutimo implements the static parameter natively:
`AllTheories3Bodies.cpp:65` — "In newtonian: M0_gravitational =
(1 + SEP_D) M0_inertial"; lines 958-961 apply `(1 + SEP_D)` to the pulsar's
acceleration AND to the pulsar-as-source terms on both companions, i.e. the
full nonlinear three-body polarization dynamics of the published `Delta`
(Voisin et al. 2020: `|Delta| ~< 2e-6`, 95%) is inside the integrator.
`SEP_D` is parameter index 20, currently 0 and unfitted, parfile scale 1e-6.

## Model

Status: Counterexample candidate. Promote the static `Delta` to the shared
one-state transfer of the 10.x chain, acting in the free-fall sector:

```text
Delta_eff(t) = Delta_static + Re sum_w [ beta_ff/(1 + i w tau_chi) ] d_w e^{i(w t + ph_w)},
```

where `{w}` is the drive-frequency set of the Request-8 invariant at the
pulsar (primary: `Omega_out`; secondary: `Omega_in`, `Omega_dif`; drive
amplitudes `d_w` and phases `ph_w` in the unit-drive and dictionary-anchored
conventions of 10.7e). `Delta_static` (= the `c_Y`-analogue plus the true
static piece) is left free and is NOT the target; the target is the pole pair
`(beta_ff, tau_chi)`.

## Response templates (two routes)

Status: Counterexample candidate. Route T2 (primary — exact): the timing
response is linear in `Delta` at the 1e-6 scale, so exact oscillating
response columns

```text
R_{w,q}(t) = d res / d A   at   SEP_D(t) = A cos(w t + ph_q),  A -> 0
```

are measurable by central differences after a surgical source patch that
makes `SEP_D` time-dependent in the integrator (`dxdt` lines only + parameter
plumbing; patch archived under `request10_external/sep_dynamic/patches/` and
built in a SEPARATE build directory in the external WSL environment,
respecting the 10.7a contract boundary). Any `(beta_ff, tau_chi)` template is
then the linear combination `sum_w Re[G_pole(i w) d_w e^{i ph_w}] R_{w,.}`.
Grid: `w in {Omega_out, Omega_in, Omega_dif}` x two quadratures = 6 columns
(12 evaluations); T2 captures any internal resonance of the forced
polarization response (no adiabatic assumption).

Status: Counterexample candidate. Route T1 (cross-check only — adiabatic):
`template ~= j_Delta(t) * envelope(t)` with `j_Delta = d res/d SEP_D`
(static column), valid only where the drive period far exceeds the
polarization response time; T1-vs-T2 agreement is itself a diagnostic of that
response time.

## Feasibility gates (all data-independent: design-matrix only, no new
real-data statistic; the real-data fit is a separately pre-registered 10.8b)

Status: Counterexample candidate.

- Gate F0 (static column): compute `j_Delta` by central difference at
  `h = 1e-6` (adapt by x10 if maxdev < 0.02 us or /10 if > 30 us); half-step
  linearity < 5%.
- Gate F1 (lever-arm reproduction): the joint-machinery marginal
  `sigma_Delta` (28 params + offset + 30 Fourier pairs, truncated SVD 1e-3,
  noise scale as in 10.7e) must land within a factor 3 of the published
  static sensitivity scale: PASS iff `sigma_Delta <= 5e-6`. This validates
  the channel's 1e-6 lever arm AND our estimator normalization end-to-end
  against the known published result.
- Gate F2 (dynamic non-degeneracy): survival fraction of the T2 dynamic
  templates against the same nuisance span must be >= 3% somewhere in
  `tau_chi in [5, 500] d`; else the channel closes on this dataset
  (collapse-style verdict, reported).
- Gate F3 (projected sensitivity / promotion): the projected
  `sigma_beta_ff(tau)` curve from the joint machinery; promote to the 10.8b
  real-data pre-registration iff `min sigma_beta_ff <= 1e-5` (>= ~3e4 gain
  over the 10.7e clock bound).

## Stop rules

Status: Proven (as rules). Stop if: the T2 patch cannot be confined to the
`SEP_D` application sites plus parameter plumbing (~20 lines); the patched
build at `A = 0` fails to reproduce the unpatched baseline residuals to
max|diff| < 1e-3 us (build-integrity rule); F0 linearity fails at both
adapted steps; F2 fails everywhere (channel closes); any gate would require a
real-data look (all 10.8 gates are design-matrix-only by construction).

## Risk register

Status: Note. R1: the quadrature/lagged polarization may be more degenerate
with orbital elements than the in-phase static signature — that is exactly
what F2 measures, not assumes. R2: the forced-polarization response time may
suppress (or resonantly enhance) the dynamic response — T2 measures it
exactly; T1/T2 disagreement localizes it. R3: patched-build integrity —
covered by the A=0 stop rule. R4: `SEP_D` also rescales the pulsar-as-source
terms; this matches the published `Delta` definition (it is the same
parameter), so F1 compares like with like. R5: `sigma_Delta` may come out
far better than 1.8e-6 (we fix more parameters than the full Voisin
posterior) — F1's factor-3 window is one-sided (an over-tight sigma triggers
a nuisance-completeness review, not a pass).

## Outputs

Status: Note. This note; then `request10_external/sep_dynamic/` with
`jac_sep/col_SEP_D.npz` (F0), `sep_static_sensitivity.json` (F1),
`patches/` + `sep_dynamic_columns.npz` (T2), `sep_feasibility_gates.json`
(F2/F3), and — only if F3 passes — a separate 10.8b pre-registration note
before any real-data fit.

## F0/F1 execution record (2026-07-11)

Status: Proven (F0 — PASS). `j_Delta` computed via parfile edit +
reconstruction (SEP_D is outside the fitted map; no full-parameter shift API
in the interface). Adaptive stepping revealed two facts: (i) the raw lever
arm is enormous — `|d res/d Delta| ~ 3.4e10 us` (a `Delta` of 1e-6 produces
~34 ms residuals); (ii) at `h = 1e-6` and `1e-7` the response saturates at
maxdev ~1.4e3 us ~ P_spin/2 = 1.37 ms — TURN WRAPPING, not orbital
nonlinearity (turns are assigned at construction). In the wrap-free regime
`h <= 1e-8` the response is linear to numerical precision (half-step
deviation 0.0000 at h = 1e-8, maxdev 340 us). The template column is taken
there.

Status: Counterexample candidate (F1 — R5 REVIEW TRIGGERED, as designed).
Against the 10.7e headline nuisance the marginal is
`sigma_Delta = 1.9e-9` (survival 0.1%), one thousand times tighter than the
published `~1.8e-6` — far outside the trust window, so the over-tight branch
fires: DO NOT accept. Working hypothesis: `Delta` rescales the pulsar's
effective GM, whose natural absorbers include parameters FIXED in this
parfile's fit set (`masspar_i`, `mass_o`, `apcosi_i`, possibly `distance`);
our 28-parameter span cannot absorb what the full published analysis
marginalized, leaving a spuriously large "surviving" direction.

Status: Counterexample candidate. R5 resolution plan (next bounded step,
still data-independent): derive central-difference columns for the fixed
physically-degenerate parameters (`masspar_i`, `mass_o`, `apcosi_i`,
`distance`; parfile-edit route, 2 reconstructions each), extend the nuisance
span, re-adjudicate F1. Expected outcomes: `sigma_Delta` relaxes toward the
published scale (lever arm confirmed; proceed to T2), or it stays orders
tighter (then our conditional-sensitivity claim must be explicitly downgraded
to the published marginal before any dynamic-channel projection is quoted).
F2/F3 remain blocked on the R5-resolved span.

## R5 resolution (executed 2026-07-11) — DOWNGRADE branch taken

Status: Proven. Of the four candidate absorber columns, THREE are null in
this parfile's active parametrization: `masspar_i`, `mass_o`, `apcosi_i`
each produce max|d res| = 0.000 us even at 30x-inflated steps (they are
either redundant given the fitted mass/geometry coordinates — the
construction log shows "Parameter set ... : 6" selecting a specific mass
parametrization — or not reached by the value-edit route). Only `distance`
moves residuals (0.211 us at 900 pc). Adding the live `distance` direction to
the 10.7e nuisance span changes the marginal by 0.2% (`sigma_Delta`
1.922e-9 -> 1.926e-9). The fixed-absorber hypothesis is therefore REFUTED:
the 934x gap between our point-estimate Fisher sensitivity (1.9e-9) and the
published marginal (1.8e-6, Voisin+20) is NOT recoverable by adding these
columns. It reflects covariance the full published Bayesian analysis carries
(mass/inclination priors, Kopeikin/parallax terms) that a local
finite-difference span at one fixed parfile cannot reconstruct.

Status: Counterexample candidate. DECISION (pre-registered downgrade branch):
our Fisher `sigma_Delta` is untrustworthy as a physical sensitivity; the
channel's static anchor is set to the PUBLISHED marginal
`sigma_Delta = 1.8e-6` for all downstream projection. What survives from F0/F1
as solid: (i) the lever arm is real and enormous — `|d res/d Delta| ~ 3.4e10 us`,
so `Delta = 1.8e-6` corresponds to ~61 ms of orbital-polarization signal
absorbed by the orbital fit; (ii) the free-fall channel therefore carries
~10^6 x more raw signal per unit coupling than the clock channel's us-scale
Einstein-delay drive. The dynamic projection `sigma_beta_ff` will be quoted
RELATIVE to the 1.8e-6 anchor times the static/dynamic non-degeneracy ratio
measured by F2 (T2 templates). Even at the conservative anchor, a plausible
`sigma_beta_ff ~ 1e-6..1e-5` would still be 10^4-10^5 tighter than the 10.7e
clock bound (`beta_phys ~ 0.4`). F2/F3 require the T2 integrator patch.
