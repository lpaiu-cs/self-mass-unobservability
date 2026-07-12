# Dynamic sensitivity as the residual free-fall observable: a first upper limit on lag-responding equivalence-principle violation from the pulsar triple PSR J0337+1715

**Status:** draft manuscript (Paper B of the two-paper split; dynamic free-fall sector)
**Repository:** `lpaiu-cs/self-mass-unobservability`, branch `lpaiu/minimal-nonlinear-sideband`
**Author:** Juneyoung, Kim
**Note:** Draft 2026-07-12, revision 2 (post 10.8f review response: a quadrature-convention error found in independent review was corrected and the full measurement chain re-run; see `notes/REQUEST10_8F_REVIEW_RESPONSE.md`). Companion to Paper A (the static finite-family collapse theorem). Source of record for all numbers: the `request10_external/` artifact tree and the `notes/REQUEST10_*` request chain on branch `lpaiu/minimal-nonlinear-sideband` (commit-chain audit in Appendix D).

---

## Abstract

The companion paper proves that, inside a fixed-order parity-even free-fall
worldline EFT domain, every observable consequence of a compact body's
self-gravity collapses into finitely many *static* sensitivity coordinates —
coordinates that a pulsar-timing orbital fit absorbs essentially completely.
The theorem's boundary map identifies exactly one escape that generates
qualitatively new free-fall data at fixed order: an orbital-timescale
internal state (boundary `A4`), whose salvaged bookkeeping is a finite state
pair `(beta, tau_chi)` — the amplitude and relaxation lag of a *dynamically
responding* sensitivity. This paper instruments that pair directly. We
promote the strong-equivalence-principle (SEP) parameter `Delta` of the
pulsar triple PSR J0337+1715 — statically constrained at the `2 x 10^-6`
level and statically 99.94% degenerate with the timing model — to a
lag-responding transfer sampled at the system's three orbital carriers, and
we construct *exact* dynamic response templates by making `Delta`
time-dependent inside the published nonlinear three-body timing integrator.
The oscillating-`Delta` signature proves structurally non-degenerate with
the timing model (8-27% survival per column, versus 0.06% for static
`Delta`). The analysis is validated end-to-end against the live
integrator: null calibration, detection-regime and limit-amplitude
injection recovery all pass at every anchor *including the headline lag*
(live-refit displacement `<= 0.035 sigma`, template recovery to
`<= 0.13%` at limit amplitude), and the result is robust to pulse-number
(turn) ambiguity over the tested alias lattice and to the drive phase.
(An intermediate run of these gates failed and was diagnosed as a defect
of the gate harness itself — the audit trail, including the failure, is
part of the record; Section 5.) A pre-registered, gate-audited pipeline
on 12,474 public Nancay TOAs (2013-2021) yields no detection and a first
upper limit on the amplitude of any lag-responding SEP oscillation:

```text
|delta Delta| < 1.7 x 10^-9   (tau_chi = 2 d, worst drive phase;
                               95% statistical CL x K_dyn = 10 anchor),
rising to 7.2 x 10^-9 at tau_chi = 200 d and 1.7 x 10^-8 at the 500 d
window edge,
```

with a purely statistical floor of `2.8 x 10^-10`, a full-rank
(no-truncation) bracket of `3.5 x 10^-9`, and an ultra-conservative bracket
of `1.6 x 10^-7` under an unresolved-anchor hypothesis that the data
themselves disfavor. In amplitude the dynamic channel probes SEP structure
about three orders of magnitude below the published static bound for lags
up to ~50 d (two orders at the window edge); in coupling units it is
complementary to (not competitive with) the clock-sector bound derived from
the same data. The result closes the loop the collapse theorem opened: the
only fixed-order free-fall observable the theorem leaves open is now
measured, bounded, and reproducible from committed artifacts.

---

## 1. Introduction

Free-fall tests of the strong equivalence principle ask whether a body's
gravitational self-energy alters the ratio of its gravitational to inertial
mass. The pulsar triple system PSR J0337+1715 provides the cleanest
strong-field laboratory: a millisecond pulsar in a 1.63-day inner orbit,
with an outer white dwarf in a 327-day orbit supplying the polarizing field.
Published analyses bound the static violation parameter `Delta`
(`m_g = (1 + Delta) m_i` for the pulsar) at the few-parts-per-million level.

The companion paper (Paper A) reframes what such static bounds can, even in
principle, contain. Inside a fixed-order (`Delta_op <= 4`), parity-even,
nonspinning, local free-fall worldline EFT domain, the admissible operator
space collapses onto finitely many static sensitivity coordinates; a
timing-model fit repackages those coordinates into its standard parameters.
Concretely, in the present data set the derivative of the residuals with
respect to static `Delta` is enormous — `|d res / d Delta| ~ 3.4 x 10^10 us`
(a `Delta` of `1.8 x 10^-6` corresponds to 61 ms of raw orbital-polarization
signal) — yet 99.94% of that column lies inside the span of the 28 fitted
timing parameters. The static SEP test is, in the theorem's language, a
measurement of the small non-absorbable remainder of one collapsed
coordinate.

Paper A's boundary map also identifies where the collapse *stops*: dropping
assumption `A4` (no orbital-timescale internal state) breaks the theorem's
Y-only reduction, and the sharp salvage (`F-A4+`) shows that what survives
is not an uncontrolled tower but a finite state datum — a response amplitude
and a relaxation time `(beta, tau_chi)` that must be carried separately from
all Wilson coefficients. That datum is the *only* qualitatively new
free-fall observable the fixed-order theory permits. This paper measures it.

Two facts make the measurement more than a reparameterized static test.
First, a *dynamically responding* `Delta(t)` — the sensitivity lag-tracking
the driving potential through a one-state transfer — produces timing
structure (growing orbital-phase modulation and carrier-quadrature shifts)
that the timing model cannot reproduce: the per-column nuisance survival is
8-27%, versus 0.06% for the static column. The degeneracy that makes the
static test hard is *absent* for the dynamic signature. Second, the response
can be computed exactly: `Delta` is implemented natively in the published
nonlinear three-body integrator, and a ten-line, integrity-gated extension
makes it time-dependent, turning the integrator itself into the template
generator — no adiabatic or secular approximation enters the analysis.

The remainder of the paper: Section 2 defines the observable and its
carrier/comparator structure; Section 3 constructs the exact templates;
Section 4 specifies the pre-registered estimation pipeline and its anchor
brackets; Section 5 reports the validation-gate record, including two
measured findings of independent interest (the turn-fixing behavior of the
integrator interface, and a quantification of chi2-free gap-hidden
pulse-number slips); Section 6 gives the results; Sections 7 and 8 state
scope and outlook.

## 2. The dynamic sensitivity observable

### 2.1 One-state transfer

The minimal `A4`-escape model attaches to the pulsar one internal state
`chi` with linear relaxation dynamics driven by an ambient-field invariant
`F(t)`:

```text
tau_chi  d chi/dt + chi = alpha F(t),
Delta_eff(t) = Delta_static + c_Y F(t) + c_chi chi(t).
```

For a drive component at angular frequency `w`, the readout transfer is

```text
G(i w) = c_Y + beta / (1 + i w tau_chi),        beta = alpha c_chi :
```

an instantaneous term plus a single-pole lag. The pole is the payload: its
frequency dependence across distinct carriers is what no finite instantaneous
(derivative-polynomial) comparator can reproduce.

### 2.2 Carriers and comparator pressure

The triple system supplies three nonresonant drive carriers,

```text
Omega_in = 3.856137 rad/d,  Omega_out = 0.019200 rad/d,
Omega_dif = |Omega_in - Omega_out| = 3.836937 rad/d
(ratio 200.8; 1:1, 2:1, 1:2 excluded),
```

and the request-chain counting theorems fix what three samples of `G` can
distinguish: a real shared-coefficient derivative comparator is pressured
through degree 4, and the deliberately generous complex-coefficient
polynomial comparator through degree 1 (`K_required = N + 2`). The carrier
projection is not arbitrary per-carrier freedom: the physical projection
manifold is phase-locked (five real parameters in the six-dimensional
carrier space), and the timing model's finite span realizes precisely such a
manifold — measured rank `<= 5/6` at every threshold, unchanged when the
planet-block Jacobian columns are re-derived at perturbative steps.

### 2.3 Free-fall implementation

In the integrator of record (Nutimo; Zenodo 13899771), `Delta` enters as a
rescaling of the pulsar's pairwise gravitational coupling,
`G_eff(pulsar, j) = (1 + Delta) G`, applied symmetrically in the
accelerations of the pulsar and of both companions — the full nonlinear
three-body-plus-planet polarization dynamics of the published static test.
The dynamic observable promotes this coupling:

```text
Delta(t) = Delta_static
         + Re sum_w [ (c_Y + beta/(1 + i w tau_chi)) d_w e^{i(w t + phi_w)} ],
```

with unit drive (`d_w = 1`) as the convention-clean primary normalization
(Section 6.3 gives the dictionary-anchored coupling version). Because
`Delta` itself is the dimensionless SEP parameter, an amplitude bound on
`delta Delta` is a physical statement independent of the drive model.

## 3. Exact response templates from the nonlinear integrator

### 3.1 The time-dependent-Delta extension

A surgical patch (seven count-asserted edits, one file) makes the coupling
time-dependent at its unique active application site,

```text
G_eff(i, j; t) = Gg[i][j] + A cos(w t + phi)   on pulsar pairs only,
```

read from the environment at each initialization, built in a separate source
tree. The integrity stop-rule holds exactly: at `A = 0` the patched build
reproduces the baseline residuals with `max |diff| = 0.0 us` (bit-exact).

### 3.2 Time-convention calibration and the integral-kernel signature

The integrator's internal time unit was calibrated empirically after two
line-search designs were stopped by their own pre-registered rules (the
response to an oscillating coupling is broadband-dominated, defeating
single-line and naive-spectral tests). The decisive probe is quasi-static: a
ramp drive (period 30,000 d, phase pi/2) fitted against the *measured*
static column. Both unit hypotheses fit with `R^2 = 0.999` (deep-adiabatic
shape degeneracy) but the amplitudes decide, and they decide coherently:
`a_hat(timedays) = 0.502` and `a_hat(days) = 0.0219`, exactly the `1/2` and
`1/(2 x 24.077)` predicted if (i) the integrator time unit is 24.0774 days
and (ii) the secular-dominated response accumulates as the *integral* of
`Delta(t)` (a phase-drift kernel), for which a ramp yields half the
instantaneous-factorization template. One measurement therefore fixes the
convention and validates the patch normalization end-to-end.

### 3.3 Production columns

Six exact response columns were derived — cosine and sine drive quadratures
at each carrier — by central differences in `A` within the verified linear
regime (adaptive amplitudes down to `3.3 x 10^-10`, where the inner-carrier
response reaches `8.5 x 10^10 us` per unit oscillation amplitude; the
turn-wrap onset at `|Delta| ~ 10^-7` bounds the usable regime and was
mapped). Any `(c_Y, beta, tau_chi, phase)` template is a closed-form linear
combination of the six columns; no kernel or adiabatic assumption enters.

Convention of record (fixed in the 10.8f review response after an
independent review caught a sign swap): with `col_s` defined as the response
to the `+pi/2` (`-sin`) drive quadrature, the *causal* (lagging-pole)
template at angular frequency `w` is `g1 col_c - g2 col_s` with
`g1 = 1/(1 + w^2 tau^2)`, `g2 = w tau/(1 + w^2 tau^2)`; the `+g2`
combination is the *advanced* (anticausal) response and serves throughout as
the adversarial control template.

The templates' key structural property, measured against the full nuisance
span of Section 4: per-column survival 8-27%, causal shared-`tau` template
survival peaking at 21.1% (`tau = 5 d`) — a factor of 130-440 less
absorbable than the static column (0.061%). The timing model cannot imitate
an oscillating coupling.

## 4. Pre-registered estimation pipeline

### 4.1 Estimator

All quoted statistics derive from a deterministic, artifacts-only joint
weighted fit (Frisch-Waugh-Lovell block elimination; seed 20260710):

- **Data**: 12,474 public Nancay TOAs (2013-2021, span 2987.9 d), baseline
  residuals of the published maximum-posterior solution (wrms 1.887 us,
  EFAC 1.1225).
- **Nuisance block**: the 28 fitted-parameter response columns (with the
  seven planet-block columns re-derived at perturbative steps after the
  original one-sigma-scale steps were shown to be nonperturbative secants),
  an offset column, thirty low-frequency Fourier pairs (prior-free red-noise
  marginalization), and the static-`Delta` column (static-leakage guard);
  truncated at relative singular value `10^-3` *with the guard column's
  truncated residual retained as an explicit extra nuisance direction*
  (rank 71 of 90; guard residual norm `1.6 x 10^-7`), so static-`Delta`
  leakage cannot re-enter through the truncation. Noise scale `s = 1.1149`
  (`chi2_red = 1.243`), computed per configuration, never reused across
  spans. A full-rank (rank-90, no-truncation) variant of the entire pipeline
  is carried as a bracket on every quoted limit.
- **Signal block**: `[T_cY, T_beta(tau)]` from the measured columns;
  `tau` grid of 65 points on `[2, 500] d`.
- **Detection**: `Z = max |z|` over the full grid, calibrated by 500
  noise simulations through the identical (trials-inclusive) statistic; an
  anti-causal control (the unphysical advanced-response template,
  `g2 -> -g2`) is scanned identically. The pre-registered detection rule
  requires global `p < 0.003` *and* a quiet control.
- **Limits**: flat-prior Gaussian-posterior 95% amplitudes `u95(tau)`.
  These are *statistical* 95% intervals; the anchor factors of Section 4.2
  are systematic safety multipliers applied on top and are not themselves
  probabilistic confidence statements. Quoted limits name their tier
  explicitly.
- **Phase marginalization**: the one unresolved convention is a common
  time-origin `t_off`; carrier phases rotate as `w t_off` (independent
  per-carrier phases belong to the excluded arbitrary-projection collapse
  class). The quoted limit is the *worst phase* per `tau` and per anchor
  tier over the full outer-period domain — 4,821 origins covering
  `[0, P_out = 327.26 d)` in inner-period/24 steps — valid whatever the
  true phase. The marginalization costs 23% at `tau = 2 d`, confirming
  phase-robust non-degeneracy.

### 4.2 Anchor brackets

The Fisher machinery's absolute scale was audited rather than assumed. A
ratio spectrum `rho_j = sigma_MCMC,j / sigma_Fisher,j` over the twenty core
timing parameters (using the released chain widths) rejects global Fisher
optimism: `rho <= 1` everywhere (median 0.005, maximum 0.84), with
`rho = 0.8-0.94` precisely where the comparison is convention-clean. The
`~10^3` gap between the Fisher and published widths of the *static* `Delta`
direction (`K_static = 934.4 = 1.8 x 10^-6 / 1.926 x 10^-9`, both factors
measured) is therefore direction-specific — the turn-wrap manifold that
inflates the static posterior demonstrably does not act on the dynamic
templates (Section 5). Per a pre-registered rule the working anchor is
`K_dyn = 10` (a safety floor, not a measurement), and every quoted quantity
carries the full bracket: the statistical Fisher floor, the `K_dyn = 10`
headline, and the `K = 934` ultra-conservative fallback. The live-refit
displacement measured by gate G1 (`<= 0.035 sigma_F` at every anchor
including the headline lag; Section 5) is negligible against this
bracket.

## 5. Validation-gate record

Every stage ran under pre-registered, committed rules; each amendment was
itself committed before the stage it governed. After an independent review
of the first complete run found a quadrature-convention swap (the advanced
template had been fit as "causal"), the convention was fixed (Section 3.3),
the correction plan was itself committed first
(`notes/REQUEST10_8F_REVIEW_RESPONSE.md`), and every data-touching stage
below was **re-run on the corrected causal branch**. The record below is
the corrected record; verdicts are reported under the original registered
criteria, unchanged.

| Gate | Question | Result |
| --- | --- | --- |
| F0 | is the static column a true derivative? | half-step linearity `0.0000` at `h = 10^-8`; wrap onset mapped at `~10^-7` |
| A=0 | does the patched integrator reproduce baseline? | `max diff = 0.0 us` (bit-exact) |
| Ramp | time convention + normalization | integral-kernel `1/2` signature, both hypotheses coherent (Sec. 3.2) |
| F2 | dynamic non-degeneracy | PASS: survival 8-27% per column; 21.1% causal shared-`tau` peak (threshold 3%) |
| F3 | projected sensitivity | PASS: anchored `6.9 x 10^-8 <= 10^-5` threshold (x146 margin) |
| G1 | live null calibration: does a live-integrator refit of the *null* data leave the estimator where the frozen linear pipeline puts it? | PASS at `tau = 2 / 18 / 52 d`: offsets `-0.002 / -0.015 / -0.035 sigma_F` (tolerance 0.3) — after one recorded harness failure, diagnosed below |
| G2a | estimator calibration, detection regime (`4 sigma_F` injections, undamped live GN) | PASS at `2 / 18 / 52 d`: `+0.009 / +0.056 / +0.057 sigma_F` (tolerance 0.5) |
| G2b | estimator calibration, limit amplitude (injections at the `K = 934`-anchored limit, step-capped live GN) | PASS at `2 / 18 / 52 d`: `-0.11% / -0.02% / +0.13%` relative (tolerance 2%) — the live model absorbs essentially none of the template |
| Turn V | does the interface re-assign pulse numbers? | PASS: alias-shift response exceeds `P/2` (1522 vs 1366 us) — turns fixed, wrap arithmetic definitional |
| Turn L | can any turn-aliased solution mimic/hide the template? | PASS on the tested lattice (`|n1| <= 12`, `|n2| <= 6`, 48 phase/reference combos per cell, tie-broken wrap): 7 chi2-viable alternative turn solutions per amplitude (integer slips hidden in observing gaps), estimator stable across *all* of them — worst `beta` deviation `1.0 x 10^-10` vs `2.1 x 10^-10` tolerance |

The G-gate row deserves its full history, because an intermediate run
*failed* and the failure is part of the record. The first causal-branch
gate run recorded G1 FAIL at both of its anchors (`+0.34 / +0.60
sigma_F`) and G2a FAIL at 52 d; a first adjudication attributed this to
the live refit converging to the unregularized least-squares point. That
attribution was wrong, and a pre-committed amendment
(`sep_gateG_adjudication.json`, amendment G-2) records the correct
decomposition: the gate harness's *measurement span* had omitted the
static-guard residual direction that the C4 correction added to every
pipeline script — measuring the untouched baseline residuals (no live
refit at all) in that defective span reproduces the entire offset, and
the guard-kept span reproduces the pipeline `z_lin` to machine precision.
With the harness completed (and the headline anchor `tau = 2 d` added),
the full battery was re-run under the unchanged registered criteria:
**9/9 PASS**, and the genuine live-refit displacement is
`<= 0.035 sigma_F` — carried at that measured size on Fisher-tier
numbers, negligible against the anchor bracket. The failed run and its
superseded adjudication remain committed; no criterion was rescored at
any point.

Two by-products deserve independent note. First, the turn-fixing behavior of
the interface is now a measured fact, not a documentation claim. Second, the
lattice scan is, to our knowledge, the first quantification of the
gap-hidden pulse-number-slip degeneracy against a full timing-model span for
this data set — and simultaneously a proof that the degeneracy is harmless
for the dynamic observable *on the tested lattice*: pulse-number freedom is
chi2-free but orthogonal to the template in the estimator's measurement
space. Robustness against arbitrary turn re-assignments beyond the lattice
is not claimed.

## 6. Results

### 6.1 No detection

No pipeline mode meets the pre-registered detection rule (global
`p < 0.003` with a quiet control). Zero-phase: `Z = 1.328`, `p = 0.39`
(anti-causal control `Z = 0.598`, `p = 0.87`). Phase-marginalized over the
full origin domain: `Z = 2.285`, `p = 0.26` — and the anti-causal control
is *equally* elevated (`Z = 2.340`, `p = 0.28`), marking the mild elevation
as noise/structure common to both quadrature conventions rather than a
causal-response feature. Dictionary-anchored: `Z = 1.815`, `p = 0.47`
(reported, not gating). Linear injection-recovery is exact at all anchors;
live-model injection-recovery is exact to `<= 0.06 sigma_F` at detection
amplitudes and `<= 0.13%` at limit amplitudes (Section 5).

### 6.2 The upper limit

```text
Upper limit on the amplitude of a lag-responding SEP oscillation:
95% statistical CL x K_dyn = 10 systematic anchor, phase-marginalized
(worst phase), window tau_chi in [2, 500] d:

  tau_chi        2 d       5 d      18 d      52 d      200 d
  |dDelta| <   1.68e-9   1.85e-9  1.98e-9   2.60e-9   7.16e-9

  statistical-only (Fisher) floor:   2.79e-10  (tau = 2 d)
  full-rank (no-truncation) bracket: 3.53e-9   (tau = 2 d)
  ultra-conservative (K = 934):      1.57e-7   (tau = 2 d)
  window edge (tau = 500 d):         1.74e-8   (K_dyn = 10)
  zero-phase convention curve:       1.36e-9   (tau = 2 d; recorded)
```

Interpretation: any mechanism — an internal mode, a light dynamical degree
of freedom, a relaxation response of the pulsar's effective self-gravity —
that makes the SEP parameter of this neutron star *oscillate* at the triple
system's carriers is excluded at amplitudes about three orders of magnitude
below the published static bound for lags up to ~50 d, shrinking to about
two orders at the 500 d window edge; the exclusion is robust to the timing
model, to pulse-number ambiguity on the tested lattice, to the drive phase,
and to the noise model at the levels recorded in Section 5.

### 6.3 Coupling units and channel complementarity

Under the leading-order dictionary drive (`F = U/c^2` at the pulsar;
dimensionless carrier amplitudes `f_w ~ 10^-11 - 10^-10` with tasc-locked
phases; O(1) geometric factors set to unity), the same fit gives the
Conjectural coupling bound `beta_phys < 18-23` for `tau in [2, 52] d`
(rising to 60 at 200 d). The clock-sector analysis of the same data (the
10.7 chain) bounds its rate-coupling composite at `0.4-0.6`: the free-fall
drive lacks the clock channel's `1/Omega` integration amplification, so in
coupling units the two channels bound *different composites* of the
underlying state data, while in raw `Delta` amplitude the free-fall channel
is deeper by `~3 x 10^8`. The two bounds are complementary faces of the
same `(beta, tau_chi)` datum.

## 7. Scope: what is and is not claimed

Claimed: a first, reproducible upper limit — a 95% *statistical* interval
carrying a stated systematic anchor factor — on the amplitude of
lag-responding SEP violation of PSR J0337+1715 over
`tau_chi in [2, 500] d`, in the stated conventions, with the full anchor
bracket (Fisher floor, full-rank bracket, `K = 10` headline, `K = 934`
fallback) published.

Not claimed: any detection; a probabilistic confidence statement for the
`K`-scaled tiers (the anchors are systematic safety factors, not
coverage-calibrated multipliers); live-refit agreement beyond the measured
`<= 0.035 sigma_F` (gate G1, all three anchors including the headline
lag); turn robustness beyond the
tested alias lattice; any statement below `tau_chi = 2 d` (the
`c_Y`-degeneracy grows toward the instantaneous limit; sub-edge rows are
diagnostic); any statement for other systems or other state models (the
one-pole transfer is the minimal `F-A4+` branch; higher-order state dynamics
would need more carriers per the counting boundaries); an absolute
maximum-likelihood red-noise treatment (the prior-free Fourier
marginalization shifted the clock-channel limits by only 6-10%); physical
drive factors beyond leading order.

## 8. Discussion

The measurement completes a loop that began as a slogan. Paper A converts
"matter cannot sense its own mass" into a conditional collapse theorem whose
sharp boundary map does two things at once: it explains why static
free-fall SEP structure is largely invisible (absorbed into the fit), and it
points to the unique fixed-order escape carrying genuinely new data. This
paper builds the instrument for that escape — exact dynamic templates from
the published integrator itself — and returns the first number:
`|delta Delta| < 1.7 x 10^-9`. The method (pre-registration with committed
amendments, live-model gates, adversarial controls, bracketed anchors)
travelled two failed template-calibration designs, one integrator segfault,
one refuted absorber hypothesis, one mis-scaled gate, one
quadrature-convention swap caught by independent review after the first
complete run (its correction and full re-run pre-registered and committed,
10.8f), one live-gate failure on the corrected branch, and one *wrong
diagnosis of that failure* — superseded, under a further pre-committed
amendment, by the measurement that located the defect in the gate harness
itself and by a full re-run that passes at every anchor. Nothing here was
quoted unvalidated; we consider that audit trail as much a result as the
limit.

Three continuations are natural. A second, independent data set (the
published multi-telescope TOAs) would test the limit's pipeline-dependence
and could tighten the anchor bracket from the data side. The sub-2-day
window — where the pole approaches instantaneity — needs the two-quadrature
degeneracy treated as a joint confidence region rather than a profiled
bound. And the anchor question itself (`K_dyn`) reduces to reproducing the
published static posterior within a linearized-plus-turn-manifold model;
the lattice machinery of Section 5 is the starting point.

## Appendix A: request-chain map

`notes/REQUEST10_DYNAMIC_CHI_MULTIFREQUENCY.md` (10.1, transfer/comparator
counting) -> `REQUEST10_TRIPLE_SHARED_TAU_CARRIER_BRIDGE.md` (10.2b) ->
`REQUEST10_TRIPLE_GR_CARRIER_INVENTORY.md` (10.3, three-carrier boundary) ->
`REQUEST10_TRIPLE_PROJECTION_NUISANCE_GATE.md` /
`REQUEST10_PHYSICAL_PROJECTION_MANIFOLD_GATE.md` (10.4/10.5, projection
classes) -> `REQUEST10_NAMED_TIMING_MODEL_PROJECTION_AUDIT.md` (10.6) ->
`REQUEST10_EXTERNAL_NUTIMO_HANDOFF_PACKET.md` + `request10_external/`
(10.7a pilot) -> `REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md` (10.7b) ->
`..._AMENDMENT_10_7C.md` -> `..._PROMOTION_10_7D.md` ->
`..._PHYSICAL_ANCHOR_REDNOISE_10_7E.md` (clock-channel chain) ->
`REQUEST10_8_DYNAMIC_SEP_POLARIZATION_CHANNEL.md` (10.8 design + F-gates) ->
`REQUEST10_8B_REALDATA_DYNAMIC_SEP_FIT.md` -> `REQUEST10_8C_ANCHOR_RESOLUTION.md`
-> `REQUEST10_8D_TURN_SEARCH_ROBUSTNESS.md` ->
`REQUEST10_8E_PHASE_MARG_PHYSICAL_ANCHOR.md` -> `REQUEST10_8_REVIEW.md` ->
`REQUEST10_8F_REVIEW_RESPONSE.md` (independent-review corrections C1-C12 and
the R1-R7 re-run record).

## Appendix B: artifacts and reproduction

All quoted numbers reproduce deterministically from
`request10_external/` with stock numpy (shared estimator core
`scripts/sep_common.py`; fixed seed 20260710): the quoted curves from
`scripts/sep_phase_marg_10_8e.py` (headline) and
`scripts/sep_joint_fit_10_8b.py` (zero-phase), gates from
`scripts/sep_feasibility_gates.py` and `scripts/turn_search_10_8d.py`,
anchor diagnosis from `scripts/anchor_resolution_10_8c.py`
(-> `sep_dynamic/anchor_resolution_10_8c.json`), static sensitivity and
`K_static` from `scripts/sep_static_sensitivity.py`, live-gate inputs and
adjudication from `scripts/make_gateG2_inputs.py` and
`scripts/gateG_adjudicate.py` (-> `sep_dynamic/sep_gateG2.json`,
`sep_dynamic/sep_gateG_adjudication.json`). Inputs of record:
`baseline_planetGR.npz`, `finite_jacobian_v2.npy` (+ `jac_v2/`),
`sep_dynamic/col_SEP_D.npz`, `sep_dynamic/sep_dynamic_columns.npz`, and the
integrator patch `sep_dynamic/patches/sepdyn.patch`. Live-model
measurements ran in the WSL runtime described in the 10.7a packet (Zenodo
13899771 build); the runtime-side scripts are archived verbatim with an
artifact map in `scripts/wsl/` (`README.md` there lists external inputs:
the runtime build, the TOA/parfile pair, and the released MCMC chain whose
consumed values are cached in `sep_dynamic/steps_provenance.json`).

## Appendix C: gate-number table

All rows from the corrected (causal-branch) re-run of record. F0:
`h = 10^-8`, maxdev 340 us, half-step dev 0.0000. A=0: 0.0 us. Ramp:
`a_hat = 0.502 / 0.0219` (predicted 0.5 / 0.0208). F2: per-column survival
0.082-0.270; best causal shared-`tau` survival 0.211 at `tau = 5 d`. F3:
anchored `6.87 x 10^-8` (Fisher `7.3 x 10^-11`). 10.8b fit: `Z = 1.328`,
`p = 0.39`, anti-causal control `Z = 0.598` (`p = 0.87`), `s = 1.1149`,
nuisance rank 71/90 (guard direction kept; guard residual
`1.65 x 10^-7`). Live gates of record (C4-complete harness, anchors
`tau = 2 / 18 / 52 d`): G1 offsets `-0.002 / -0.015 / -0.035 sigma_F`
(tol 0.3); G2a `+0.009 / +0.056 / +0.057 sigma_F` (tol 0.5); G2b
`-0.109% / -0.024% / +0.128%` (tol 2%) — 9/9 PASS. Superseded first
gate run (defective harness span, recorded): G1 `+0.342 / +0.599` FAIL,
G2a `0.346 / 0.619` PASS/FAIL, G2b PASS; span-only delta on `res0`
reproduces those offsets (`+0.344 / +0.607`), guard-kept span reproduces
`z_lin` to `~1e-16`. Turn V: 1522 us > P/2 = 1366 us. Turn L (lattice
`|n1| <= 12`, `|n2| <= 6`, 48 combos/cell): 7 viable cells; worst `beta`
dev `1.01 x 10^-10` (tol `2.07 x 10^-10`). 10.8c: core `rho` median 0.005,
max 0.838; `K_static = 934.4`. 10.8e (4,821 origins): `Z = 2.285`
(`p = 0.26`), control `Z = 2.340` (`p = 0.28`), dictionary-anchored
`Z = 1.815` (`p = 0.47`); worst-phase cost +23% at `tau = 2 d`; physical
coupling `beta_phys < 18.0 / 18.1 / 18.5 / 23.0 / 59.7` at
`tau = 2 / 5 / 18 / 52 / 200 d`.

## Appendix D: commit-chain audit (pre-registrations and amendments)

Clock chain: `3d9337c` (10.7b pre-reg) -> `92667cf` (10.7c) -> `a88d27f` ->
`e7bf8e9` (10.7d) -> `43a6618` -> `e7f5f5a` (10.7e) -> `aa50c9b`.
Free-fall chain: `ac13406` (10.8 design) -> `69dd1c1` (F0/F1) -> `fa1fb02`
(R5 refuted) -> `2515481` (patch + A=0) -> `b3633cf` (T2/F2/F3) ->
`2282924` (10.8b pre-reg) -> `6e55569` (G2 amendment) -> `426ad2c` (10.8b)
-> `26abfa1` (10.8c pre-reg) -> `99f0ab3` (10.8c) -> `349c70c` (10.8d
pre-reg) -> `5df4870`/`94c04c3` (amendments) -> `14c950d` (10.8d) ->
`3827649` (10.8e pre-reg) -> `3556363` (10.8e) -> `52153b3` (10.8f
review-response pre-reg: corrections C1-C12 committed before the re-runs)
-> the R1-R7 correction commits of this revision (recorded in git history
immediately after `52153b3`). Every amendment precedes the stage it
governs.

## References (draft)

1. S. M. Ransom et al., "A millisecond pulsar in a stellar triple system,"
   Nature 505, 520 (2014).
2. A. M. Archibald et al., "Universality of free fall from the orbital
   motion of a pulsar in a stellar triple system," Nature 559, 73 (2018).
3. G. Voisin et al., "An improved test of the strong equivalence principle
   with the pulsar in a triple star system," A&A 638, A24 (2020).
4. G. Voisin et al., data and code release, Zenodo 13899771 (Nutimo; A&A
   2024, doi:10.1051/0004-6361/202452100).
5. T. Damour, "Binary systems as test-beds of gravity theories," in
   Physics of Relativistic Objects in Compact Binaries (2007) — the
   `(Delta, gammabar, betabar)` strong-field parameterization.
6. C. M. Will, "The confrontation between general relativity and
   experiment," Living Rev. Relativity 17, 4 (2014).
7. W. D. Goldberger and I. Z. Rothstein, "An effective field theory of
   gravity for extended objects," Phys. Rev. D 73, 104029 (2006).
8. K. Nordtvedt, "Equivalence principle for massive bodies," Phys. Rev.
   169, 1014 (1968).
9. J. Kim, "Finite-Family Collapse of Free-Fall Self-Energy Couplings: A
   Fixed-Order Worldline-EFT Theorem, its Uniqueness No-Go, and Sharp
   Boundary Escapes" (Paper A, companion), this repository,
   `paper/paper-A-collapse-theorem.md` (branch `main`).
