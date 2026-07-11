# Dynamic sensitivity as the residual free-fall observable: a first upper limit on lag-responding equivalence-principle violation from the pulsar triple PSR J0337+1715

**Status:** draft manuscript (Paper B of the two-paper split; dynamic free-fall sector)
**Repository:** `lpaiu-cs/self-mass-unobservability`, branch `lpaiu/minimal-nonlinear-sideband`
**Author:** Juneyoung, Kim

Draft: 2026-07-12. Companion to Paper A (the static finite-family collapse
theorem). Source of record for all numbers: the `request10_external/`
artifact tree and the `notes/REQUEST10_*` request chain on branch
`lpaiu/minimal-nonlinear-sideband` (commit-chain audit in Appendix D).

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
the timing model (11-28% survival per column, versus 0.06% for static
`Delta`), a property we validate end-to-end against the live integrator,
including pulse-number (turn) ambiguity and phase-convention robustness.
A pre-registered, gate-audited pipeline on 12,474 public Nancay TOAs
(2013-2021) yields no detection and a first upper limit on the amplitude of
any lag-responding SEP oscillation:

```text
|delta Delta| < 1.3 x 10^-9   (95% CL, tau_chi = 2 d, worst drive phase),
rising to 5.0 x 10^-9 at tau_chi = 200 d over the window [2, 500] d,
```

with a purely statistical floor of `1.9 x 10^-10` and an ultra-conservative
bracket of `1.2 x 10^-7` under an unresolved-anchor hypothesis that the data
themselves disfavor. In amplitude the dynamic channel probes SEP structure
three orders of magnitude below the static bound's absorbable scale and
eight orders below the raw static lever arm; in coupling units it is
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
11-28%, versus 0.06% for the static column. The degeneracy that makes the
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
response reaches `2.6 x 10^12 us` per unit oscillation amplitude; the
turn-wrap onset at `|Delta| ~ 10^-7` bounds the usable regime and was
mapped). Any `(c_Y, beta, tau_chi, phase)` template is a closed-form linear
combination of the six columns; no kernel or adiabatic assumption enters.

The templates' key structural property, measured against the full nuisance
span of Section 4: per-column survival 11-28%, shared-`tau` template
survival peaking at 43.8% (`tau = 17.7 d`) — three orders of magnitude less
absorbable than the static column (0.06%). The timing model cannot imitate
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
  truncated at relative singular value `10^-3`. Noise scale
  `s = 1.115` (`chi2_red = 1.243`).
- **Signal block**: `[T_cY, T_beta(tau)]` from the measured columns;
  `tau` grid of 62 points on `[2, 500] d`.
- **Detection**: `Z = max |z|` over the full grid, calibrated by 500
  noise simulations through the identical (trials-inclusive) statistic; an
  anti-causal control (the unphysical advanced-response template,
  `g2 -> -g2`) is scanned identically and must stay quiet.
- **Limits**: flat-prior Gaussian-posterior 95% amplitudes `u95(tau)`.
- **Phase marginalization**: the one unresolved convention is a common
  time-origin `t_off`; carrier phases rotate as `w t_off` (independent
  per-carrier phases belong to the excluded arbitrary-projection collapse
  class). The quoted limit is the *worst phase* per `tau` over 24 origins
  spanning the inner period — valid whatever the true phase. The
  marginalization costs only 12%, confirming phase-robust non-degeneracy.

### 4.2 Anchor brackets

The Fisher machinery's absolute scale was audited rather than assumed. A
ratio spectrum `rho_j = sigma_MCMC,j / sigma_Fisher,j` over the twenty core
timing parameters (using the released chain widths) rejects global Fisher
optimism: `rho <= 1` everywhere, with `rho = 0.8-0.94` precisely where the
comparison is convention-clean. The `~10^3` gap between the Fisher and
published widths of the *static* `Delta` direction is therefore
direction-specific — the turn-wrap manifold that inflates the static
posterior demonstrably does not act on the dynamic templates (Section 5).
Per a pre-registered rule the working anchor is `K_dyn = 10` (a safety
floor, not a measurement), and every quoted quantity carries the full
bracket: the statistical Fisher floor, the `K_dyn = 10` headline, and the
`K = 934` ultra-conservative fallback.

## 5. Validation-gate record

Every stage ran under pre-registered, committed rules; each amendment was
itself committed before the stage it governed. The complete record:

| Gate | Question | Result |
| --- | --- | --- |
| F0 | is the static column a true derivative? | half-step linearity `0.0000` at `h = 10^-8`; wrap onset mapped at `~10^-7` |
| A=0 | does the patched integrator reproduce baseline? | `max diff = 0.0 us` (bit-exact) |
| Ramp | time convention + normalization | integral-kernel `1/2` signature, both hypotheses coherent (Sec. 3.2) |
| F2 | dynamic non-degeneracy | survival 11-28% per column; 43.8% shared-`tau` peak |
| F3 | projected sensitivity | anchored `6 x 10^-8 <= 10^-5` threshold (x160 margin) |
| G1 | live null calibration | one Gauss-Newton fit of the real residuals with live integrator evaluations; `|z_GN - z_lin| <= 0.019` (tolerance 0.3) |
| G2a | estimator calibration, detection regime | `4 sigma` injections recovered within `0.06 sigma` (undamped live GN) |
| G2b | estimator calibration, limit amplitude | injections at the conservative limit (~48,000 whitened units) recovered to 0.16-0.23% relative through step-capped live GN; the live model absorbs essentially nothing |
| Turn V | does the interface re-assign pulse numbers? | no: alias-shift response exceeds `P/2` (1522 vs 1366 us) — turns fixed, wrap arithmetic definitional |
| Turn L | can any turn-aliased solution mimic/hide the template? | 14 chi2-viable alternative turn solutions exist per amplitude (integer slips hidden in observing gaps — the phase-connection ambiguity, quantified against the full span), and the estimator is stable across *all* of them: worst `beta` deviation `5.6 x 10^-11` vs `1.8 x 10^-10` tolerance |

Two by-products deserve independent note. First, the turn-fixing behavior of
the interface is now a measured fact, not a documentation claim. Second, the
lattice scan is, to our knowledge, the first quantification of the
gap-hidden pulse-number-slip degeneracy against a full timing-model span for
this data set — and simultaneously a proof that the degeneracy is harmless
for the dynamic observable: pulse-number freedom is chi2-free but orthogonal
to the template in the estimator's measurement space.

## 6. Results

### 6.1 No detection

Across every pipeline mode — zero-phase, phase-marginalized
(trials-corrected `Z = 1.91`, global `p = 0.30`), and dictionary-anchored
(`p = 0.53`) — the carriers are noise-consistent, and the anti-causal
control is quiet (`p = 0.27`). Linear injection-recovery is exact at all
anchors; live-model injection-recovery is exact to sub-sigma (Section 5).

### 6.2 The upper limit

```text
95% CL upper limit on the amplitude of a lag-responding SEP oscillation,
phase-marginalized (worst phase), K_dyn = 10, window tau_chi in [2, 500] d:

  tau_chi        2 d       5 d      18 d      52 d      200 d
  |dDelta| <   1.30e-9   1.36e-9  1.42e-9   1.84e-9   4.95e-9

  statistical (Fisher) floor:     1.9e-10  (tau = 2 d)
  ultra-conservative (K = 934):   1.2e-7   (tau = 2 d)
  zero-phase convention curve:    1.17e-9  (tau = 2 d; recorded)
```

Interpretation: any mechanism — an internal mode, a light dynamical degree
of freedom, a relaxation response of the pulsar's effective self-gravity —
that makes the SEP parameter of this neutron star *oscillate* at the triple
system's carriers with lags between 2 and 500 days is excluded at amplitudes
three orders of magnitude below the static bound, and the exclusion is
robust to the timing model, pulse-number ambiguity, drive phase, and noise
model at the levels recorded in Section 5.

### 6.3 Coupling units and channel complementarity

Under the leading-order dictionary drive (`F = U/c^2` at the pulsar;
dimensionless carrier amplitudes `f_w ~ 10^-11 - 10^-10` with tasc-locked
phases; O(1) geometric factors set to unity), the same fit gives the
Conjectural coupling bound `beta_phys < 14-17` for `tau in [2, 52] d`. The
clock-sector analysis of the same data (the 10.7 chain) bounds its
rate-coupling composite at `0.4-0.6`: the free-fall drive lacks the clock
channel's `1/Omega` integration amplification, so in coupling units the two
channels bound *different composites* of the underlying state data, while in
raw `Delta` amplitude the free-fall channel is deeper by `3 x 10^8`. The two
bounds are complementary faces of the same `(beta, tau_chi)` datum.

## 7. Scope: what is and is not claimed

Claimed: a first, validated, reproducible 95% upper limit on the amplitude
of lag-responding SEP violation of PSR J0337+1715 over
`tau_chi in [2, 500] d`, in the stated conventions, with the full anchor
bracket published.

Not claimed: any detection; any statement below `tau_chi = 2 d` (the
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
`|delta Delta| < 1.3 x 10^-9`. The method (pre-registration with committed
amendments, live-model gates, adversarial controls, bracketed anchors)
travelled two failed template-calibration designs, one integrator segfault,
one refuted absorber hypothesis, and one mis-scaled gate without ever
quoting an unvalidated number; we consider that audit trail as much a result
as the limit.

Three continuations are natural. A second, independent data set (the
published multi-telescope TOAs) would test the limit's pipeline-dependence
and could tighten the anchor bracket from the data side. The sub-2-day
window — where the pole approaches instantaneity — needs the two-quadrature
degeneracy treated as a joint confidence region rather than a profiled
bound. And the anchor question itself (`K_dyn`) reduces to reproducing the
published static posterior within a linearized-plus-turn-manifold model; the
lattice machinery of Section 5 is the starting point.

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
`REQUEST10_8E_PHASE_MARG_PHYSICAL_ANCHOR.md` -> `REQUEST10_8_REVIEW.md`.

## Appendix B: artifacts and reproduction

All quoted numbers reproduce deterministically from
`request10_external/` with stock numpy: the quoted curves from
`scripts/sep_phase_marg_10_8e.py` (headline) and
`scripts/sep_joint_fit_10_8b.py` (zero-phase), gates from
`scripts/sep_feasibility_gates.py` and `scripts/turn_search_10_8d.py`,
anchor diagnosis from `sep_dynamic/anchor_resolution_10_8c.json`. Inputs of
record: `baseline_planetGR.npz`, `finite_jacobian_v2.npy` (+ `jac_v2/`),
`sep_dynamic/col_SEP_D.npz`, `sep_dynamic/sep_dynamic_columns.npz`, and the
integrator patch `sep_dynamic/patches/sepdyn.patch`. Live-model gates ran in
the WSL runtime described in the 10.7a packet (Zenodo 13899771 build).

## Appendix C: gate-number table

F0: `h = 10^-8`, maxdev 340 us, half-step dev 0.0000. A=0: 0.0 us. Ramp:
`a_hat = 0.502 / 0.0219` (predicted 0.5 / 0.0208). F2: best survival 0.438
at `tau = 17.7 d`. F3: anchored `6.1 x 10^-8`. 10.8b fit: `Z = 0.318`,
`p = 0.95`, anti-causal `p = 0.48`, `s = 1.1148`, nuisance rank 70/90.
G1: `|dz| <= 0.019`. G2a: dev `<= 0.062 sigma_F`. G2b: rel. dev
`<= 0.23%`. Turn V: 1522 us > P/2 = 1366 us. Turn L: 14 viable cells;
worst `beta` dev `5.6 x 10^-11` (tol `1.8 x 10^-10`). 10.8c: core `rho`
max 0.84. 10.8e: trials-corrected `Z = 1.91 (p = 0.30)`; worst-phase cost
12%.

## Appendix D: commit-chain audit (pre-registrations and amendments)

Clock chain: `3d9337c` (10.7b pre-reg) -> `92667cf` (10.7c) -> `a88d27f` ->
`e7bf8e9` (10.7d) -> `43a6618` -> `e7f5f5a` (10.7e) -> `aa50c9b`.
Free-fall chain: `ac13406` (10.8 design) -> `69dd1c1` (F0/F1) -> `fa1fb02`
(R5 refuted) -> `2515481` (patch + A=0) -> `b3633cf` (T2/F2/F3) ->
`2282924` (10.8b pre-reg) -> `6e55569` (G2 amendment) -> `426ad2c` (10.8b)
-> `26abfa1` (10.8c pre-reg) -> `99f0ab3` (10.8c) -> `349c70c` (10.8d
pre-reg) -> `5df4870`/`94c04c3` (amendments) -> `14c950d` (10.8d) ->
`3827649` (10.8e pre-reg) -> `3556363` (10.8e). Every amendment precedes
the stage it governs.

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
