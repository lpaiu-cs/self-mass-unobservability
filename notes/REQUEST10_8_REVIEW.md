# Request 10.8 Chain: Adversarial Self-Review (Referee Pass)

Status: Note. Requested review of the 10.8/10.8b/10.8c chain before Paper B.
Format: the strongest referee attacks we can construct, our answers, and
what must happen before publication-grade claims.

## Result under review

Dynamic free-fall SEP channel on PSR J0337+1715: no detection; quoted 95%
upper limits on the pole-residue amplitude `beta_ff` of an oscillating
`Delta(t)` (window `tau_chi in [2, 500] d`): Fisher floor `1.2e-10`,
re-anchored headline `1.17e-9` (K_dyn = 10, 10.8c), ultra-conservative
bracket `1.09e-7` (K = 934, 10.8b).

## Attacks and answers

- **A1 (turn re-assignment).** "Your live gates fix pulse-number assignments
  at construction; the published pipeline re-searches turns. A dynamic
  template at 300 us amplitudes might be partially absorbable by a different
  turn solution + parameter shifts." — RESOLVED by 10.8d (was conceded).
  Measured: the runtime fixes turns (alias maxdev > P/2); the turn-alias
  lattice scan found 14 chi2-viable alternative turn solutions per amplitude
  (slips hidden in observing gaps — the phase-connection ambiguity,
  quantified), and the estimator's beta-hat is stable across ALL of them
  (worst deviation 5.6e-11 vs 1.8e-10 tolerance; sub-sigma_F). Turn freedom
  cannot mimic or hide the dynamic template. Two documented amendments were
  needed en route (the wrap-validation target was ill-posed; the chi2-margin
  criterion conflated existence with absorption) — both committed before the
  stages they governed.

- **A2 (K_dyn = 10 is a rule artifact).** "The rho spectrum is muddied by
  convention differences (your wider Fourier marginalization), and the rule
  branch fired on semantics designed for rho >= 1." — Partially conceded.
  What is solid: rho <= 1 everywhere refutes global Fisher optimism, and the
  convention-clean directions (abcosi_o 0.84, delta_i 0.94, tasc_extra1 0.8)
  show near-exact calibration. The floor of 10 is a safety factor, not a
  measurement; the full bracket (Fisher / K10 / K934) is published in every
  artifact so no reader is captive to the choice.

- **A3 (drive phase convention).** "Zero drive phase in integrator time is
  arbitrary; a real signal has orbit-locked phase, and your limit is
  per-convention." — True and inherited from the 10.7 chain conventions. For
  an upper limit this is a definitional choice (documented); a
  phase-marginalized two-quadrature estimator (fit both `T(tau)` and its
  quadrature partner) is a mechanical extension planned for the Paper B
  version. Sensitivity loss for a detection would be ~cos(dphi).

- **A4 (unit drive).** "beta_ff at unit drive per carrier is not a physical
  coupling." — Same status as the 10.7e clock quote before E2; the Request-8
  dictionary anchoring (D_k weights, tasc-locked phases) is a mechanical
  re-weighting of the same measured columns. Planned for Paper B; the
  structural conclusions (non-degeneracy, gate record) are unaffected.

- **A5 (noise model).** "White + EFAC + 30 Fourier pairs is not a
  maximum-likelihood red-noise treatment." — Same convention as the quoted
  10.7e headline, where the red-marginalization changed limits by only
  6-10%; chi2_red = 1.24. A spectral-model treatment is a refinement, not a
  blocker.

- **A6 (template basis completeness).** "You test the one-pole shared-tau
  transfer only." — By design: the 10.1-10.3 chain's counting theorems fix
  the comparator classes three carriers can support; the one-state pole is
  the minimal dynamic branch the theorem track's F-A4+ salvage names. Scope
  statement, not a gap.

- **A7 (planet-column linearity at GN caps).** "30x-step GN caps extrapolate
  beyond the half-step linearity checks (0.2-2.5% at h)." — The gates are
  empirical: G2b recovered the injection to 0.16-0.23% THROUGH those
  excursions, which is the quantity of record.

- **A8 (single dataset / single pipeline).** "Nancay-only public release,
  one timing code." — Correct; framed as a first limit and a
  method-validation, not a last word. The entire chain is deterministic and
  artifacts-only reproducible from the committed inputs.

- **A9 (window edges).** "tau < 2 d unexplored; tau -> 0 degenerates with
  c_Y." — Correct and stated; the quoted window edge carries the minimum,
  as in the 10.7 chain, and sub-edge rows are diagnostic.

- **A10 (integral-kernel calibration).** "Your time-unit calibration relied
  on an adiabatic factorization you yourselves showed fails (the 1/2)." —
  The 1/2 IS the integral-kernel prediction, measured at 0.4% agreement
  (0.502), with the alternative hypothesis simultaneously consistent
  (0.0219 vs 0.0208 predicted); production columns are direct measurements
  with no kernel assumption.

## Verdict

Status: Note. Quotable now (with caveats stated): the bracketed limit
`beta_ff < 1.2e-9 (K_dyn = 10) .. 1.1e-7 (K = 934)`, Fisher floor 1.2e-10,
no detection, full gate record, and (10.8d) turn-search robustness RESOLVED
(beta-hat stable across all chi2-viable turn-alias solutions). Remaining
publication items, in order: (1) phase-marginalized estimator variant;
(2) physical drive anchoring (Request-8 weights). Paper B structure stands:
theorem frame (F-A4+) -> methods (10.7 chain + gate discipline) -> headline
(10.8b/c/d) -> loophole map and remaining items as future work.
