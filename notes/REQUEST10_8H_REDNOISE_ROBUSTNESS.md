# Request 10.8h: Red-Noise Robustness Of The Free-Fall Limits (Pre-Registration)

Status: Note. This note registers, BEFORE the run, a robustness scan of the
free-fall-channel limits against the red-noise model order. It responds to
second-round review finding B#5 (recorded in
`REQUEST10_8G_REVIEW_RESPONSE_2.md` as out-of-scope for the text round): the
outer carrier `Omega_out` (period `P_out = 327.26 d`) sits inside the
low-frequency band that the prior-free Fourier red-noise block occupies, and
the pole response is increasingly outer-carrier-dominated at large
`tau_chi` — exactly where the quoted limits are weakest. The existing
"6-10%" red-noise sensitivity statement is a clock-channel result and does
not cover the free-fall rows.

## Design facts (computable without data; fixed here)

- The red-noise block is `m_rn` Fourier pairs at frequencies `j/T`,
  `j = 1..m_rn`, `T = t_max - t_min ~= 2987.9 d` (`sep_common.build_nuisance`;
  of record `m_rn = 30`, relative-SV cut `1e-3`, guard kept per C4).
- Outer-carrier proximity: `T / P_out ~= 9.13`, so the `j = 9` mode
  (period 332.0 d) lies 1.4% from the outer carrier. Hence:
  `m_rn = 5` is a below-carrier control (deepest mode 597.6 d, no mode near
  the carrier); every `m_rn >= 10` brackets the carrier; larger `m_rn` adds
  modes below it. The inner carriers (`P ~ 1.63 d`, `j ~ 1834`) are
  untouched by any tested `m_rn`.
- `m_rn = 0` removes the block entirely (no-red-noise control).

## Registered plan

- Harness: `scripts/sep_rn_robustness_10_8h.py` (committed with this note,
  before the run) -> `sep_dynamic/sep_rn_robustness_10_8h.json`.
  It mirrors the 10.8e factorized worst-phase machinery exactly
  (`sep_common`: truncated+guard span, per-span noise scale per C10,
  6-projected-column scan over the full common-origin domain
  `t_off in [0, P_out)` at inner-period/24 steps, per-K worst-phase maxima
  per C3), varying only the red-noise order.
- Grid: `m_rn in {0, 5, 10, 15, 20, 30, 45, 60}` (30 = of record).
- Anchors: `tau_chi in {2, 5, 18, 52, 200, 500} d` (the quoted rows plus the
  window edge). No new tau grid; no full-rank re-scan (the full-rank bracket
  of record is a separate, already-published axis).
- Statistic of record: per anchor and per tier
  (`fisher`, `K10`, `K934`), the worst-phase `u95` on the `m_rn` span, and
  the ratio `R(tau, m) = u95pm_K10(tau; m) / u95pm_K10(tau; 30)`.
  Descriptive additions (no criteria attached): span rank, noise scale
  `s(m)`, worst-phase causal-template survival, and `max|z|(tau, m)`.
- Harness integrity stop-rule: the `m_rn = 30` row must reproduce the
  committed 10.8e values — `u95pm_K10` at the five JSON anchors to relative
  `< 1e-6`, and the `tau = 500 d` tsv row to `< 0.5%` (tsv precision). If it
  does not, the harness is defective: STOP, fix, re-commit; no
  interpretation of any other row.
- Robustness criterion (registered): anchor `tau` is declared
  red-noise-stable iff `max over m in {45, 60} of R(tau, m) <= 1.5`.
  Rationale: only a *richer* noise model weakening the limit threatens the
  quoted rows; poorer-model rows (`m < 30`) and strengthening are reported
  but carry no flag.
- Registered consequences (fixed before the run):
  1. The of-record headline and limit table (m_rn = 30, pre-registered in
     the 10.8b/10.8e chain) do NOT change regardless of outcome. No
     rescoring of any detection verdict; `max|z|` is reported descriptively
     only.
  2. If any anchor fails the criterion, the manuscript quotes, at those
     rows, the worst `u95pm_K10` over `m in {45, 60}` as an explicit
     red-noise bracket (alongside the existing full-rank bracket), and the
     Section 7 scope list drops the clock-channel-derived reassurance in
     favor of the measured free-fall numbers.
  3. If all anchors pass, the manuscript states the measured maximal
     deviation and the outer-carrier/Fourier-grid proximity fact, replacing
     the clock-channel citation.
- Registered expectation (honesty statement, not a criterion): weak `m_rn`
  dependence at `tau = 2-18 d` (inner-carrier-dominated templates); visible
  weakening possible at `tau >= 52 d` where the pole shifts weight to the
  outer carrier; the scan measures its magnitude.
- Out of scope (stated): DM/chromatic noise modeling and an
  enterprise-style GP cross-check (external tooling; listed as future
  work); the sub-2-day window; any change to the estimator or spans other
  than `m_rn`.

## Outcomes (to be recorded post-run; nothing below this line is
## pre-registered except the section's existence)

- (pending)
