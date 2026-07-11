# Request 10.8d: Turn-Search Robustness Of The Dynamic-SEP Limit

Status: Counterexample candidate. The 10.8 review conceded ONE unprobed
absorber: pulse-number (turn) re-assignment. Our live gates fix turns; the
published pipeline re-searches them, and the static-Delta posterior width is
plausibly wrap-manifold-dominated (F0 saw wrap onset at |Delta| ~ 1e-7).
This note pre-registers a bounded quantitative probe of whether ANY
turn-aliased solution can absorb the dynamic template at the quoted
amplitudes. Committed before either stage is run.

## Physical framing

Turn re-assignment moves residuals in discrete quanta of P_spin = 2733 us.
A wrap-assisted absorber is a smooth parameter excursion large enough to
advance the predicted phase by integer turns somewhere (the turn-alias
lattice, e.g. `delta_f = n / T_span`), whose WRAPPED residual pattern (a
sawtooth with P-quantum teeth) would have to mimic the smooth
oscillatory-growing dynamic template (~1-3 us at the K10 limit). The probe
quantifies the chi2 cost of every lattice solution.

## Stage V (pre-registered): live wrap-model validation (WSL, ~3 evals)

Shift `spinfreq` (fitted index 4) by `delta_f in {0.5, 1.0}/T_span`
cycles/day via the fitted-parameter interface; the predicted analytic
response is a phase ramp `delta_f * (t - t_ref) * P` wrapped per TOA to
(-P/2, P/2] (with a free offset+slope absorbed). PASS iff the measured
residual change matches the analytic wrapped-ramp model to < 5% relative in
W-norm for both shifts; else STOP (offline wrap modeling invalid; 10.8d
must move fully live).

## Stage L (pre-registered): lattice absorption scan (offline, exact wrap model)

Synthetic data `d = res0 + beta_inj * T_beta(tau)`. Candidate solutions:
turn-alias lattice `s_(n1,n2)(t) = wrap( n1/T * (t - t0) * P + n2 * 2/T^2 *
(t - t0)^2 / 2 * P )` (frequency and frequency-derivative aliases;
`n1 in [-6, 6]`, `n2 in [-3, 3]`, excluding (0,0)) plus the full smooth span
(the 10.8b nuisance block AND the signal templates, i.e. the fitter may also
use beta to help absorb). For each lattice point:

```text
chi2(n1, n2) = min_{c, cY, beta'} || w (d - s_(n1,n2) - B0 c - cY T_cY - beta' T_beta) ||^2
```

PASS iff for EVERY nonzero lattice point and for all three test amplitudes
`beta_inj in {u95_fisher(2 d), u95_K10(2 d), u95_K934(18 d)}`:

```text
chi2(n1, n2) - chi2(0, 0) >= 25       (no wrap solution within 5 sigma)
```

and at (0,0) the recovered `beta'` is within 10% of `beta_inj` (sanity).
Report the minimum lattice margin. FAIL in any cell -> the K10 quote reverts
to K934 pending a fully live turn-marginalized analysis.

## Quote-upgrade rule

Status: Counterexample candidate. V PASS and L PASS -> the 10.8c K10 quote
drops its "unprobed turn-reassignment absorber" caveat (the caveat text is
replaced by a pointer to this note's margins). The K934 bracket and Fisher
floor remain recorded. No numerical change to the quoted curve.

## Outputs

Status: Note. `sep_dynamic/turn_search_10_8d.json` (stage V measurements,
lattice margins), updated `joint_fit_summary.md` caveat, updated review note.
