# Request 4: MLRS Recalc Seam Probe

## Scope

This stage is the direct follow-on to the bounded `MLRS` interface probe.

It still does **not** implement the real `delta_SEP` physics.
It asks two narrower software-interface questions:

1. does the current recalc seam behave linearly and with odd symmetry under
   small constant perturbations,
2. does a slow synodic-like time-varying perturbation survive the same seam
   without breaking the downstream batch stack.

## Why This Stage Exists

The previous interface probe established that:

- the sample replay mismatch is a real recalc-state microdrift, not just a
  final formatting artifact,
- and a local correction inserted at the recalc `OMC` seam propagates through
  normalpoint generation without breaking the downstream executables.

That was enough to justify one more bounded test before touching the deeper
`jjreadnp.f` / `jeulpkg.f` prediction-state seams.

## Implementation

The probe is implemented in
[request4_llr_mlrs_recalc_seam_probe.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_recalc_seam_probe.py:1).

It does three things:

1. rebuilds the lab `llr_npt_crd` with the current
   [smu_probe_adjust_omc.f](/Users/lpaiu/vs/lab/self-mass-unobservability/data/request4_llr/mlrs_handshake_lab/src/llr_npt/smu_probe_adjust_omc.f:1),
2. runs a constant-offset linearity suite on the first bundled lunar sample
   with offsets
   `{-0.005, -0.002, -0.001, +0.001, +0.002, +0.005} ns`,
3. runs a slow synthetic sine perturbation with:
   - amplitude `0.005 ns`
   - period `29.530588 d`
   - anchor placed between the two bundled lunar sample sessions

The helper now supports both:

- `SMU_OMC_OFFSET_NS`
- `SMU_OMC_SINE_AMPL_NS`, `SMU_OMC_SINE_PERIOD_S`,
  `SMU_OMC_SINE_T0_JD`, `SMU_OMC_SINE_T0_SOD`,
  `SMU_OMC_SINE_PHASE_RAD`

## Outputs

- probe summary:
  [request4_llr_mlrs_recalc_seam_probe_summary.json](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_recalc_seam_probe_summary.json:1)
- probe figure:
  [request4_llr_mlrs_recalc_seam_probe_summary.svg](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_recalc_seam_probe_summary.svg)
- probe script:
  [request4_llr_mlrs_recalc_seam_probe.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_recalc_seam_probe.py:1)

## Result 1: Constant-Probe Linearity Is Clean

For the constant perturbation suite on `s25y11d016t0231#103`:

- all six constant runs change the same `1456 / 1761` `93` residual lines,
- the mean residual shift `Δ93` tracks the injected offset with gain
  `1.0000000000003` to `1.0000000000038`,
- the normal-point `PmM` response scales with gain exactly `1000 ps / ns`,
- the normal-point time-of-flight shift scales with gain
  `1.0000000827e-9` to `1.0000889006e-9 s / ns`,
- `frfin` changes only `12` lines in every constant run,
- `frd` remains an exact match in every constant run,
- `npt` changes only one `11` record in every constant run,
- and `+ε` / `-ε` runs cancel to numerical precision:
  `frrec` odd-symmetry sums stay at `O(10^-15 ns)`, while `npt` TOF and `PmM`
  sums are exactly `0` for `ε = 0.001, 0.002, 0.005 ns`.

So the bounded seam is not only alive. It is also numerically well-behaved as a
small-signal interface.

## Result 2: A Slow Time-Varying Perturbation Also Survives

The synthetic slow sine probe is not the real physics, but it matters for the
next gate.

Using a `29.530588 d` period and `0.005 ns` amplitude:

- the batch chain still runs end-to-end on both bundled lunar sample sessions,
- the perturbation appears in the `93` residual stream on both sessions,
- the downstream filtering/normalpoint stack remains stable,
- `frd` stays exact-match on both sessions,
- `s25y11d016t0231#103` shows a compressed negative response:
  `1456 / 1761` changed `93` lines, mean `Δ93 = -2.415e-3 ns`,
  `npt ΔTOF = -3.000e-12 s`, `npt ΔPmM = -2.4 ps`,
- `s25y11d042t0234#103` shows a smaller positive response:
  `256 / 714` changed `93` lines, mean `Δ93 = +1.250e-3 ns`,
  `npt ΔTOF = +1.000e-12 s`, `npt ΔPmM = 0.0 ps`,
- and the sign/magnitude split across the two sessions confirms that the seam
  can carry a slow time-varying perturbation class rather than only a constant
  bias.

So the seam is not restricted to constant offsets only. It can also carry a
slowly varying perturbation class, although the batch stack clearly compresses
and projects that signal in a sample-dependent way.

## Interpretation

This stage strengthens the `MLRS` reading again, but only in a narrow sense.

What is now justified:

- `MLRS` exposes a recalc-layer seam that behaves like a legitimate
  small-signal interface candidate.

What is **not** yet justified:

- that the real `delta_SEP` physics can be represented as a recalc-layer
  residual tweak,
- or that the final weak-field estimator should be built around this seam.

That distinction matters.

The right reading is therefore:

- this is a bounded software-interface survival test,
- not a physical validation of `delta_SEP`,
- and not yet a claim that the weak-field branch is solved.

The real next gate is one level deeper:

- can the true prediction/state perturbation be inserted at the
  [jjreadnp.f](/Users/lpaiu/vs/lab/self-mass-unobservability/data/request4_llr/mlrs_handshake_lab/src/llr_npt/jjreadnp.f:30)
  / [jeulpkg.f](/Users/lpaiu/vs/lab/self-mass-unobservability/data/request4_llr/mlrs_handshake_lab/src/llr_npt/jeulpkg.f:341)
  seam without broad legacy-code surgery.

If yes, `MLRS` remains a live hand-off path.
If no, this branch should stop here and Request 4 should fall through to a
mature external estimator.
