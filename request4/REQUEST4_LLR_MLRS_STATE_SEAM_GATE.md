# Request 4: MLRS State Seam Gate

## Scope

This stage asks the narrow question left open by the recalc-seam probe:

- if the insertion point moves from the recalc `OMC` layer down to the
  `jjreadnp.f` ephemeris-state bridge, does the work still stay bounded,
  or does it immediately widen into broad legacy-code surgery.

It still does **not** implement the real `delta_SEP` physics.
It is a software-architecture gate, not a weak-field fit.

## Implementation

The gate is implemented in
[request4_llr_mlrs_state_seam_gate.py](request4_llr_mlrs_state_seam_gate.py#L1).

The probe adds exactly one new local helper,
[smu_probe_adjust_states.f](../data/request4_llr/mlrs_handshake_lab/src/llr_npt/smu_probe_adjust_states.f#L1),
and one new call site in
[jjreadnp.f](../data/request4_llr/mlrs_handshake_lab/src/llr_npt/jjreadnp.f#L1).

The helper is synthetic and bounded:

- it reads `SMU_STATE_EM_SHIFT_M`,
- converts the requested displacement from meters to AU,
- and applies that shift to the Moon position only, along the instantaneous
  Earth-to-Moon direction, immediately after `PLEPH` has filled
  `TABOUT(:,3)` and `TABOUT(:,11)`.

No changes were made to:

- `PLEPH`,
- `jeulpkg.f`,
- the batch scripts,
- or the downstream filtering executables.

Only these build steps changed:

- recompile `jjreadnp.o`,
- compile `smu_probe_adjust_states.o`,
- relink `llr_npt_crd`.

## Probe Configuration

The gate used one bundled public sample:

- `s25y11d016t0231#103`

and two symmetric synthetic perturbations:

- `SMU_STATE_EM_SHIFT_M = -0.01 m`
- `SMU_STATE_EM_SHIFT_M = +0.01 m`

This is not a physical `delta_SEP` waveform.
It is a bounded state-level seam test.

## Outputs

- summary:
  [request4_llr_mlrs_state_seam_gate_summary.json](request4_llr_mlrs_state_seam_gate_summary.json#L1)
- figure:
  [request4_llr_mlrs_state_seam_gate_summary.svg](request4_llr_mlrs_state_seam_gate_summary.svg)
- runner:
  [request4_llr_mlrs_state_seam_gate.py](request4_llr_mlrs_state_seam_gate.py#L1)

## Result

The deeper state seam stays bounded.

For both `+/- 0.01 m` runs:

- the pipeline still runs end-to-end,
- the same `1456 / 1761` `93` residual lines respond,
- `frfin` changes only `12` lines,
- `frd` remains an exact match,
- and `npt` changes only one `11` record.

Numerically:

- `-0.01 m` gives mean `Δ93 = +6.671497e-2 ns`, `ΔPmM = +66.7 ps`
- `+0.01 m` gives mean `Δ93 = -6.671909e-2 ns`, `ΔPmM = -66.7 ps`

The odd-symmetry check is also clean enough for a gate decision:

- `frrec` mean-sum is `-4.12e-6 ns`
- `npt` time-of-flight sum is `0`
- `npt PmM` sum is `2.27e-13 ps`

So the synthetic state perturbation survives the `jjreadnp` bridge without
forcing edits into `PLEPH`, `jeulpkg.f`, or the batch stack.

## Interpretation

This stage answers the stop-rule question.

At the **software-interface** level:

- moving from the recalc seam to the deeper prediction/state seam does
  **not** widen into broad surgery.

At the **physics-model** level:

- nothing here proves that the real weak-field `delta_SEP` effect is correctly
  represented by this synthetic state displacement,
- and nothing here proves that a final weak-field posterior can be obtained
  without a stronger dynamical model.

That distinction should stay explicit.

The correct reading is:

- `MLRS` remains a live hand-off candidate,
- the broad-surgery stop rule is **not** triggered by the `jjreadnp` bridge
  itself,
- but the next escalation must still be judged on physics adequacy.

If a real `delta_SEP` insertion later demands:

- editing `PLEPH`,
- regenerating the ephemeris externally,
- or widening the change set far beyond the current local bridge,

then the branch should still stop and fall through to a mature external
estimator.
